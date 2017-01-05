#!/usr/bin/env python
"""Downloads and parses various Entrez databases."""
import os
import sys
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from urllib.error import HTTPError
from http.client import IncompleteRead
from xml.etree import ElementTree
from logging import INFO, DEBUG
import pandas as pd
import time
from shutil import rmtree

from Bio import Entrez
from mongoengine import connect

from sramongo.logger import logger
from sramongo.mongo import MongoDB
from sramongo.sra import SraExperiment, XMLSchemaException
from sramongo.mongo_schema import Submission, Organization,  Study, \
    Sample, Experiment, Run

STUDIES_ADDED = set()
EXPERIMENTS_ADDED = set()
RUNS_ADDED = set()

_DEBUG = False

class Cache(object):
    def __init__(self, directory='./'):
        cdir = os.path.join(directory, '.cache/sra2mongo')
        if not os.path.exists(cdir):
            os.makedirs(cdir)
        self.cachedir = os.path.abspath(cdir)
        self.cached = set([x.replace('.xml', '') for x in os.listdir(path=self.cachedir) if x.endswith('.xml')])

    def add_to_cache(self, name, _type, data):
        self.cached.add(name)

        if _type == 'xml':
            ext = 'xml'
        else:
            ext = 'csv'

        fname = '{}.{}'.format(name, ext)

        with open(os.path.join(self.cachedir, fname), 'w') as fh:
            fh.write(data)

    def get_cache(self, name, _type):
        if _type == 'xml':
            ext = 'xml'
        else:
            ext = 'csv'

        fname = os.path.join(self.cachedir, '{}.{}'.format(name, ext))

        if os.path.exists(fname):
            with open(fname, 'r') as fh:
                return fh.read()
        else:
            return None

    def __iter__(self):
        for f in sorted(self.cached):
            xml = os.path.join(self.cachedir, '{}.{}'.format(f, 'xml'))
            csv = os.path.join(self.cachedir, '{}.{}'.format(f, 'csv'))
            yield xml, csv

    def clear(self):
        rmtree(self.cachedir)
        self.cached = set()


def arguments():
    """Pulls in command line arguments."""

    DESCRIPTION = """\
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=Raw)

    parser.add_argument("--email", dest="email", action='store', required=True,
                        help="An email address is required for querying Entrez databases.")

    parser.add_argument("--query", dest="query", action='store', required=True,
                        help="Query to submit to Entrez.")

    parser.add_argument("--dbDir", dest="dbDir", action='store', required=True,
                        help="Folder containing mongo database.")

    parser.add_argument("--logDir", dest="logDir", action='store', required=False, default=None,
                        help="Folder containing mongo database.")

    parser.add_argument("--port", dest="port", action='store', type=int, required=False, default=27017,
                        help="Mongo database port.")

    parser.add_argument("--db", dest="db", action='store', required=False, default='sra',
                        help="Name of the database.")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False,
                        help="Turn on debug output.")

    parser.add_argument("--force", dest="force", action='store_true', required=False,
                        help="Forces clearing the cache.")

    return parser.parse_args(['--dbDir', '/data/db', '--logDir', '/data/log', '--email', 'justin.fear@nih.gov',
                              '--query', '"Drosophila melanogaster"[Orgn]', '--debug'])


def query_sra(query, **kwargs):
    options = {'term': query, 'db': 'sra', 'usehistory': 'y'}
    options.update(kwargs)
    handle = Entrez.esearch(**options)
    records = Entrez.read(handle)
    handle.close()

    try:
        records['Count'] = int(records['Count'])
        logger.info('There are {:,} records.'.format(records['Count']))
    except:
        pass

    return records


def fetch_sra(records, cache, runinfo_retmode='text', **kwargs):
    count = records['Count']
    batch_size = 500

    options = {'db': 'sra', 'webenv': records['WebEnv'],
               'query_key': records['QueryKey'], 'retmax': batch_size}
    options.update(kwargs)

    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)

        cache_xml = cache.get_cache(start, 'xml')
        cache_runinfo = cache.get_cache(start, 'runinfo')

        if (cache_xml is not None) & (cache_runinfo is not None):
            logger.info('Using cached version for: {} to {}.'.format(start+1, end))
        else:
            logger.info("Downloading record %i to %i" % (start+1, end))
            attempt = 0
            while attempt < 3:
                attempt += 1
                try:
                    xml_handle = Entrez.efetch(retstart=start, **options)
                    xml_data = xml_handle.read()
                    cache.add_to_cache(start, 'xml', xml_data)
                    xml_handle.close()

                    runinfo_handle = Entrez.efetch(rettype="runinfo", retmode=runinfo_retmode,
                                                   retstart=start, **options)
                    runinfo_data = runinfo_handle.read()
                    cache.add_to_cache(start, 'runinfo', runinfo_data)
                    runinfo_handle.close()

                except HTTPError as err:
                    if (500 <= err.code <= 599) & (attempt < 3):
                        logger.warn("Received error from server %s" % err)
                        logger.warn("Attempt %i of 3" % attempt)
                        time.sleep(15)
                    else:
                        logger.error("Received error from server %s" % err)
                        logger.error("Please re-run command latter.")
                        sys.exit(1)
                except IncompleteRead as err:
                    if attempt < 3:
                        logger.warn("Received error from server %s" % err)
                        logger.warn("Attempt %i of 3" % attempt)
                        time.sleep(15)
                    else:
                        logger.error("Received error from server %s" % err)
                        logger.error("Please re-run command latter.")
                        sys.exit(1)


def add_pkg_to_database(pkg, runinfo):

    sraExperiment = SraExperiment(pkg)

    # Add Experiment Package to database
    # Build study document
    submission = Submission.build_from_SraExperiment(sraExperiment)
    organization = Organization.build_from_SraExperiment(sraExperiment)
    study = Study.build_from_SraExperiment(sraExperiment, submission=submission,
                                           organization=organization)
    # Build sample document
    sample = Sample.build_from_SraExperiment(sraExperiment)

    # Build experiment document
    experiment = Experiment.build_from_SraExperiment(
        sraExperiment, study=study)

    # Add experiment id to study
    if study is not None:
        study.modify(push__experiments=experiment.experiment_id)
        STUDIES_ADDED.add(study.study_id)

    # Build run documents
    runs = Run.build_from_SraExperiment(sraExperiment, runinfo)
    RUNS_ADDED.update(set([run.run_id for run in runs]))

    if experiment is not None:
        experiment.modify(push_all__runs=[run.run_id for run in runs])
        EXPERIMENTS_ADDED.add(experiment.study_id)


def main():
    # Import commandline arguments.
    args = arguments()

    # Set logging level
    if args.debug:
        logger.setLevel(DEBUG)
        global _DEBUG
        _DEBUG = True
    else:
        logger.setLevel(INFO)

    cache = Cache()
    if args.force:
        logger.info('Clearing cache.')
        cache.clear()

    # Set email for entrez so they can contact if abbuse. Required by
    # biopython.
    Entrez.email = args.email

    logger.info('Starting MongoDB')
    with MongoDB(dbDir=args.dbDir, logDir=args.logDir, port=args.port):
        logger.info('Connecting to: {}'.format(args.db))
        connect(args.db)

        logger.info('Querying SRA: {}'.format(args.query))
        sra_query = query_sra(args.query)

        logger.info('Downloading documents')
        logger.info('Saving to cache: {}'.format(cache.cachedir))
        fetch_sra(sra_query, cache)

        logger.info('Adding documents to database')
        for xml, runinfo in cache:
            logger.debug('Parsing: {}'.format(xml))
            tree = ElementTree.parse(xml)
            ri = pd.read_csv(runinfo, index_col='Run')
            for exp_pkg in tree.findall('EXPERIMENT_PACKAGE'):
                add_pkg_to_database(exp_pkg, ri)

    logger.info('Studies Added: {:,}'.format(len(STUDIES_ADDED)))
    logger.info(STUDIES_ADDED)

    logger.info('Experiments Added: {:,}'.format(len(EXPERIMENTS_ADDED)))
    logger.info(EXPERIMENTS_ADDED)

    logger.info('Runs Added: {:,}'.format(len(RUNS_ADDED)))
    logger.info(RUNS_ADDED)
