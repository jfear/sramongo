#!/usr/bin/env python
"""Downloads and parses various Entrez databases."""
import os
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from urllib.error import HTTPError
from http.client import IncompleteRead
from xml.etree import ElementTree
from logging import INFO, DEBUG
import pandas as pd
from io import StringIO
from tempfile import NamedTemporaryFile
import time

from Bio import Entrez
from mongoengine import connect

from sramongo.logger import logger
from sramongo.mongo import MongoDB
from sramongo.sra import SraExperiment, XMLSchemaException
from sramongo.mongo_schema import Submission, Organization,  Study, \
    Sample, Pool, Experiment, Run

STUDIES_ADDED = set()
EXPERIMENTS_ADDED = set()
RUNS_ADDED = set()

_DEBUG = False

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

    return parser.parse_args(['--dbDir', '/data/db', '--logDir', '/data/log', '--email', 'justin.fear@nih.gov',
                              '--query', '"Drosophila melanogaster"[Orgn]', '--debug'])


class Cache(object):
    def __init__(self):
        if not os.path.exists('.cache/sra2mongo'):
            os.makedirs('.cache/sra2mongo')
        self.cachedir = os.path.abspath('.cache/sra2mongo')
        self.cached = set([x for x in os.listdir(path=self.cachedir) if x.endswith('.txt')])

    def add_to_cache(self, name, _type, data):
        fname = '{}_{}.txt'.format(name, _type)
        with open(os.path.join(self.cachedir, fname), 'w') as fh:
            fh.write(data)
        self.cached.add(fname)

    def get_cache(self, name, _type):
        fname = '{}_{}.txt'.format(name, _type)
        if fname in self.cached:
            with open(os.path.join(self.cachedir, fname), 'r') as fh:
                return fh.read()
        else:
            return None


def fetch_sra(records, runinfo_retmode='text', **kwargs):
    count = int(records['Count'])
    batch_size = 500

    options = {'db': 'sra', 'webenv': records['WebEnv'],
               'query_key': records['QueryKey'], 'retmax': batch_size}
    options.update(kwargs)

    cache = Cache()
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)

        cache_xml = cache.get_cache(start, 'xml')
        cache_runinfo = cache.get_cache(start, 'runinfo')

        if (cache_xml is not None) and (cache_runinfo is not None):
            logger.warn('Using cached for {} to {} instead of downloading. '
                        'To force download please remove: {}'.format(start+1, end, cache.cachedir))
            xml_data = cache_xml
            runinfo_data = cache_runinfo
        else:
            logger.info("Downloading record %i to %i" % (start+1, end))
            attempt = 0
            while attempt < 3:
                attempt += 1
                try:
                    xml_handle = Entrez.efetch(retstart=start, **options)

                    runinfo_handle = Entrez.efetch(rettype="runinfo", retmode=runinfo_retmode,
                                                   retstart=start, **options)

                    xml_data = xml_handle.read()
                    xml_handle.close()
                    cache.add_to_cache(start, 'xml', xml_data)

                    runinfo_data = runinfo_handle.read()
                    runinfo_handle.close()
                    cache.add_to_cache(start, 'runinfo', runinfo_data)

                except HTTPError as err:
                    if 500 <= err.code <= 599:
                        logger.warn("Received error from server %s" % err)
                        logger.warn("Attempt %i of 3" % attempt)
                        time.sleep(15)
                    else:
                        raise
                except IncompleteRead:
                    if attempt < 3:
                        logger.warn('Incomplete read, trying again. %s' % start)
                        logger.warn("Attempt %i of 3" % attempt)
                        time.sleep(15)
                    else:
                        raise

        yield xml_data, runinfo_data


def query_sra(query, **kwargs):
    options = {'term': query, 'db': 'sra', 'usehistory': 'y'}
    options.update(kwargs)
    handle = Entrez.esearch(**options)
    records = Entrez.read(handle)
    handle.close()
    return records


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
    pool = Pool.build_from_SraExperiment(sraExperiment)
    experiment = Experiment.build_from_SraExperiment(sraExperiment,
                                                     study=study, samples=pool)

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

    # Set email for entrez so they can contact if abbuse. Required by
    # biopython.
    Entrez.email = args.email

    logger.info('Starting MongoDB')
    with MongoDB(dbDir=args.dbDir, logDir=args.logDir, port=args.port):
        client = connect(args.db)
        try:
            logger.info('Querying SRA: {}'.format(args.query))
            sra_query = query_sra(args.query)

            logger.info('Downloading documents')
            for xml, runinfo in fetch_sra(sra_query):
                tree = ElementTree.fromstring(xml)
                ri = pd.read_csv(StringIO(runinfo), index_col='Run')
                for exp_pkg in tree:
                    add_pkg_to_database(exp_pkg, ri)
        except:
            raise
        finally:
           client.close()

           logger.info('Studies Added')
           logger.info(STUDIES_ADDED)

           logger.info('Experiments Added')
           logger.info(EXPERIMENTS_ADDED)

           logger.info('Runs Added')
           logger.info(RUNS_ADDED)

