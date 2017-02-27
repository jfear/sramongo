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
from mongoengine.errors import ValidationError

from sramongo.logger import logger
from sramongo.mongo import MongoDB
from sramongo.sra import SraExperiment, XMLSchemaException
from sramongo.biosample import BioSampleParse
from sramongo.mongo_schema import Ncbi

from ipdb import slaunch_ipdb_on_exception

_DEBUG = False

class Cache(object):
    def __init__(self, directory='', clean=False):
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.cachedir = os.path.abspath(directory)

        if clean:
            logger.info('Clearing cache: {}.'.format(cache.cachedir))
            self.clear()

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
            if os.path.exists(csv):
                yield xml, csv
            else:
                yield xml

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
                        help="Folder in which to store log.")

    parser.add_argument("--port", dest="port", action='store', type=int, required=False, default=27017,
                        help="Mongo database port.")

    parser.add_argument("--db", dest="db", action='store', required=False, default='sra',
                        help="Name of the database.")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False,
                        help="Turn on debug output.")

    parser.add_argument("--force", dest="force", action='store_true', required=False,
                        help="Forces clearing the cache.")

    return parser.parse_args()


def iter_query(query):
    batch_size = 5000
    for start in range(0, len(query), batch_size):
        end = min(len(query), start+batch_size)
        yield '[accn] OR '.join(query[start:end]) + '[accn]'


def ncbi_query(query, **kwargs):
    options = {'term': query, 'db': 'sra'}
    options.update(kwargs)

    if isinstance(query, str):
        # Simple query
        options['usehistory'] = 'y'
        handle = Entrez.esearch(**options)
        records = Entrez.read(handle)
        records['Count'] = int(records['Count'])
        logger.info('There are {:,} records.'.format(records['Count']))
        return records
    elif isinstance(query, list):
        logger.debug('Using list of accessions.')
        # This is used when there is a list of accession numbers to query.
        options['retmax'] = 100000          # the highest valid value

        # Some of the biosample_accn are actually biosamples IDs lets pull
        # those out.
        ids = []
        accn = []
        for _id in query:
            if _id.startswith('SAM'):
                accn.append(_id)
            else:
                ids.append(_id)

        for q in iter_query(accn):
            options['term'] = q
            handle = Entrez.esearch(**options)
            records = Entrez.read(handle)
            ids.extend(records['IdList'])
            handle.close()

        handle = Entrez.epost(options['db'], id=','.join(ids))
        records = Entrez.read(handle)
        handle.close()
        records['Count'] = len(ids)
        logger.info('There are {:,} records.'.format(records['Count']))
        return records
    else:
        logger.debug(query[:10])
        raise ValueError('Query should be a string or list of Accession numbers.')


def fetch_sra(records, cache, runinfo_retmode='text', **kwargs):
    count = records['Count']
    batch_size = 500

    options = {'db': 'sra', 'webenv': records['WebEnv'],
               'query_key': records['QueryKey'], 'retmax': batch_size}
    options.update(kwargs)

    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)

        cache_xml = cache.get_cache(start, 'xml')
        cache_runinfo = cache.get_cache(start, 'csv')

        if (cache_xml is not None):
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

                    if runinfo_retmode is not None:
                        runinfo_handle = Entrez.efetch(rettype="runinfo", retmode=runinfo_retmode,
                                                       retstart=start, **options)
                        runinfo_data = runinfo_handle.read()
                        cache.add_to_cache(start, 'csv', runinfo_data)
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


def parse_sra(cache):
    for xml, runinfo in cache:
        logger.debug('Parsing: {}'.format(xml))
        # Import XML
        tree = ElementTree.parse(xml)

        # Import runinfo
        runinfo = pd.read_csv(runinfo, index_col='Run')

        # Iterate over experiments and parse
        for exp_pkg in tree.findall('EXPERIMENT_PACKAGE'):
            sraExperiment = SraExperiment(exp_pkg)
            # Iterate over runs and update missing fields from runinfo
            for i, run in enumerate(sraExperiment.sra['run']):
                srr = run['run_id']
                try:
                    rinfo = runinfo.loc[srr].dropna().to_dict()
                    if 'ReleaseDate' in rinfo:
                        sraExperiment.sra['run'][i]['release_date'] = rinfo['ReleaseDate']
                    if 'LoadDate' in rinfo:
                        sraExperiment.sra['run'][i]['load_date'] = rinfo['LoadDate']
                    if 'size_MB' in rinfo:
                        sraExperiment.sra['run'][i]['size_MB'] = rinfo['size_MB']
                    if 'download_path' in rinfo:
                        sraExperiment.sra['run'][i]['download_path'] = rinfo['download_path']
                except KeyError:
                    pass

            with slaunch_ipdb_on_exception():
                Ncbi.objects(pk=sraExperiment.srx).modify(upsert=True, sra=sraExperiment.sra)


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

    # Connect to mongo server
    logger.info('Starting MongoDB')
    with MongoDB(dbDir=args.dbDir, logDir=args.logDir, port=args.port):
        logger.info('Connecting to: {}'.format(args.db))
        connect(args.db)

        # SRA
        logger.info('Querying SRA: {}'.format(args.query))
        sra_query = ncbi_query(args.query)

        logger.info('Downloading documents')
        cache = Cache(directory='.cache/sra2mongo/sra', clean=args.force)

        logger.info('Saving to cache: {}'.format(cache.cachedir))
        fetch_sra(sra_query, cache)

        logger.info('Parsing SRA XML')
        parse_sra(cache)

#         # Query BioSample
#         bs_cache = Cache(directory='.cache/sra2mongo/biosample')
#         biosample_accn = list(BioSample.objects.distinct('biosample_id'))
#
#         logger.info('Querying BioSample: with {:,} ids'.format(len(biosample_accn)))
#         logger.info('Saving to cache: {}'.format(bs_cache.cachedir))
#         bs_query = ncbi_query(biosample_accn, db='biosample')
#
#         logger.info('Downloading BioSample documents')
#         logger.info('Saving to cache: {}'.format(bs_cache.cachedir))
#         fetch_sra(bs_query, bs_cache, runinfo_retmode=None, db='biosample')
#
#         logger.info('Adding documents to database')
#         for xml in bs_cache:
#             logger.debug('Parsing: {}'.format(xml))
#             tree = ElementTree.parse(xml)
#             for pkg in tree.findall('BioSample'):
#                 bs = BioSampleParse(pkg)
#                 b = BioSample.build_from_BioSample(bs)


if __name__ == '__main__':
    main()
