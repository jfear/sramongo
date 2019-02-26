#!/usr/bin/env python
"""Downloads and parses various Entrez databases.

import os
import sys
from urllib.error import HTTPError, URLError
from http.client import IncompleteRead
from xml.etree import ElementTree
import pandas as pd
import time
from shutil import rmtree
from glob import glob
from io import StringIO

from Bio import Entrez

from sramongo.mongo import MongoDB2
from sramongo.sra import SraExperiment, XMLSchemaException
from sramongo.biosample import BioSampleParse
from sramongo.bioproject import BioProjectParse
from sramongo.pubmed import PubmedParse
from sramongo.mongo_schema import Ncbi
"""

import sys
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from logging import INFO, DEBUG

import pymongo
from mongoengine import connect
from mongoengine.errors import ValidationError

from sramongo import parsers_sra_xml, parsers_bioproject_xml, parsers_biosample_xml, parsers_pubmed_xml
from sramongo.xml_helpers import xml_to_root
from sramongo.logger import logger
from sramongo.services import entrez

_DEBUG = False



def arguments():
    """Pulls in command line arguments."""

    DESCRIPTION = """\
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=Raw)

    parser.add_argument("--email", dest="email", action='store', required=False, default=False,
                        help="An email address is required for querying Entrez databases.")

    parser.add_argument("--api", dest="api_key", action='store', required=False, default=False,
                        help="A users ENTREZ API Key. Will speed up download.")

    parser.add_argument("--query", dest="query", action='store', required=True,
                        help="Query to submit to Entrez.")

    parser.add_argument("--host", dest="host", action='store', required=False, default='localhost',
                        help="Location of an already running database. If provided ignores --dbpath and --logDir.")

    parser.add_argument("--port", dest="port", action='store', type=int, required=False, default=27017,
                        help="Mongo database port.")

    parser.add_argument("--db", dest="db", action='store', required=False, default='sramongo',
                        help="Name of the database.")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False,
                        help="Turn on debug output.")

    parser.add_argument("--force", dest="force", action='store_true', required=False,
                        help="Forces clearing the cache.")

    args = parser.parse_args()

    if not (args.email or args.api_key):
        logger.error('You must provide either an `--email` or `--api`.')
        sys.exit()

    return args


# def iter_query(query):
#     batch_size = 5000
#     for start in range(0, len(query), batch_size):
#         end = min(len(query), start+batch_size)
#         if query[0].startswith('SAM'):
#             yield '[accn] OR '.join(query[start:end]) + '[accn]'
#         elif query[0].startswith('PRJ'):
#             # For some reason you cannot query BioProject using the [accn]
#             # keyword.
#             yield ' OR '.join(query[start:end])
#
#
# def ncbi_query(query, **kwargs):
#     options = {'term': query, 'db': 'sra'}
#     options.update(kwargs)
#
#     if isinstance(query, str):
#         # Simple query
#         options['usehistory'] = 'y'
#         handle = Entrez.esearch(**options)
#         records = Entrez.read(handle)
#         records['Count'] = int(records['Count'])
#         logger.info('There are {:,} records.'.format(records['Count']))
#         return records
#     elif isinstance(query, list):
#         logger.debug('Using list of accessions.')
#         # This is used when there is a list of accession numbers to query.
#         options['retmax'] = 100000          # the highest valid value
#
#         # Some of accn are actually IDs lets pull those out.
#         ids = []
#         accn = []
#         for _id in set(query):
#             if _id.startswith('SAM') or _id.startswith('PRJ'):
#                 accn.append(_id)
#             elif _id.startswith('PMID'):
#                 ids.append(_id.replce('PMID', ''))
#             else:
#                 ids.append(_id)
#
#         for q in iter_query(accn):
#             options['term'] = q
#             handle = Entrez.esearch(**options)
#             records = Entrez.read(handle)
#             ids.extend(records['IdList'])
#             handle.close()
#
#         ids = list(set(ids))
#         handle = Entrez.epost(options['db'], id=','.join(ids))
#         records = Entrez.read(handle)
#         handle.close()
#         records['Count'] = len(ids)
#         logger.info('There are {:,} records.'.format(records['Count']))
#         return records
#     else:
#         logger.debug(query[:10])
#         raise ValueError('Query should be a string or list of Accession numbers.')
#
#
# def catch_xml_error(xml):
#     try:
#         tree = ElementTree.parse(xml)
#     except ElementTree.ParseError:
#         logger.debug('Current XML file appears to be empty; re-download.')
#         return True
#
#     res = tree.find('ERROR')
#     if res is None:
#         return False
#     else:
#         logger.debug('Current XML has an ERROR tag; re-download.')
#         return True
#
#
# def fetch_ncbi(records, cache, runinfo_retmode='text', **kwargs):
#     count = records['Count']
#     batch_size = 500
#
#     options = {'db': 'sra', 'webenv': records['WebEnv'],
#                'query_key': records['QueryKey'], 'retmax': batch_size}
#     options.update(kwargs)
#
#     for start in range(0, count, batch_size):
#         end = min(count, start+batch_size)
#
#         cache_xml = cache.get_cache(start, 'xml')
#         cache_runinfo = cache.get_cache(start, 'csv')
#
#         if (cache_xml is not None) & (catch_xml_error(StringIO(cache_xml)) is False):
#             logger.info('Using cached version for: {} to {}.'.format(start+1, end))
#         else:
#             logger.info("Downloading record %i to %i" % (start+1, end))
#             attempt = 0
#             while attempt < 3:
#                 attempt += 1
#                 try:
#                     xml_handle = Entrez.efetch(retmode='xml', retstart=start, **options)
#                     xml_data = xml_handle.read()
#                     cache.add_to_cache(start, 'xml', xml_data)
#                     xml_handle.close()
#
#                     if runinfo_retmode is not None:
#                         runinfo_handle = Entrez.efetch(rettype="runinfo", retmode=runinfo_retmode,
#                                                        retstart=start, **options)
#                         runinfo_data = runinfo_handle.read()
#                         cache.add_to_cache(start, 'csv', runinfo_data)
#                         runinfo_handle.close()
#
#                 except HTTPError as err:
#                     if (500 <= err.code <= 599) & (attempt < 3):
#                         logger.warning("Received error from server %s" % err)
#                         logger.warning("Attempt %i of 3" % attempt)
#                         time.sleep(15)
#                     else:
#                         logger.error("Received error from server %s" % err)
#                         logger.error("Please re-run command latter.")
#                         cache.remove_from_cache(start)
#                         sys.exit(1)
#                 except IncompleteRead as err:
#                     if attempt < 3:
#                         logger.warning("Received error from server %s" % err)
#                         logger.warning("Attempt %i of 3" % attempt)
#                         time.sleep(15)
#                     else:
#                         logger.error("Received error from server %s" % err)
#                         logger.error("Please re-run command latter.")
#                         cache.remove_from_cache(start)
#                         sys.exit(1)
#                 except URLError as err:
#                     if (err.code == 60) & (attempt < 3):
#                         logger.warning("Received error from server %s" % err)
#                         logger.warning("Attempt %i of 3" % attempt)
#                         cache.remove_from_cache(start)
#                         time.sleep(15)
#                     else:
#                         logger.error("Received error from server %s" % err)
#                         logger.error("Please re-run command latter.")
#                         cache.remove_from_cache(start)
#                         sys.exit(1)
#
#
# def parse_sra(cache):
#     for xml, runinfo in cache:
#         logger.debug('Parsing: {}'.format(xml))
#         # Import XML
#         tree = ElementTree.parse(xml)
#
#         # Import runinfo
#         runinfo = pd.read_csv(runinfo, index_col='Run')
#
#         # Iterate over experiments and parse
#         for exp_pkg in tree.findall('EXPERIMENT_PACKAGE'):
#             sraExperiment = SraExperiment(exp_pkg)
#             # Iterate over runs and update missing fields from runinfo
#             for i, run in enumerate(sraExperiment.sra['run']):
#                 srr = run['run_id']
#                 try:
#                     rinfo = runinfo.loc[srr].dropna().to_dict()
#                     if 'ReleaseDate' in rinfo:
#                         sraExperiment.sra['run'][i]['release_date'] = rinfo['ReleaseDate']
#                     if 'LoadDate' in rinfo:
#                         sraExperiment.sra['run'][i]['load_date'] = rinfo['LoadDate']
#                     if 'size_MB' in rinfo:
#                         sraExperiment.sra['run'][i]['size_MB'] = rinfo['size_MB']
#                     if 'download_path' in rinfo:
#                         sraExperiment.sra['run'][i]['download_path'] = rinfo['download_path']
#                 except KeyError:
#                     pass
#
#             Ncbi.objects(pk=sraExperiment.srx).modify(upsert=True, sra=sraExperiment.sra)
#
#
# def parse_biosample(cache):
#     for xml in cache:
#         logger.debug('Parsing: {}'.format(xml))
#         tree = ElementTree.parse(xml)
#         for pkg in tree.findall('BioSample'):
#             bs = BioSampleParse(pkg)
#             try:
#                 Ncbi.objects(sra__sample__BioSample=bs.biosample['biosample_accn']).update(add_to_set__biosample=bs.biosample)
#             except KeyError:
#                 pass
#
#             try:
#                 Ncbi.objects(sra__sample__BioSample=bs.biosample['biosample_id']).update(add_to_set__biosample=bs.biosample)
#             except KeyError:
#                 pass
#
#
# def parse_bioproject(cache):
#     for xml in cache:
#         logger.debug('Parsing: {}'.format(xml))
#         tree = ElementTree.parse(xml)
#         for pkg in tree.findall('DocumentSummary'):
#             bs = BioProjectParse(pkg)
#             Ncbi.objects(sra__study__BioProject=bs.bioproject['bioproject_accn']).update(bioproject=bs.bioproject)
#             Ncbi.objects(sra__study__BioProject=bs.bioproject['bioproject_id']).update(bioproject=bs.bioproject)
#
#
# def parse_pubmed(cache):
#     for xml in cache:
#         logger.debug('Parsing: {}'.format(xml))
#         tree = ElementTree.parse(xml)
#         for pkg in tree.findall('PubmedArticle/MedlineCitation'):
#             pm = PubmedParse(pkg)
#             Ncbi.objects(sra__study__pubmed=pm.pubmed['pubmed_id']).update(add_to_set__pubmed=pm.pubmed)
#             Ncbi.objects(sra__sample__pubmed=pm.pubmed['pubmed_id']).update(add_to_set__pubmed=pm.pubmed)
#             Ncbi.objects(sra__experiment__pubmed=pm.pubmed['pubmed_id']).update(add_to_set__pubmed=pm.pubmed)
#
#
# def main_old():
#     # Import commandline arguments.
#     args = arguments()
#
#     # Set logging level
#     if args.debug:
#         logger.setLevel(DEBUG)
#         global _DEBUG
#         _DEBUG = True
#     else:
#         logger.setLevel(INFO)
#
#     # Set email for entrez so they can contact if abbuse. Required by
#     # biopython.
#     Entrez.email = args.email
#
#     # Connect to mongo server
#     if args.host:
#         logger.info('Connecting to running server at: mongodb://{}:{}'.format(args.host, args.port))
#         host = args.host
#         mongodb = False
#     else:
#         logger.info('Starting MongoDB')
#         mongodb = MongoDB2(dbDir=args.dbDir, logDir=args.logDir, port=args.port)
#         host = '127.0.0.1'
#         logger.info('Running server at: mongodb://{}:{}'.format(host, args.port))
#
#     logger.info('Connecting to: {}'.format(args.db))
#     client = connect(args.db, host=host, port=args.port)
#
#     try:
#         # SRA
#         logger.info('Querying SRA: {}'.format(args.query))
#         sra_query = ncbi_query(args.query)
#
#         logger.info('Downloading documents')
#         cache = Cache(directory='.cache/sra2mongo/sra', clean=args.force)
#
#         logger.info('Saving to cache: {}'.format(cache.cachedir))
#         fetch_ncbi(sra_query, cache)
#
#         logger.info('Parsing SRA XML')
#         parse_sra(cache)
#
#         # Query BioSample
#         cache = Cache(directory='.cache/sra2mongo/biosample')
#         accn = list(Ncbi.objects.distinct('sra.sample.BioSample'))
#
#         logger.info('Querying BioSample: with {:,} ids'.format(len(accn)))
#         logger.info('Saving to cache: {}'.format(cache.cachedir))
#         query = ncbi_query(accn, db='biosample')
#
#         logger.info('Downloading BioSample documents')
#         logger.info('Saving to cache: {}'.format(cache.cachedir))
#         fetch_ncbi(query, cache, runinfo_retmode=None, db='biosample')
#
#         logger.info('Adding documents to database')
#         parse_biosample(cache)
#
#         # Query BioProject
#         cache = Cache(directory='.cache/sra2mongo/bioproject')
#         accn = list(Ncbi.objects.distinct('sra.study.BioProject'))
#
#         logger.info('Querying BioProject: with {:,} ids'.format(len(accn)))
#         logger.info('Saving to cache: {}'.format(cache.cachedir))
#         query = ncbi_query(accn, db='bioproject')
#
#         logger.info('Downloading BioProject documents')
#         logger.info('Saving to cache: {}'.format(cache.cachedir))
#         fetch_ncbi(query, cache, runinfo_retmode=None, db='bioproject')
#
#         logger.info('Adding documents to database')
#         parse_bioproject(cache)
#
#         # Query Pubmed
#         cache = Cache(directory='.cache/sra2mongo/pubmed')
#         accn = list(Ncbi.objects.distinct('sra.study.pubmed'))
#         accn.extend(list(Ncbi.objects.distinct('sra.experiment.pubmed')))
#         accn.extend(list(Ncbi.objects.distinct('sra.sample.pubmed')))
#         accn = list(set(accn))
#
#         logger.info('Querying Pubmed: with {:,} ids'.format(len(accn)))
#         logger.info('Saving to cache: {}'.format(cache.cachedir))
#         query = ncbi_query(accn, db='pubmed')
#
#         logger.info('Downloading Pubmed documents')
#         logger.info('Saving to cache: {}'.format(cache.cachedir))
#         fetch_ncbi(query, cache, runinfo_retmode=None, db='pubmed')
#
#         logger.info('Adding documents to database')
#         parse_pubmed(cache)
#     except Exception as err:
#         raise err
#     finally:
#         client.close()
#         if mongodb:
#             mongodb.close()


def get_ids_to_download(collection, esummary_results):
    """Check sramongo database and see if we need to download record."""
    srxs = list(esummary_results.keys())
    sramongo_results = {
        rec['srx']: rec['sra_update_date']
        for rec in collection.find({'srx': {'$in': srxs}},
                                   {'srx': True, 'sra_update_date': True})
    }

    for srx in srxs:
        sra_date = esummary_results[srx].update_date
        sramongo_date = sramongo_results.get(srx, None)
        if sra_date != sramongo_date:
            yield esummary_results[srx].id


def iter_sra_results(epost_result, **kwargs):
    for result in entrez.efetch('sra', webenv=epost_result.webenv, query_key=epost_result.query_key, **kwargs):
        yield from parsers_sra_xml.parse_sra_experiment_set(result)


def process_sra(collection, esummary_results, epost_result, **kwargs):
    for srx, experiment_xml in iter_sra_results(epost_result, **kwargs):
        root = xml_to_root(experiment_xml)
        doc = parsers_sra_xml.parse_sra_experiment(root)
        pass



def process_bioproject(collection, epost_result, **kwargs):
    links = entrez.elink(db='bioproject', dbfrom='sra', webenv=epost_result.webenv, query_key=epost_result.query_key,
                         **kwargs)
    if links.query_key:
        for result in entrez.efetch('bioproject', webenv=links.webenv, query_key=links.query_key, **defaults):
            yield from parsers_bioproject_xml.parse_bioproject_set(result)


def process_biosample(collection, epost_result, **kwargs):
    links = entrez.elink(db='biosample', dbfrom='sra', webenv=epost_result.webenv, query_key=epost_result.query_key,
                         **kwargs)
    if links.query_key:
        for result in entrez.efetch('biosample', webenv=links.webenv, query_key=links.query_key, **defaults):
            yield from parsers_biosample_xml.parse_biosample_set(result)


def process_pubmed(collection, epost_result, **kwargs):
    links = entrez.elink(db='pubmed', dbfrom='sra', webenv=epost_result.webenv, query_key=epost_result.query_key,
                         **kwargs)
    if links.query_key:
        for result in entrez.efetch('pubmed', webenv=links.webenv, query_key=links.query_key, **defaults):
            yield from parsers_pubmed_xml.parse_pubmed_set(result)


def main(args, collection):
    defaults = dict(api_key=args.api_key, email=args.email)

    if _DEBUG:
        defaults['retmax'] = 10

    # Search Entrez
    esearch_result = entrez.esearch('sra', args.query, **defaults)
    esummary_results = {
        result.srx: result
        for result in entrez.esummary('sra', webenv=esearch_result.webenv,
                                      query_key=esearch_result.query_key, **defaults)
    }

    # Figure out what records we need to download
    ids_to_download = list(get_ids_to_download(collection, esummary_results))
    epost_result = entrez.epost('sra', ids=ids_to_download, **defaults)
    logger.debug(epost_result)

    process_sra(collection, esummary_results, epost_result, **defaults)
    # process_bioproject(epost_result, **defaults)
    # biosample_results = process_biosample(epost_result, **defaults)
    # pubmed_results = process_pubmed(epost_result, **defaults)


if __name__ == '__main__':
    args = arguments()

    # Set logging levelA
    if args.debug:
        _DEBUG = True
        logger.setLevel(DEBUG)
        logger.debug(args)
    else:
        logger.setLevel(INFO)

    try:
        client = connect(args.db, host=args.host, port=args.port)
        db = client[args.db]
        collection: pymongo.mongo_client.database.Collection = db['sra']
        main(args, collection)
    finally:
        client.close()
