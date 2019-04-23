#!/usr/bin/env python
"""Downloads and parses various Entrez databases."""
import sys
import time
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from logging import INFO, DEBUG
from datetime import datetime

import pymongo
from mongoengine import connect

from sramongo import parsers_sra_xml, parsers_bioproject_xml, parsers_biosample_xml, parsers_pubmed_xml
from sramongo.xml_helpers import xml_to_root
from sramongo.logger import logger
from sramongo.services import entrez
from sramongo.utils import chunked

_DEBUG = False
DEBUG_SIZE = 2_000
BATCH_SIZE = 200


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
                        help="Location of an already running database.")

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


def run_sra2mongo(args, collection):
    defaults = dict(api_key=args.api_key, email=args.email)

    if _DEBUG:
        defaults['retmax'] = DEBUG_SIZE

    ids_to_update = check_sra_for_updated_ids(args.query, collection, defaults)
    docs = download_sra_xml(ids_to_update, defaults)
    update_sramongo_sra_records(docs, collection)

    ids = get_bioproject_ids(collection)
    docs = download_bioproject_xml(ids, defaults)
    update_sramongo_bioproject_records(docs, collection)

    ids = get_biosample_ids(collection)
    docs = download_biosample_xml(ids, defaults)
    update_sramongo_biosample_records(docs, collection)

    ids = get_pubmed_ids(collection)
    docs = download_pubmed_xml(ids, defaults)
    update_sramongo_pubmed_records(docs, collection)


def check_sra_for_updated_ids(query, collection, defaults):
    logger.info(f'SRA - Querying SraMongo for last update')
    last_srx_update = get_sramongo_last_srx_update(collection)

    logger.info(f'SRA - Querying SRA for: {query}')
    esearch_result = entrez.esearch('sra', query, **defaults)
    webenv = esearch_result.webenv
    query_key = esearch_result.query_key

    if _DEBUG:
        count = DEBUG_SIZE
    else:
        count = esearch_result.count

    logger.info('SRA - Checking for Updates')
    ids_to_update = []
    for esummary_result in sra_esummary(webenv, query_key, count, defaults):
        srx = esummary_result.accn
        if esummary_result.update_date != last_srx_update.get(srx, None):
            ids_to_update.append(esummary_result.id)

    logger.info(f'SRA - {len(ids_to_update):,} IDs to Update')
    return ids_to_update


def download_sra_xml(ids_to_update, defaults):
    logger.info(f'SRA - Downloading XML Records')
    for i, ids in enumerate(chunked(ids_to_update, BATCH_SIZE)):
        start, end = i * BATCH_SIZE, (i + 1) * BATCH_SIZE
        logger.info(f'SRA - Downloading XML Records [Batch {start:,}-{end:,}]')

        epost_result = entrez.epost('sra', ids=ids, **defaults)
        webenv = epost_result.webenv
        query_key = epost_result.query_key
        count = len(ids)

        esummary_results = {
            result.accn: result
            for result in sra_esummary(webenv, query_key, count, defaults)
        }

        for srx, xml in sra_efetch(webenv, query_key, count, defaults):
            root = xml_to_root(xml)
            doc = parsers_sra_xml.parse_sra_experiment(root)
            doc.sra_id = int(esummary_results[srx].id)
            doc.sra_create_date = esummary_results[srx].create_date
            doc.sra_update_date = esummary_results[srx].update_date
            yield doc
        time.sleep(1)


def update_sramongo_sra_records(docs, collection):
    db_operations = []
    for doc in docs:
        db_operations.append(pymongo.ReplaceOne({'srx': doc.srx}, doc.to_mongo(), upsert=True))

        # Write intermediate results
        if len(db_operations) > 500:
            res = collection.bulk_write(db_operations)
            logger.debug(res.bulk_api_result)
            db_operations = []

    if len(db_operations) > 0:
        res = collection.bulk_write(db_operations)
        logger.debug(res.bulk_api_result)


def sra_esummary(webenv, query_key, count, defaults):
    for docs in entrez.esummary('sra', webenv=webenv, query_key=query_key, count=count, **defaults):
        yield from parsers_sra_xml.parse_sra_esummary_result(docs)


def sra_efetch(webenv, query_key, count, defaults):
    for docs in entrez.efetch('sra', webenv=webenv, query_key=query_key, count=count, **defaults):
        yield from parsers_sra_xml.parse_sra_efetch_result(docs)


def get_sramongo_last_srx_update(collection):
    return {
        record['srx']: record['sra_update_date']
        for record in collection.find({}, {'srx': True, 'sra_update_date': True})
    }


def get_bioproject_ids(collection):
    logger.info('BioProject - Getting IDs from SraMongo')
    bioproject_ids = set()
    for record in collection.find({'study.bioproject': {'$exists': True}, 'BioProject': {'$exists': False}},
                                  {'study.bioproject': True}):
        bioproject_ids.add(record['study']['bioproject'])

    # If BioProject is already there, don't update if it was added today.
    now = datetime.now()
    for record in collection.find({'BioProject': {'$exists': True}}, {'BioProject': True}):
        accn = record['BioProject']['accn']
        dt = record['BioProject']['sramongo_last_updated']
        if (now - dt).days > 1:
            bioproject_ids.add(accn)

    logger.info(f'BioProject - {len(bioproject_ids):,} IDs.')
    return bioproject_ids


def download_bioproject_xml(bioproject_ids, defaults):
    logger.info(f'BioProject - Downloading XML Records')
    for i, ids in enumerate(chunked(bioproject_ids, BATCH_SIZE)):
        start, end = i * BATCH_SIZE, (i + 1) * BATCH_SIZE
        logger.info(f'BioProject - Downloading XML Records [Batch {start:,}-{end:,}]')
        query = '+OR+'.join(ids)
        esearch_result = entrez.esearch('bioproject', query, **defaults)
        webenv = esearch_result.webenv
        query_key = esearch_result.query_key
        count = esearch_result.count
        for accn, xml in bioproject_efetch(webenv, query_key, count, defaults):
            root = xml_to_root(xml)
            doc = parsers_bioproject_xml.parse_bioproject(root)
            yield doc
        time.sleep(1)


def bioproject_efetch(webenv, query_key, count, defaults):
    for result in entrez.efetch('bioproject', webenv=webenv, query_key=query_key, count=count, **defaults):
        yield from parsers_bioproject_xml.parse_bioproject_efetch_result(result)


def update_sramongo_bioproject_records(docs, collection):
    db_operations = []
    for doc in docs:
        db_operations.append(
            pymongo.UpdateMany(
                {
                    'study.bioproject': doc.accn,
                    '$or': [
                        {'BioProject': {'$exists': False}},
                        {'BioProject.last_update': {'$ne': doc.last_update}},
                    ]
                },
                {'$set': {'BioProject': doc.to_mongo()}}
            )
        )

        # Write intermediate results
        if len(db_operations) > 500:
            res = collection.bulk_write(db_operations)
            logger.debug(res.bulk_api_result)
            db_operations = []

    if len(db_operations) > 0:
        res = collection.bulk_write(db_operations)
        logger.debug(res.bulk_api_result)


def get_biosample_ids(collection):
    logger.info('BioSample - Getting IDs from SraMongo')
    ids = set()
    for record in collection.find({'sample.biosample': {'$exists': True}, 'BioSample': {'$exists': False}},
                                  {'sample.biosample': True}):
        ids.add(record['sample']['biosample'])

    # If BioSample is already there, don't update if it was added today.
    now = datetime.now()
    for record in collection.find({'BioSample': {'$exists': True}}, {'BioSample': True}):
        accn = record['BioSample']['accn']
        dt = record['BioSample']['sramongo_last_updated']
        if (now - dt).days > 1:
            ids.add(accn)
    logger.info(f'BioSample - {len(ids):,} IDs.')
    return ids


def download_biosample_xml(biosample_ids, defaults):
    logger.info(f'BioSample - Downloading XML Records')
    for i, ids in enumerate(chunked(biosample_ids, BATCH_SIZE)):
        start, end = i * BATCH_SIZE, (i + 1) * BATCH_SIZE
        logger.info(f'BioSample - Downloading XML Records [Batch {start:,}-{end:,}]')
        query = '+OR+'.join(ids)
        esearch_result = entrez.esearch('biosample', query, **defaults)
        webenv = esearch_result.webenv
        query_key = esearch_result.query_key
        count = esearch_result.count
        for accn, xml in biosample_efetch(webenv, query_key, count, defaults):
            root = xml_to_root(xml)
            doc = parsers_biosample_xml.parse_biosample(root)
            yield doc
        time.sleep(1)


def biosample_efetch(webenv, query_key, count, defaults):
    for result in entrez.efetch('biosample', webenv=webenv, query_key=query_key, count=count, **defaults):
        yield from parsers_biosample_xml.parse_biosample_efetch_result(result)


def update_sramongo_biosample_records(docs, collection):
    db_operations = []
    for doc in docs:
        db_operations.append(
            pymongo.UpdateMany(
                {
                    'sample.biosample': doc.accn,
                    '$or': [
                        {'BioSample.accn': {'$exists': False}},
                        {'BioSample.last_update': {'$ne': doc.last_update}},
                    ]
                },
                {'$set': {'BioSample': doc.to_mongo()}}
            )
        )

        # Write intermediate results
        if len(db_operations) > 500:
            res = collection.bulk_write(db_operations)
            logger.debug(res.bulk_api_result)
            db_operations = []

    if len(db_operations) > 0:
        res = collection.bulk_write(db_operations)
        logger.debug(res.bulk_api_result)


def get_pubmed_ids(collection):
    logger.info('Pubmed - Getting IDs from SraMongo')
    ids = set()
    for record in collection.find({'study.pubmed': {'$exists': True}, 'papers': {'$eq': []}},
                                  {'study.pubmed': True}):
        ids |= set(record['study']['pubmed'])

    # If Pubmed is already there, don't update if it was added today.
    now = datetime.now()
    for record in collection.find({'papers': {'$ne': []}}, {'papers': True}):
        for rec in record['papers']:
            accn = rec['accn']
            dt = rec['sramongo_last_updated']
            if (now - dt).days > 1:
                ids.add(accn)

    logger.info(f'Pubmed - {len(ids):,} IDs.')
    return [str(_id) for _id in ids]


def download_pubmed_xml(pubmed_ids, defaults):
    logger.info(f'Pubmed - Downloading XML Records')
    for i, ids in enumerate(chunked(pubmed_ids, BATCH_SIZE)):
        start, end = i * BATCH_SIZE, (i + 1) * BATCH_SIZE
        logger.info(f'Pubmed - Downloading XML Records [Batch {start:,}-{end:,}]')
        epost_result = entrez.epost('pubmed', ids=ids, **defaults)
        webenv = epost_result.webenv
        query_key = epost_result.query_key
        count = len(ids)
        for accn, xml in pubmed_efetch(webenv, query_key, count, defaults):
            root = xml_to_root(xml)
            doc = parsers_pubmed_xml.parse_pubmed(root)
            yield doc
        time.sleep(1)


def pubmed_efetch(webenv, query_key, count, defaults):
    for result in entrez.efetch('pubmed', webenv=webenv, query_key=query_key, count=count, **defaults):
        yield from parsers_pubmed_xml.parse_pubmed_efetch_result(result)


def update_sramongo_pubmed_records(docs, collection):
    db_operations = []
    for doc in docs:
        db_operations.append(
            pymongo.UpdateMany(
                {
                    'study.pubmed': doc.accn,
                    '$or': [
                        {'papers': {'$eq': []}},
                        {'papers.date_revised': {'$elemMatch': {'$ne': doc.date_revised}}}
                    ]
                },
                {'$addToSet': {'papers': doc.to_mongo()}}
            )
        )

        # Write intermediate results
        if len(db_operations) > 500:
            res = collection.bulk_write(db_operations)
            logger.debug(res.bulk_api_result)
            db_operations = []

    if len(db_operations) > 0:
        res = collection.bulk_write(db_operations)
        logger.debug(res.bulk_api_result)


def main():
    args = arguments()

    # Set logging levelA
    if args.debug:
        global _DEBUG
        _DEBUG = True
        logger.setLevel(DEBUG)
        logger.debug(args)
    else:
        logger.setLevel(INFO)

    try:
        client = connect(args.db, host=args.host, port=args.port)
        db = client[args.db]
        collection: pymongo.mongo_client.database.Collection = db['ncbi']
        logger.info('Connecting to database: %s', str(collection))
        run_sra2mongo(args, collection)
    finally:
        logger.info('Closing Database Connection')
        client.close()


if __name__ == '__main__':
    main()
