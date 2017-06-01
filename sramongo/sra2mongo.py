#!/usr/bin/env python
"""Downloads and parses various Entrez databases."""
import os
import sys
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from urllib.error import HTTPError, URLError
from http.client import IncompleteRead
from xml.etree import ElementTree
from logging import INFO, DEBUG
import pandas as pd
import time
from shutil import rmtree
from glob import glob
from io import StringIO

from Bio import Entrez
from mongoengine import connect
from mongoengine.errors import ValidationError

from sramongo.logger import logger
from sramongo.mongo import MongoDB2
from sramongo.sra import SraExperiment, XMLSchemaException
from sramongo.biosample import BioSampleParse
from sramongo.bioproject import BioProjectParse
from sramongo.pubmed import PubmedParse
from sramongo.mongo_schema import Ncbi

_DEBUG = False

class Cache(object):
    def __init__(self, directory='', clean=False):

        # Clean cache if option given
        if clean & os.path.exists(directory):
            logger.info('Clearing cache: {}.'.format(directory))
            self.clear()

        # Make cach directory
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.cachedir = os.path.abspath(directory)

        # Scan cache and make list of names.
        self.cached = set([x.replace('.xml', '') for x in os.listdir(path=self.cachedir) if x.endswith('.xml')])

    def add_to_cache(self, name, _type, data):
        """Add downloaded data to the cache."""
        self.cached.add(str(name))

        if _type == 'xml':
            ext = 'xml'
        else:
            ext = 'csv'

        fname = '{}.{}'.format(name, ext)

        with open(os.path.join(self.cachedir, fname), 'w') as fh:
            fh.write(data)

    def remove_from_cache(self, name):
        """Remove file from cache."""
        self.cached.discard(str(name))
        for fn in glob(os.path.join(self.cachedir, str(name) + '*')):
            if os.path.exists(fn):
                os.remove(fn)

    def get_cache(self, name, _type):
        """Get file contents from cache."""
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

    parser.add_argument("--host", dest="host", action='store', required=False,
                        help="Location of an already running database. If provided ignores --dbpath and --logDir.")

    parser.add_argument("--dbpath", dest="dbDir", action='store', required=False,
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

    args = parser.parse_args()

    if args.host and (args.dbDir or args.logDir):
        logger.warning('Both --host and --dbpath/--logDir were provided; going to try using running mongo server.')

    if not (args.host or args.dbDir or args.logDir):
        logger.warning('You must provide either a `--host hostname` with a running database or '
                '`--dbpath <path to db> --logDir <path to log>`')

    return args


def iter_query(query):
    batch_size = 5000
    for start in range(0, len(query), batch_size):
        end = min(len(query), start+batch_size)
        if query[0].startswith('SAM'):
            yield '[accn] OR '.join(query[start:end]) + '[accn]'
        elif query[0].startswith('PRJ'):
            # For some reason you cannot query BioProject using the [accn]
            # keyword.
            yield ' OR '.join(query[start:end])


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

        # Some of accn are actually IDs lets pull those out.
        ids = []
        accn = []
        for _id in set(query):
            if _id.startswith('SAM') or _id.startswith('PRJ'):
                accn.append(_id)
            elif _id.startswith('PMID'):
                ids.append(_id.replce('PMID', ''))
            else:
                ids.append(_id)

        for q in iter_query(accn):
            options['term'] = q
            handle = Entrez.esearch(**options)
            records = Entrez.read(handle)
            ids.extend(records['IdList'])
            handle.close()

        ids = list(set(ids))
        handle = Entrez.epost(options['db'], id=','.join(ids))
        records = Entrez.read(handle)
        handle.close()
        records['Count'] = len(ids)
        logger.info('There are {:,} records.'.format(records['Count']))
        return records
    else:
        logger.debug(query[:10])
        raise ValueError('Query should be a string or list of Accession numbers.')


def catch_xml_error(xml):
    try:
        tree = ElementTree.parse(xml)
    except ElementTree.ParseError:
        logger.debug('Current XML file appears to be empty; re-download.')
        return True

    res = tree.find('ERROR')
    if res is None:
        return False
    else:
        logger.debug('Current XML has an ERROR tag; re-download.')
        return True


def fetch_ncbi(records, cache, runinfo_retmode='text', **kwargs):
    count = records['Count']
    batch_size = 500

    options = {'db': 'sra', 'webenv': records['WebEnv'],
               'query_key': records['QueryKey'], 'retmax': batch_size}
    options.update(kwargs)

    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)

        cache_xml = cache.get_cache(start, 'xml')
        cache_runinfo = cache.get_cache(start, 'csv')

        if (cache_xml is not None) & (catch_xml_error(StringIO(cache_xml)) is False):
            logger.info('Using cached version for: {} to {}.'.format(start+1, end))
        else:
            logger.info("Downloading record %i to %i" % (start+1, end))
            attempt = 0
            while attempt < 3:
                attempt += 1
                try:
                    xml_handle = Entrez.efetch(retmode='xml', retstart=start, **options)
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
                        logger.warning("Received error from server %s" % err)
                        logger.warning("Attempt %i of 3" % attempt)
                        time.sleep(15)
                    else:
                        logger.error("Received error from server %s" % err)
                        logger.error("Please re-run command latter.")
                        cache.remove_from_cache(start)
                        sys.exit(1)
                except IncompleteRead as err:
                    if attempt < 3:
                        logger.warning("Received error from server %s" % err)
                        logger.warning("Attempt %i of 3" % attempt)
                        time.sleep(15)
                    else:
                        logger.error("Received error from server %s" % err)
                        logger.error("Please re-run command latter.")
                        cache.remove_from_cache(start)
                        sys.exit(1)
                except URLError as err:
                    if (err.code == 60) & (attempt < 3):
                        logger.warning("Received error from server %s" % err)
                        logger.warning("Attempt %i of 3" % attempt)
                        cache.remove_from_cache(start)
                        time.sleep(15)
                    else:
                        logger.error("Received error from server %s" % err)
                        logger.error("Please re-run command latter.")
                        cache.remove_from_cache(start)
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

            Ncbi.objects(pk=sraExperiment.srx).modify(upsert=True, sra=sraExperiment.sra)


def parse_biosample(cache):
    for xml in cache:
        logger.debug('Parsing: {}'.format(xml))
        tree = ElementTree.parse(xml)
        for pkg in tree.findall('BioSample'):
            bs = BioSampleParse(pkg)
            try:
                Ncbi.objects(sra__sample__BioSample=bs.biosample['biosample_accn']).update(add_to_set__biosample=bs.biosample)
            except KeyError:
                pass

            try:
                Ncbi.objects(sra__sample__BioSample=bs.biosample['biosample_id']).update(add_to_set__biosample=bs.biosample)
            except KeyError:
                pass


def parse_bioproject(cache):
    for xml in cache:
        logger.debug('Parsing: {}'.format(xml))
        tree = ElementTree.parse(xml)
        for pkg in tree.findall('DocumentSummary'):
            bs = BioProjectParse(pkg)
            Ncbi.objects(sra__study__BioProject=bs.bioproject['bioproject_accn']).update(bioproject=bs.bioproject)
            Ncbi.objects(sra__study__BioProject=bs.bioproject['bioproject_id']).update(bioproject=bs.bioproject)


def parse_pubmed(cache):
    for xml in cache:
        logger.debug('Parsing: {}'.format(xml))
        tree = ElementTree.parse(xml)
        for pkg in tree.findall('PubmedArticle/MedlineCitation'):
            pm = PubmedParse(pkg)
            Ncbi.objects(sra__study__pubmed=pm.pubmed['pubmed_id']).update(add_to_set__pubmed=pm.pubmed)
            Ncbi.objects(sra__sample__pubmed=pm.pubmed['pubmed_id']).update(add_to_set__pubmed=pm.pubmed)
            Ncbi.objects(sra__experiment__pubmed=pm.pubmed['pubmed_id']).update(add_to_set__pubmed=pm.pubmed)


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
    if args.host:
        logger.info('Connecting to running server at: mongodb://{}:{}'.format(args.host, args.port))
        host = args.host
        mongodb = False
    else:
        logger.info('Starting MongoDB')
        mongodb = MongoDB2(dbDir=args.dbDir, logDir=args.logDir, port=args.port)
        host = '127.0.0.1'
        logger.info('Running server at: mongodb://{}:{}'.format(host, args.port))

    logger.info('Connecting to: {}'.format(args.db))
    client = connect(args.db, host=host, port=args.port)

    try:
        # SRA
        logger.info('Querying SRA: {}'.format(args.query))
        sra_query = ncbi_query(args.query)

        logger.info('Downloading documents')
        cache = Cache(directory='.cache/sra2mongo/sra', clean=args.force)

        logger.info('Saving to cache: {}'.format(cache.cachedir))
        fetch_ncbi(sra_query, cache)

        logger.info('Parsing SRA XML')
        parse_sra(cache)

        # Query BioSample
        cache = Cache(directory='.cache/sra2mongo/biosample')
        accn = list(Ncbi.objects.distinct('sra.sample.BioSample'))

        logger.info('Querying BioSample: with {:,} ids'.format(len(accn)))
        logger.info('Saving to cache: {}'.format(cache.cachedir))
        query = ncbi_query(accn, db='biosample')

        logger.info('Downloading BioSample documents')
        logger.info('Saving to cache: {}'.format(cache.cachedir))
        fetch_ncbi(query, cache, runinfo_retmode=None, db='biosample')

        logger.info('Adding documents to database')
        parse_biosample(cache)

        # Query BioProject
        cache = Cache(directory='.cache/sra2mongo/bioproject')
        accn = list(Ncbi.objects.distinct('sra.study.BioProject'))

        logger.info('Querying BioProject: with {:,} ids'.format(len(accn)))
        logger.info('Saving to cache: {}'.format(cache.cachedir))
        query = ncbi_query(accn, db='bioproject')

        logger.info('Downloading BioProject documents')
        logger.info('Saving to cache: {}'.format(cache.cachedir))
        fetch_ncbi(query, cache, runinfo_retmode=None, db='bioproject')

        logger.info('Adding documents to database')
        parse_bioproject(cache)

        # Query Pubmed
        cache = Cache(directory='.cache/sra2mongo/pubmed')
        accn = list(Ncbi.objects.distinct('sra.study.pubmed'))
        accn.extend(list(Ncbi.objects.distinct('sra.experiment.pubmed')))
        accn.extend(list(Ncbi.objects.distinct('sra.sample.pubmed')))
        accn = list(set(accn))

        logger.info('Querying Pubmed: with {:,} ids'.format(len(accn)))
        logger.info('Saving to cache: {}'.format(cache.cachedir))
        query = ncbi_query(accn, db='pubmed')

        logger.info('Downloading Pubmed documents')
        logger.info('Saving to cache: {}'.format(cache.cachedir))
        fetch_ncbi(query, cache, runinfo_retmode=None, db='pubmed')

        logger.info('Adding documents to database')
        parse_pubmed(cache)
    except Exception as err:
        raise err
    finally:
        client.close()
        if mongodb:
            mongodb.close()


if __name__ == '__main__':
    main()
