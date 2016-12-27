#!/usr/bin/env python
"""Downloads and parses various Entrez databases."""
import os
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from urllib.error import HTTPError
from xml.etree import ElementTree

from Bio import Entrez
from mongoengine import connect

from sramongo.logger import logger
from sramongo.mongo import start_mongo
from sramongo.sra import SraExperiment
from sramongo.mongo_schema import Experiment, Sample, Organization, Run, Study


def arguments():
    """Pulls in command line arguments."""

    DESCRIPTION = """\
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=Raw)

    parser.add_argument("--email", dest="email", action='store', required=True,
                        help="An email address is required for querying Entrez databases.")

    parser.add_argument("--query", dest="query", action='store', required=True,
                        help="Query to submit to Entrez.")

    parser.add_argument("--db", dest="db", action='store', required=True,
                        help="Folder containing mongo database.")

    parser.add_argument("--dbLog", dest="dbLog", action='store', required=False, default=None,
                        help="Folder containing mongo database.")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False,
                        help="Turn on debug output.")

    return parser.parse_args(['--db', '/data/db', '--dbLog', '/data/log', '--email', 'justin.fear@nih.gov', '--query', '"Drosophila melanogaster"[Orgn]', '--debug'])

def download_sra(records):
    records = sra_records

    batch_size = 3
    out_handle = open("orchid_rpl16.fasta", "w")
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                fetch_handle = Entrez.efetch(db="nucleotide",
                                             rettype="fasta", retmode="text",
                                             retstart=start, retmax=batch_size,
                                             webenv=webenv, query_key=query_key)
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(15)
                else:
                    raise
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()

def main():
    # Import commandline arguments.
    args = arguments()
    Entrez.email = args.email
    query = args.query

    # Start Mongo
    db = args.db
    dblog = args.dbLog
    if not os.path.exists(db):
        os.mkdir(db)

    dblog = args.dbLog
    if not dblog is None:
        if not os.path.exists(dblog):
            os.mkdir(dblog)

    pid = start_mongo(dbDir=db, logDir=dblog)
    client = connect('sra')

    try:
        # query SRA
        search_handle = Entrez.esearch(db='sra', retmax=10, term=query, usehistory="y")
        search_records = Entrez.read(search_handle)
        search_handle.close()

        webenv = search_records['WebEnv']
        query_key = search_records['QueryKey']

        fetch_handle = Entrez.efetch(db='sra', retmax=10, webenv=webenv, query_key=query_key)
        fetch_records = fetch_handle.read()
        fetch_handle.close()

        tree = ElementTree.fromstring(fetch_records)
        for pkg in tree:
            sra = SraExperiment(pkg)
            study = Study(**sra.study).save()
            org = Organization(**sra.organization).save()
    except:
        pass
    finally:
       client.close()
       pid.kill()
