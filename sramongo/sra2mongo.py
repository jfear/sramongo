#!/usr/bin/env python
"""Downloads and parses various Entrez databases."""
import os
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from urllib.error import HTTPError
from xml.etree import ElementTree
from logging import INFO, DEBUG

from Bio import Entrez
from mongoengine import connect

from sramongo.logger import logger
from sramongo.mongo import MongoDB
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

    return parser.parse_args(['--db', '/data/db', '--dbLog', '/data/log', '--email', 'justin.fear@nih.gov',
                              '--query', '"Drosophila melanogaster"[Orgn]', '--debug'])


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

def query_sra(query):
    handle = Entrez.esearch(db='sra', term=query, usehistory="y")
    records = Entrez.read(handle)
    handle.close()
    return records

def fetch_sra(records):
    handle = Entrez.efetch(db='sra', retmax=10,
                           webenv=records['WebEnv'],
                           query_key=records['QueryKey'])

    records = handle.read()
    handle.close()
    return records

def add_pkg_to_database(pkg):
    sraExperiment = SraExperiment(pkg))

    # Add Experiment Package to database
    # Build study document
    submission = Submission(**sraExperiment.submission)
    organization = Organization(**sraExperiment.organization)
    study = Study(submission=submission, organization=organization, **sraExperiment.study)
    study.save()

    # Build sample document
    sample = Sample(**sraExperiment.sample)
    sample.save()

    # Build experiment document
    experiment = Experiment(study=study, **sraExperiment.experiment)
    experiment.samples = sraExperiment.pool
    experiment.save()

    # Add experiment id to study
    study.experiments.append(experiment.id)

    # Build run documents
    for r in sraExperiment.run:
        run = Run(experiment=experiment, **r)
        run.save()
        experiment.runs.append(run)
    experiment.save()

def main():
    # Import commandline arguments.
    args = arguments()

    # Set logging level
    if args.debug:
        logger.setLevel(DEBUG)
    else:
        logger.setLevel(INFO)

    Entrez.email = args.email

    logger.info('Starting MongoDB')
    with MongoDB(dbDir=args.dbDir, logDir=args.logDir, port=args.port):
        client = connect('sra')
        try:
            sra_query = query_sra(args.query)
            sra_records = fetch_sra(sra_query)

            tree = ElementTree.fromstring(sra_records)
            sras = []
            # Iterate over experiment packages
            for exp_pkg in tree:
        except:
            pass
        finally:
           client.close()
