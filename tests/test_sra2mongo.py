"""Test for the sra2mongo program."""
import pytest
import pandas as pd
from io import StringIO
from xml.etree import ElementTree

from mongoengine import connect
from Bio import Entrez

from sramongo.mongo_schema import Experiment, Run
from sramongo.sra2mongo import query_sra, fetch_sra, add_pkg_to_database

# TODO: Probably should not set my email permanently.
Entrez.email = 'justin.fear@nih.gov'


@pytest.fixture(scope='session')
def sra_query():
    return query_sra('ERR1751044')


@pytest.fixture(scope='session')
def sra_fetch(sra_query):
    for xml, runinfo in fetch_sra(sra_query, runinfo_retmode='text'):
        ri = pd.read_csv(StringIO(runinfo), index_col='Run')
        tree = ElementTree.fromstring(xml)
        return tree, ri


def test_query_sra(sra_query):
    assert sra_query['RetMax'] == '1'
    assert sra_query['TranslationStack'][0]['Term'] == 'ERR1751044[All Fields]'


def test_fetch_sra(sra_fetch):
    xml, runinfo = sra_fetch
    assert xml.find('*//EXPERIMENT').get('accession') == 'ERX1819343'
    assert runinfo.index[0] == 'ERR1751044'


def test_add_pkg_to_database(sra_fetch, mongoDB):
    client = connect('test_sra')
    xml, runinfo = sra_fetch

    for pkg in xml:
        add_pkg_to_database(pkg, runinfo)

    run = Run.objects(run_id='ERR1751044').first().select_related(max_depth=9)
    assert run.run_id == 'ERR1751044'
    assert run.experiment_id == 'ERX1819343'
    assert run.experiment.experiment_id == 'ERX1819343'
