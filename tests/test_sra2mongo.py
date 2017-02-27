"""Test for the sra2mongo program."""
import os
import pytest
import pandas as pd
from xml.etree import ElementTree

from mongoengine import connect
from Bio import Entrez

from sramongo.mongo_schema import Ncbi
from sramongo.sra2mongo import Cache, ncbi_query, fetch_sra, parse_sra

# TODO: Probably should not set my email permanently.
Entrez.email = 'justin.fear@nih.gov'

@pytest.fixture(scope='session')
def cache(tmpdir_factory):
    tmpdir = str(tmpdir_factory.mktemp('sra2mongo'))
    return Cache(tmpdir)


@pytest.fixture(scope='session')
def sra_query():
    return ncbi_query('ERR1751044')


@pytest.fixture(scope='session')
def sra_fetch(sra_query, cache):
    fetch_sra(sra_query, cache, runinfo_retmode='text')
    return cache


def test_ncbi_query(sra_query):
    assert sra_query['RetMax'] == '1'
    assert sra_query['TranslationStack'][0]['Term'] == 'ERR1751044[All Fields]'


def test_fetch_sra(sra_fetch):
    assert 0 in sra_fetch.cached
    assert os.path.exists(os.path.join(sra_fetch.cachedir, '0.xml'))
    assert os.path.exists(os.path.join(sra_fetch.cachedir, '0.csv'))


def test_parse_sra(sra_fetch, mongoDB):
    connect('test_sra')
    parse_sra(sra_fetch)
    ncbi = Ncbi.objects(sra__run__run_id='ERR1751044').first()
    assert ncbi.sra.run[0].run_id == 'ERR1751044'
    assert ncbi.sra.experiment.experiment_id == 'ERX1819343'
