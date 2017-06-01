"""Test for the sra2mongo program."""
import os
import pytest
import pandas as pd
from xml.etree import ElementTree
from io import StringIO

from mongoengine import connect
from Bio import Entrez

from sramongo.mongo_schema import Ncbi
from sramongo.sra2mongo import Cache, ncbi_query, fetch_ncbi, parse_sra, catch_xml_error

# TODO: Probably should not set my email permanently.
Entrez.email = 'justin.fear@nih.gov'

@pytest.fixture
def tmpdir(tmpdir_factory):
    return str(tmpdir_factory.mktemp('sra2mongo'))


@pytest.fixture(scope='session')
def sra_query():
    return ncbi_query('ERR1751044')


@pytest.fixture
def sra_fetch(sra_query, tmpdir):
    cache = Cache(tmpdir)
    fetch_ncbi(sra_query, cache, runinfo_retmode='text')
    return cache


def test_ncbi_query(sra_query):
    assert sra_query['RetMax'] == '1'
    assert sra_query['TranslationStack'][0]['Term'] == 'ERR1751044[All Fields]'


def test_fetch_sra(sra_fetch):
    print(sra_fetch.cached)
    assert '0' in sra_fetch.cached
    assert os.path.exists(os.path.join(sra_fetch.cachedir, '0.xml'))
    assert os.path.exists(os.path.join(sra_fetch.cachedir, '0.csv'))


def test_fetch_sra_remove(sra_fetch):
    assert '0' in sra_fetch.cached
    sra_fetch.remove_from_cache(0)
    assert os.path.exists(os.path.join(sra_fetch.cachedir, '0.xml')) is False
    assert os.path.exists(os.path.join(sra_fetch.cachedir, '0.csv')) is False
    assert '0' not in sra_fetch.cached


def test_parse_sra(sra_fetch, mongoDB):
    connect('test_sra')
    parse_sra(sra_fetch)
    ncbi = Ncbi.objects(sra__run__run_id='ERR1751044').first()
    assert ncbi.sra.run[0].run_id == 'ERR1751044'
    assert ncbi.sra.experiment.experiment_id == 'ERX1819343'


def test_catch_xml_error():
    XMLERROR = """<?xml version="1.0" encoding="UTF-8" ?>
    <!DOCTYPE eEfetchResult PUBLIC "-//NLM//DTD efetch 20131226//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20131226/efetch.dtd">
    <eFetchResult>
        <ERROR>Unable to obtain query #1</ERROR>
    </eFetchResult>"""

    assert catch_xml_error(StringIO(XMLERROR)) is True

def test_catch_xml_empty():
    XMLERROR = ""
    assert catch_xml_error(StringIO(XMLERROR)) is True
