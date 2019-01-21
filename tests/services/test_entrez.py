import os
import datetime
from textwrap import dedent

import pytest

from sramongo.services import entrez

DB = 'sra'
QUERY = '"Drosophila melanogaster"[orgn]'
RETMAX = 2
API_KEY = os.environ.get('ENTREZ_API_KEY', False)


@pytest.fixture(scope='module')
def small_esearch_results() -> entrez.EsearchResult:
    esearch_results = entrez.esearch(DB, QUERY, retmax=RETMAX, api_key=API_KEY)
    return esearch_results


@pytest.fixture(scope='module')
def experiment_set_xml() -> str:
    xml = dedent("""<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
          <EXPERIMENT_PACKAGE>
            <EXPERIMENT accession="SRX5231949" alias="Library_2">
              <IDENTIFIERS>
                <PRIMARY_ID>SRX5231949</PRIMARY_ID>
              </IDENTIFIERS>
            </EXPERIMENT>
            <SUBMISSION accession="SRA832165" alias="SUB5027898">
              <IDENTIFIERS>
                <PRIMARY_ID>SRA832165</PRIMARY_ID>
              </IDENTIFIERS>
            </SUBMISSION>
          </EXPERIMENT_PACKAGE>
          <EXPERIMENT_PACKAGE>
            <EXPERIMENT accession="SRX5231948" alias="Library_1">
              <IDENTIFIERS>
                <PRIMARY_ID>SRX5231948</PRIMARY_ID>
              </IDENTIFIERS>
              <TITLE>d77ZT8</TITLE>
            </EXPERIMENT>
            <SUBMISSION accession="SRA832165" alias="SUB5027898">
              <IDENTIFIERS>
                <PRIMARY_ID>SRA832165</PRIMARY_ID>
              </IDENTIFIERS>
            </SUBMISSION>
          </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>
    """)
    return xml


def test_urlencode_query():
    import urllib.parse
    cleaned_query = urllib.parse.quote_plus(QUERY)
    assert cleaned_query == '%22Drosophila+melanogaster%22%5Borgn%5D'


def test_esearch_sra_nohistory():
    esearch_results = entrez.esearch(DB, QUERY, userhistory=False, retmax=RETMAX, api_key=API_KEY)
    assert len(esearch_results.ids) == RETMAX
    assert esearch_results.webenv == ''
    assert esearch_results.query_key == ''


def test_esearch_sra_withhistory():
    esearch_results = entrez.esearch(DB, QUERY, userhistory=True, retmax=RETMAX, api_key=API_KEY)
    assert len(esearch_results.ids) == RETMAX
    assert esearch_results.webenv != ''
    assert esearch_results.query_key != ''


def test_epost(small_esearch_results):
    ids = small_esearch_results.ids[:2]
    epost_results = entrez.epost(DB, ids=ids, api_key=API_KEY)
    assert epost_results.query_key == '1'


def test_esummary_no_history(small_esearch_results):
    ids = small_esearch_results.ids
    esummary_results = entrez.esummary(DB, ids, api_key=API_KEY)
    assert len(esummary_results) == RETMAX
    assert esummary_results[0].srx != ''
    assert esummary_results[0].id != ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime


def test_esummary_with_history_retmax(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    esummary_results = entrez.esummary(DB, webenv=webenv, query_key=query_key, retmax=600, api_key=API_KEY)
    assert len(esummary_results) == 600
    assert esummary_results[0].srx != ''
    assert esummary_results[0].id != ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime


def test_esummary_with_history_count(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    esummary_results = entrez.esummary(DB, webenv=webenv, query_key=query_key, count=600, api_key=API_KEY)
    assert len(esummary_results) == 600


def test_parse_efetch_experiment_set(experiment_set_xml):
    for experiment in entrez.parse_efetch_experiment_set(experiment_set_xml):
        if experiment.srx == 'SRX5231949':
            assert 'Library_2' in experiment.xml
        else:
            assert 'Library_1' in experiment.xml


def test_efetch_no_history(small_esearch_results):
    ids = small_esearch_results.ids
    for results in entrez.efetch(DB, ids, api_key=API_KEY):
        assert results.srx.startswith('SRX') | results.srx.startswith('DRX') | results.srx.startswith('ERX')
