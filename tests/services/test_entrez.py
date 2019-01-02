import datetime

import pytest

from sramongo.services import entrez

DB = 'sra'
QUERY = '"Drosophila melanogaster"[orgn]'
RETMAX = 2


@pytest.fixture
def small_esearch_results() -> entrez.EsearchResult:
    esearch_results = entrez.esearch(DB, QUERY, retmax=RETMAX)
    return esearch_results


def test_urlencode_query():
    import urllib.parse
    cleaned_query = urllib.parse.quote_plus(QUERY)
    assert cleaned_query == '%22Drosophila+melanogaster%22%5Borgn%5D'


def test_esearch_sra_nohistory():
    esearch_results = entrez.esearch(DB, QUERY, userhistory=False, retmax=RETMAX)
    assert len(esearch_results.ids) == RETMAX
    assert esearch_results.webenv == ''
    assert esearch_results.query_key == ''


def test_esearch_sra_withhistory():
    esearch_results = entrez.esearch(DB, QUERY, userhistory=True, retmax=RETMAX)
    assert len(esearch_results.ids) == RETMAX
    assert esearch_results.webenv != ''
    assert esearch_results.query_key != ''


def test_esummary_nohisotry(small_esearch_results):
    ids = small_esearch_results.ids
    esummary_results = entrez.esummary(DB, ids)
    assert len(esummary_results) == RETMAX
    assert esummary_results[0].srx != ''
    assert esummary_results[0].id != ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime


def test_esummary_withhisotry_retmax(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    esummary_results = entrez.esummary(DB, webenv=webenv, query_key=query_key, retmax=600)
    assert len(esummary_results) == 600
    assert esummary_results[0].srx != ''
    assert esummary_results[0].id != ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime


def test_esummary_withhisotry_count(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    esummary_results = entrez.esummary(DB, webenv=webenv, query_key=query_key, count=600)
    assert len(esummary_results) == 600


def test_efetch_nohisotry(small_esearch_results):
    ids = small_esearch_results.ids
    for results in entrez.efetch(DB, ids):
        # TODO: Add XML parsing here.
        pass
