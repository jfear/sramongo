import datetime

import pytest

from sramongo.services import entrez

db = 'sra'
query = '"Drosophila melanogaster"[orgn]'
retmax = 2


@pytest.fixture
def small_esearch_results() -> entrez.EsearchResult:
    esearch_results = entrez.esearch(db, query, retmax=retmax)
    return esearch_results


def test_urlencode_query():
    import urllib.parse
    cleaned_query = urllib.parse.quote_plus(query)
    assert cleaned_query == '%22Drosophila+melanogaster%22%5Borgn%5D'


def test_esearch_sra_nohistory():
    esearch_results = entrez.esearch(db, query, userhistory=False, retmax=retmax)
    assert len(esearch_results.ids) == retmax
    assert esearch_results.webenv == ''
    assert esearch_results.query_key == ''


def test_esearch_sra_withhistory():
    esearch_results = entrez.esearch(db, query, userhistory=True, retmax=retmax)
    assert len(esearch_results.ids) == retmax
    assert esearch_results.webenv != ''
    assert esearch_results.query_key != ''


def test_esummary_nohisotry(small_esearch_results):
    ids = small_esearch_results.ids
    esummary_results = entrez.esummary(db, ids)
    assert len(esummary_results) == retmax
    assert esummary_results[0].srx != ''
    assert esummary_results[0].id != ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime


def test_esummary_withhisotry(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    esummary_results = entrez.esummary(db, webenv=webenv, query_key=query_key)
    assert len(esummary_results) == retmax
    assert esummary_results[0].srx != ''
    assert esummary_results[0].id != ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime

