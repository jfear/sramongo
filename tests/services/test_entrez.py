import os
import datetime
from textwrap import dedent

import pytest

from sramongo import parsers_pubmed_xml, parsers_biosample_xml, parsers_bioproject_xml, parsers_sra_xml
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
    cleaned_query = urllib.parse.quote_plus(QUERY, safe='/+')
    assert cleaned_query == '%22Drosophila+melanogaster%22%5Borgn%5D'


@pytest.mark.skip
def test_esearch_sra_nohistory():
    esearch_results = entrez.esearch(DB, QUERY, userhistory=False, retmax=RETMAX, api_key=API_KEY)
    assert len(esearch_results.ids) == RETMAX
    assert esearch_results.webenv == ''
    assert esearch_results.query_key == ''


@pytest.mark.skip
def test_esearch_sra_withhistory():
    esearch_results = entrez.esearch(DB, QUERY, userhistory=True, retmax=RETMAX, api_key=API_KEY)
    assert len(esearch_results.ids) == RETMAX
    assert esearch_results.webenv != ''
    assert esearch_results.query_key != ''


@pytest.mark.skip
def test_epost(small_esearch_results):
    ids = small_esearch_results.ids[:2]
    epost_results = entrez.epost(DB, ids=ids, api_key=API_KEY)
    assert epost_results.query_key == '1'


@pytest.mark.skip
def test_esummary_sra_no_history(small_esearch_results):
    ids = small_esearch_results.ids
    esummary_results = []
    for docs in entrez.esummary(DB, ids, api_key=API_KEY):
        esummary_results.extend(list(parsers_sra_xml.parse_sra_esummary_result(docs)))
    assert len(esummary_results) == RETMAX
    assert esummary_results[0].accn != ''
    assert esummary_results[0].accn != ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime


@pytest.mark.skip
def test_esummary_sra_with_history_retmax(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    esummary_results = []
    for docs in entrez.esummary(DB, webenv=webenv, query_key=query_key, retmax=600, api_key=API_KEY):
        esummary_results.extend(list(parsers_sra_xml.parse_sra_esummary_result(docs)))
    assert len(esummary_results) == 600
    assert esummary_results[0].accn != ''
    assert esummary_results[0].accn!= ''
    assert type(esummary_results[0].create_date) == datetime.datetime
    assert type(esummary_results[0].update_date) == datetime.datetime


@pytest.mark.skip
def test_esummary_sra_with_history_count(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    esummary_results = []
    for docs in entrez.esummary(DB, webenv=webenv, query_key=query_key, count=600, api_key=API_KEY):
        esummary_results.extend(list(parsers_sra_xml.parse_sra_esummary_result(docs)))
    assert len(esummary_results) == 600


def test_parse_efetch_result(experiment_set_xml):
    for experiment in parsers_sra_xml.parse_sra_efetch_result(experiment_set_xml):
        if experiment.accn == 'SRX5231949':
            assert 'Library_2' in experiment.xml
        else:
            assert 'Library_1' in experiment.xml


@pytest.mark.skip
def test_efetch_no_history(small_esearch_results):
    ids = small_esearch_results.ids
    for result in entrez.efetch(DB, ids, api_key=API_KEY):
        for experiment in parsers_sra_xml.parse_sra_efetch_result(result):
            assert experiment.accn.startswith('SRX') | experiment.accn.startswith('DRX') | experiment.accn.startswith('ERX')


@pytest.mark.skip
def test_elink_no_hisotry(small_esearch_results):
    ids = small_esearch_results.ids
    result = entrez.elink(db='biosample', dbfrom='sra', ids=ids)
    assert result.dbfrom == 'sra'
    assert result.dbto == 'biosample'
    assert result.query_key == '1'


@pytest.mark.skip
def test_elink_with_hisotry(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    result = entrez.elink(db='biosample', dbfrom='sra', webenv=webenv, query_key=query_key)
    assert result.dbfrom == 'sra'
    assert result.dbto == 'biosample'


@pytest.mark.skip
def test_elink_no_hisotry_no_results(small_esearch_results):
    ids = small_esearch_results.ids
    result = entrez.elink(db='pubmed', dbfrom='sra', ids=ids)
    assert result.dbfrom == 'sra'
    assert result.dbto == ''
    assert result.query_key == ''


@pytest.mark.skip
def test_efetch_bioproject(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    link = entrez.elink('bioproject', 'sra', webenv=webenv, query_key=query_key, api_key=API_KEY, retmax=RETMAX)
    for result in entrez.efetch('bioproject', webenv=link.webenv, query_key=link.query_key, api_key=API_KEY, retmax=RETMAX):
        for document in parsers_bioproject_xml.parse_bioproject_efetch_result(result):
            assert document.accn.startswith('PRJ')


@pytest.mark.skip
def test_efetch_biosample(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    link = entrez.elink('biosample', 'sra', webenv=webenv, query_key=query_key, api_key=API_KEY, retmax=RETMAX)
    for result in entrez.efetch('biosample', webenv=link.webenv, query_key=link.query_key, api_key=API_KEY, retmax=RETMAX):
        for document in parsers_biosample_xml.parse_biosample_efetch_result(result):
            assert document.accn.startswith('SAMN')


@pytest.mark.skip
def test_efetch_pubmed(small_esearch_results):
    webenv = small_esearch_results.webenv
    query_key = small_esearch_results.query_key
    link = entrez.elink('pubmed', 'sra', webenv=webenv, query_key=query_key, api_key=API_KEY, retmax=RETMAX)
    for result in entrez.efetch('pubmed', webenv=link.webenv, query_key=link.query_key, api_key=API_KEY, retmax=RETMAX):
        for document in parsers_pubmed_xml.parse_pubmed_efetch_result(result):
            assert isinstance(document.accn, str)
