#!/usr/bin/env python
from xml.etree import ElementTree

from dateutil.parser import parse as dateutil_parse
import pytest

from sramongo import parsers_pubmed_xml


def test_parse_pubmed(pubmed_xml):
    root = pubmed_xml
    pubmed = parsers_pubmed_xml.parse_pubmed(root)
    assert pubmed.accn == 26732976
    assert pubmed.title[:10] == 'Comparison'
    assert pubmed.journal == 'BMC Genomics'
    assert pubmed.year == '2016'
    assert pubmed.volume == '17'
    assert pubmed.page == '28'
    assert pubmed.abstract[:11] == 'A generally'
    assert pubmed.abstract[-11:] == 'count data.'
    assert len(pubmed.authors) == 8
    assert pubmed.authors[0]['first_name'] == 'Yanzhu'
    assert pubmed.authors[0]['last_name'] == 'Lin'
    assert pubmed.authors[-1]['first_name'] == 'Susan T'
    assert pubmed.authors[-1]['last_name'] == 'Harbison'
    assert pubmed.authors[0]['affiliation'][:30] == 'Laboratory of Systems Genetics'
    assert pubmed.citation == 'Lin, Yanzhu, et al. BMC Genomics 17 (2016): 28.'
    assert pubmed.date_created == dateutil_parse('2016-01-06')
    assert pubmed.date_completed == dateutil_parse('2016-09-28')
    assert pubmed.date_revised == dateutil_parse('2016-10-19')

