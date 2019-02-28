from sramongo import parsers_biosample_xml
from dateutil.parser import parse as dateutil_parse


def test_parse_biosample(biosample_xml):
    root = biosample_xml
    biosample = parsers_biosample_xml.parse_biosample(root)
    assert biosample.accn == 'SAMN02981965'
    assert biosample.id == 2981965
    assert biosample.title == 'DGRP563 M_E3_2_L3'
    assert biosample.description == ''
    assert biosample.last_update == dateutil_parse('2015-12-22T01:19:06.303')
    assert biosample.submission_date == dateutil_parse('2014-08-11T13:36:38.303')
    assert biosample.attributes[0]['name'] == 'source_name'
    assert biosample.attributes[0]['value'] == 'Whole body'
