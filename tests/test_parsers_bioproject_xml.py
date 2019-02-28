from sramongo import parsers_bioproject_xml
from dateutil.parser import parse as dateutil_parse


def test_parse_bioproject(bioproject_xml):
    root = bioproject_xml
    bioproject = parsers_bioproject_xml.parse_bioproject(root)
    assert bioproject.accn == 'PRJNA258012'
    assert bioproject.id == 258012
    assert bioproject.name[:13] == 'mRNA sequence'
    assert bioproject.title[:13] == 'mRNA sequence'
    assert bioproject.description[:11] == 'Our primary'
    assert bioproject.last_update == dateutil_parse('2016-04-20')
    assert bioproject.submission_date == dateutil_parse('2014-08-11')
