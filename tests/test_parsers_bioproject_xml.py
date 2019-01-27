from sramongo import parsers_bioproject_xml


def test_parse_bioproject(bioproject_xml):
    root = bioproject_xml
    bioproject = parsers_bioproject_xml.parse_bioproject(root)
    assert bioproject.accn == 'PRJNA258012'
    assert bioproject.id == 258012
    assert bioproject.name[:13] == 'mRNA sequence'
    assert bioproject.title[:13] == 'mRNA sequence'
    assert bioproject.description[:11] == 'Our primary'
    # TODO make dates using datetime
    assert bioproject.last_update == '2016-04-20'
    assert bioproject.submission_date == '2014-08-11'
