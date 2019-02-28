import sramongo.xml_helpers
from pathlib import Path


def test_xml_to_root_from_file_handler():
    fname = Path(__file__).parent / 'data/ERR1662611.xml'
    with open(fname) as fh:
        root = sramongo.xml_helpers.xml_to_root(fh)
    experiment = root.find('EXPERIMENT')
    assert experiment.attrib['accession'] == 'ERX1732932'


def test_xml_to_root_from_string():
    fname = Path(__file__).parent / 'data/ERR1662611.xml'
    with open(fname) as fh:
        xml = fh.read()
    root = sramongo.xml_helpers.xml_to_root(xml)
    experiment = root.find('EXPERIMENT')
    assert experiment.attrib['accession'] == 'ERX1732932'


def test_xml_get_text(sra_xml_root):
    root = sra_xml_root
    srx = sramongo.xml_helpers.get_xml_text(root, 'EXPERIMENT/IDENTIFIERS/PRIMARY_ID')
    assert srx == 'SRX971855'


def test_xml_get_text_invalid_path(sra_xml_root):
    root = sra_xml_root
    result = sramongo.xml_helpers.get_xml_text(root, 'EXPERIMENT/IDENTIFIERS/PRIMARY_ID/NOT_REALLY_HERE')
    assert result == ''
