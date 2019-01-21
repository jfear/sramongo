import pytest

from sramongo import sra_xml
from sramongo import models


@pytest.fixture(scope='session')
def sra_xml_root():
    fname = 'data/SRX971855.xml'
    return sra_xml.xml_to_root(fname)


@pytest.fixture(scope='session')
def sra_xml_root2():
    fname = 'data/SRR3001915.xml'
    return sra_xml.xml_to_root(fname)


def test_xml_to_root_from_file_handler():
    fname = 'data/sra_ERR1662611.xml'
    with open(fname) as fh:
        root = sra_xml.xml_to_root(fh)
    experiment = root.find('EXPERIMENT_PACKAGE/EXPERIMENT')
    assert experiment.attrib['accession'] == 'ERX1732932'


def test_xml_to_root_from_string():
    fname = 'data/sra_ERR1662611.xml'
    with open(fname) as fh:
        xml = fh.read()
    root = sra_xml.xml_to_root(xml)
    experiment = root.find('EXPERIMENT_PACKAGE/EXPERIMENT')
    assert experiment.attrib['accession'] == 'ERX1732932'


def test_xml_get_text(sra_xml_root):
    root = sra_xml_root
    srx = sra_xml.get_xml_text(root, 'EXPERIMENT/IDENTIFIERS/PRIMARY_ID')
    assert srx == 'SRX971855'


def test_xml_get_text_invalid_path(sra_xml_root):
    root = sra_xml_root
    result = sra_xml.get_xml_text(root, 'EXPERIMENT/IDENTIFIERS/PRIMARY_ID/NOT_REALLY_HERE')
    assert result == ''


def test_add_library_layout(sra_xml_root):
    root = sra_xml_root
    sra = models.SraDocument()
    sra_xml.add_library_layout(root, sra)
    assert sra.library_layout == 'SINGLE'


def test_add_platform_information(sra_xml_root):
    root = sra_xml_root
    sra = models.SraDocument()
    sra_xml.add_platform_information(root, sra)
    assert sra.platform == 'ILLUMINA'
    assert sra.instrument_model == 'Illumina Genome Analyzer IIx'


def test_parse_sra_study(sra_xml_root):
    root = sra_xml_root
    study = sra_xml.parse_sra_study(root)
    assert study.accn == 'SRP056660'
    assert study.center_name == 'GEO'


def test_parse_sra_organization(sra_xml_root2):
    root = sra_xml_root2
    organization = sra_xml.parse_sra_organization(root)
    assert organization.organization_type == 'center'


def test_add_sample_attributes(sra_xml_root2):
    root = sra_xml_root2
    sample = models.Sample()
    sra_xml.add_sample_attributes(root, sample)
    assert sample.attributes[0]['name'] == 'source_name'
    assert sample.attributes[0]['value'] == 'Whole body'
