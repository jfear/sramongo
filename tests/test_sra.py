#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.sra import SraExperiment

@pytest.fixture(scope="session")
def sra_etree():
    """ Element tree of single sra experiment. """
    fname = 'data/sra_SRR3001915.xml'
    tree = ElementTree.parse(fname)
    root = tree.getroot()
    return root.find('EXPERIMENT_PACKAGE')

@pytest.fixture(scope="session")
def sraTree(sra_etree):
    return SraExperiment(sra_etree)

def test_parse_organization(sraTree):
    assert sraTree.organization['type'] == 'center'
    assert sraTree.organization['abbreviation'] == 'GEO'
    assert sraTree.organization['name'] == 'NCBI'
    assert sraTree.organization['email'] == 'geo-group@ncbi.nlm.nih.gov'
    assert sraTree.organization['first_name'] == 'Geo'
    assert sraTree.organization['last_name'] == 'Curators'

def test_parse_submission(sraTree):
    assert sraTree.submission['submission_id'] == 'SRA178685'
    assert sraTree.submission['submitter_id'][0]['db'] == 'GEO'
    assert sraTree.submission['submitter_id'][0]['id'] == 'GEO: GSE60314'

def test_parse_study(sraTree):
    assert sraTree.study['study_id'] == 'SRP045429'
    assert sraTree.study['external_id'][0]['db'] == 'BioProject'
    assert sraTree.study['external_id'][0]['id'] == 'PRJNA258012'
    assert sraTree.study['external_id'][1]['db'] == 'GEO'
    assert sraTree.study['external_id'][1]['id'] == 'GSE60314'
    assert sraTree.study['study_title'].strip().split(' ')[0] == 'mRNA'
    assert sraTree.study['study_title'].strip().split(' ')[-1] == 'environments'
    assert sraTree.study['study_type'] == 'Transcriptome Analysis'
    assert sraTree.study['study_abstract'].strip().split(' ')[0] == 'Our'
    assert sraTree.study['study_abstract'].strip().split(' ')[-1] == 'level.'
    assert sraTree.study['center_project_name'] == 'GSE60314'
    assert sraTree.study['study_xref_link'][0]['db'] == 'pubmed'
    assert sraTree.study['study_xref_link'][0]['id'] == '26732976'

def test_parse_experiment(sraTree):
    assert sraTree.experiment['experiment_id'] == 'SRX1483046'
    assert sraTree.experiment['submitter_id'][0]['db'] == 'GEO'
    assert sraTree.experiment['submitter_id'][0]['id'] == 'GSM1973128'
    assert sraTree.experiment['title'].strip() == 'GSM1973128: 563M322; Drosophila melanogaster; RNA-Seq'
    assert sraTree.experiment['library_strategy'] == 'RNA-Seq'
    assert sraTree.experiment['library_source'] == 'TRANSCRIPTOMIC'
    assert sraTree.experiment['library_selection'] == 'cDNA'
    assert sraTree.experiment['library_layout'] == 'SINGLE'
    assert sraTree.experiment['library_construction_protocol'].strip().split(' ')[0] == 'Single'
    assert sraTree.experiment['library_construction_protocol'].strip().split(' ')[-1] == 'GEO_biosample_summary.xls.'
    assert sraTree.experiment['platform'] == 'ILLUMINA'
    assert sraTree.experiment['instrument_model'] == 'Illumina HiSeq 2000'
    assert sraTree.experiment['experiment_xref_link'][0]['db'] == 'gds'
    assert sraTree.experiment['experiment_xref_link'][0]['id'] == '301973128'
    assert sraTree.experiment['experiment_attribute']['GEO Accession'] == 'GSM1973128'

