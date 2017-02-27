#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.sra import SraExperiment


class TestSRR3001915:
    @pytest.fixture(scope='class')
    def sra_etree(self):
        """ Element tree of single sra experiment. """
        fname = 'tests/data/sra_SRR3001915.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('EXPERIMENT_PACKAGE')

    @pytest.fixture(autouse=True, scope='class')
    def sraTree(self, sra_etree):
        return SraExperiment(sra_etree)

    def test_parse_organization(self, sraTree):
        assert sraTree.sra['organization']['organization_type'] == 'center'
        assert sraTree.sra['organization']['abbreviation'] == 'GEO'
        assert sraTree.sra['organization']['name'] == 'NCBI'
        assert sraTree.sra['organization']['email'] == 'geo-group@ncbi.nlm.nih.gov'
        assert sraTree.sra['organization']['first_name'] == 'Geo'
        assert sraTree.sra['organization']['last_name'] == 'Curators'

    def test_parse_submission(self, sraTree):
        assert sraTree.sra['submission']['submission_id'] == 'SRA178685'
        assert sraTree.sra['submission']['submitter_id'][0]['id'] == 'GEO: GSE60314'

    def test_parse_study(self, sraTree):
        assert sraTree.sra['study']['study_id'] == 'SRP045429'
        assert sraTree.sra['study']['BioProject'] == 'PRJNA258012'
        assert sraTree.sra['study']['GEO'] == 'GSE60314'
        assert sraTree.sra['study']['title'].strip().split(' ')[0] == 'mRNA'
        assert sraTree.sra['study']['title'].strip().split(' ')[-1] == 'environments'
        assert sraTree.sra['study']['study_type'] == 'Transcriptome Analysis'
        assert sraTree.sra['study']['abstract'].strip().split(' ')[0] == 'Our'
        assert sraTree.sra['study']['abstract'].strip().split(' ')[-1] == 'level.'
        assert sraTree.sra['study']['center_project_name'] == 'GSE60314'
        assert sraTree.sra['study']['pubmed'] == '26732976'

    def test_parse_experiment(self, sraTree):
        assert sraTree.sra['experiment']['experiment_id'] == 'SRX1483046'
        assert sraTree.sra['experiment']['submitter_id'][0]['id'] == 'GSM1973128'
        assert sraTree.sra['experiment']['title'].strip() == 'GSM1973128: 563M322; Drosophila melanogaster; RNA-Seq'
        assert sraTree.sra['experiment']['library_strategy'] == 'RNA-Seq'
        assert sraTree.sra['experiment']['library_source'] == 'TRANSCRIPTOMIC'
        assert sraTree.sra['experiment']['library_selection'] == 'cDNA'
        assert sraTree.sra['experiment']['library_layout'] == 'SINGLE'
        assert sraTree.sra['experiment']['library_construction_protocol'].strip().split(' ')[0] == 'Single'
        assert sraTree.sra['experiment']['library_construction_protocol'].strip().split(' ')[-1] == 'GEO_biosample_summary.xls.'
        assert sraTree.sra['experiment']['platform'] == 'ILLUMINA'
        assert sraTree.sra['experiment']['instrument_model'] == 'Illumina HiSeq 2000'
        assert sraTree.sra['experiment']['GEO_Dataset'] == '301973128'
        assert sraTree.sra['experiment']['attributes'][0]['name'] == 'GEO Accession'
        assert sraTree.sra['experiment']['attributes'][0]['value'] == 'GSM1973128'

    def test_parse_sample(self, sraTree):
        assert sraTree.sra['sample']['sample_id'] == 'SRS679015'
        assert sraTree.sra['sample']['BioSample'] == 'SAMN02981965'
        assert sraTree.sra['sample']['GEO'] == 'GSM1471477'
        assert sraTree.sra['sample']['title'] == 'DGRP563 M_E3_2_L3'
        assert sraTree.sra['sample']['taxon_id'] == '7227'
        assert sraTree.sra['sample']['scientific_name'] == 'Drosophila melanogaster'
        assert sraTree.sra['sample']['BioProject'] == '258012'

        attrs = {
                'source_name': 'Whole body',
                'strain': 'DGRP-563',
                'developmental stage': 'Adult',
                'Sex': 'male',
                'tissue': 'Whole body',
                }

        for attribute in sraTree.sra['sample']['attributes']:
            name = attribute['name']
            assert attrs[name] == attribute['value']

    def test_parse_pool(self, sraTree):
        assert sraTree.sra['pool'][0] == 'SRS679015'

    def test_run(self, sraTree):
        assert sraTree.sra['run'][0]['run_id'] == 'SRR3001915'
        assert sraTree.sra['run'][0]['samples'][0] == 'SRS679015'
        assert sraTree.sra['run'][0]['experiment_id'] == 'SRX1483046'
        assert sraTree.sra['run'][0]['nspots'] == 2001211
        assert sraTree.sra['run'][0]['nbases'] == 152092036
        assert sraTree.sra['run'][0]['nreads'] == 1
        assert sraTree.sra['run'][0]['read_count_r1'] == 2001211.0
        assert sraTree.sra['run'][0]['read_len_r1'] == 76.0

        assert sraTree.sra['run'][0]['tax_analysis']['tax_counts']['subclass'][0]['parent'] == 'Pterygota'
        assert sraTree.sra['run'][0]['tax_analysis']['tax_counts']['subclass'][0]['self_count'] == 95
        assert sraTree.sra['run'][0]['tax_analysis']['tax_counts']['subclass'][0]['total_count'] == 360608
        assert sraTree.sra['run'][0]['tax_analysis']['tax_counts']['subclass'][0]['tax_id'] == '33340'
        assert sraTree.sra['run'][0]['tax_analysis']['tax_counts']['subclass'][0]['name'] == 'Neoptera'


class TestSRR5100239:
    @pytest.fixture(scope='class')
    def sra_etree(self):
        """ Element tree of single sra experiment. """
        fname = 'tests/data/sra_SRR5100239.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('EXPERIMENT_PACKAGE')

    @pytest.fixture(scope='class')
    def sraTree(self, sra_etree):
        return SraExperiment(sra_etree)

    def test_parse_organization(self, sraTree):
        assert sraTree.sra['organization']['organization_type'] == 'center'
        assert sraTree.sra['organization']['abbreviation'] == 'GEO'
        assert sraTree.sra['organization']['name'] == 'NCBI'
        assert sraTree.sra['organization']['email'] == 'geo-group@ncbi.nlm.nih.gov'
        assert sraTree.sra['organization']['first_name'] == 'Geo'
        assert sraTree.sra['organization']['last_name'] == 'Curators'

    def test_parse_submission(self, sraTree):
        assert sraTree.sra['submission']['submission_id'] == 'SRA502713'
        assert sraTree.sra['submission']['submitter_id'][0]['id'] == 'GEO: GSE92305'

    def test_parse_study(self, sraTree):
        assert sraTree.sra['study']['study_id'] == 'SRP094985'
        assert sraTree.sra['study']['BioProject'] == 'PRJNA357158'
        assert sraTree.sra['study']['GEO'] == 'GSE92305'
        assert sraTree.sra['study']['title'].strip().split(' ')[0] == 'Characterization'
        assert sraTree.sra['study']['title'].strip().split(' ')[-1] == 'gametogenesis'
        assert sraTree.sra['study']['study_type'] == 'Transcriptome Analysis'
        assert sraTree.sra['study']['abstract'].strip().split(' ')[0] == 'Total'
        assert sraTree.sra['study']['abstract'].strip().split(' ')[-1] == 'testes'
        assert sraTree.sra['study']['center_project_name'] == 'GSE92305'

    def test_parse_experiment(self, sraTree):
        assert sraTree.sra['experiment']['experiment_id'] == 'SRX2416970'
        assert sraTree.sra['experiment']['submitter_id'][0]['id'] == 'GSM2425536'
        assert sraTree.sra['experiment']['title'].strip().split(' ')[1] == 'Drosophila'
        assert sraTree.sra['experiment']['title'].strip().split(' ')[-1] == 'RNA-Seq'
        assert sraTree.sra['experiment']['library_strategy'] == 'RNA-Seq'
        assert sraTree.sra['experiment']['library_source'] == 'TRANSCRIPTOMIC'
        assert sraTree.sra['experiment']['library_selection'] == 'cDNA'
        assert sraTree.sra['experiment']['library_layout'] == 'PAIRED'
        assert sraTree.sra['experiment']['library_construction_protocol'].strip().split(' ')[0] == 'Total'
        assert sraTree.sra['experiment']['library_construction_protocol'].strip().split(' ')[-1] == 'Center.'
        assert sraTree.sra['experiment']['platform'] == 'ILLUMINA'
        assert sraTree.sra['experiment']['instrument_model'] == 'Illumina HiSeq 2500'
        assert sraTree.sra['experiment']['GEO_Dataset'] == '302425536'
        assert sraTree.sra['experiment']['attributes'][0]['name'] == 'GEO Accession'
        assert sraTree.sra['experiment']['attributes'][0]['value'] == 'GSM2425536'

    def test_parse_sample(self, sraTree):
        assert sraTree.sra['sample']['sample_id'] == 'SRS1854358'
        assert sraTree.sra['sample']['BioSample'] == 'SAMN06134387'
        assert sraTree.sra['sample']['GEO'] == 'GSM2425536'
        assert sraTree.sra['sample']['title'] == 'Drosophila melanogaster w1118 testis total RNA replicate 1'
        assert sraTree.sra['sample']['taxon_id'] == '7227'
        assert sraTree.sra['sample']['scientific_name'] == 'Drosophila melanogaster'
        assert sraTree.sra['sample']['BioProject'] == '357158'

        attrs = {
                'source_name': 'dissected testes',
                'tissue': 'testes',
                'genotype': 'w1118',
                }

        for attribute in sraTree.sra['sample']['attributes']:
            name = attribute['name']
            assert attrs[name] == attribute['value']

    def test_parse_pool(self, sraTree):
        assert sraTree.sra['pool'][0] == 'SRS1854358'

    def test_run(self, sraTree):
        assert sraTree.sra['run'][0]['run_id'] == 'SRR5100239'
        assert sraTree.sra['run'][0]['samples'][0] == 'SRS1854358'
        assert sraTree.sra['run'][0]['experiment_id'] == 'SRX2416970'
        assert sraTree.sra['run'][0]['nspots'] == 25818690
        assert sraTree.sra['run'][0]['nbases'] == 2581869000
        assert sraTree.sra['run'][0]['nreads'] == 2
        assert sraTree.sra['run'][0]['read_count_r1'] == 25818690.0
        assert sraTree.sra['run'][0]['read_len_r1'] == 50.0
        assert sraTree.sra['run'][0]['read_count_r2'] == 25818690.0
        assert sraTree.sra['run'][0]['read_len_r2'] == 50.0
