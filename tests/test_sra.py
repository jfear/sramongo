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
        assert sraTree.organization['organization_type'] == 'center'
        assert sraTree.organization['abbreviation'] == 'GEO'
        assert sraTree.organization['name'] == 'NCBI'
        assert sraTree.organization['email'] == 'geo-group@ncbi.nlm.nih.gov'
        assert sraTree.organization['first_name'] == 'Geo'
        assert sraTree.organization['last_name'] == 'Curators'

    def test_parse_submission(self, sraTree):
        assert sraTree.submission['submission_id'] == 'SRA178685'
        assert sraTree.submission['submitter_id'][0]['id'] == 'GEO: GSE60314'

    def test_parse_study(self, sraTree):
        assert sraTree.study['study_id'] == 'SRP045429'
        assert sraTree.study['BioProject'] == 'PRJNA258012'
        assert sraTree.study['GEO'] == 'GSE60314'
        assert sraTree.study['title'].strip().split(' ')[0] == 'mRNA'
        assert sraTree.study['title'].strip().split(' ')[-1] == 'environments'
        assert sraTree.study['study_type'] == 'Transcriptome Analysis'
        assert sraTree.study['abstract'].strip().split(' ')[0] == 'Our'
        assert sraTree.study['abstract'].strip().split(' ')[-1] == 'level.'
        assert sraTree.study['center_project_name'] == 'GSE60314'
        assert sraTree.study['pubmed'] == '26732976'

    def test_parse_experiment(self, sraTree):
        assert sraTree.experiment['experiment_id'] == 'SRX1483046'
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
        assert sraTree.experiment['GEO_Dataset'] == '301973128'
        assert sraTree.experiment['attributes']['geo accession'] == 'GSM1973128'

    def test_parse_sample(self, sraTree):
        assert sraTree.sample['sample_id'] == 'SRS679015'
        assert sraTree.sample['BioSample'] == 'SAMN02981965'
        assert sraTree.sample['GEO'] == 'GSM1471477'
        assert sraTree.sample['title'] == 'DGRP563 M_E3_2_L3'
        assert sraTree.sample['taxon_id'] == '7227'
        assert sraTree.sample['scientific_name'] == 'Drosophila melanogaster'
        assert sraTree.sample['BioProject'] == '258012'
        assert sraTree.sample['attributes']['source_name'] == 'whole body'
        assert sraTree.sample['attributes']['strain'] == 'dgrp-563'
        assert sraTree.sample['attributes']['developmental stage'] == 'adult'
        assert sraTree.sample['attributes']['sex'] == 'male'
        assert sraTree.sample['attributes']['tissue'] == 'whole body'

    def test_parse_pool(self, sraTree):
        assert sraTree.samples[0] == 'SRS679015'

    def test_run(self, sraTree):
        assert sraTree.run[0]['run_id'] == 'SRR3001915'
        assert sraTree.run[0]['samples'][0] == 'SRS679015'
        assert sraTree.run[0]['experiment_id'] == 'SRX1483046'
        assert sraTree.run[0]['nspots'] == 2001211
        assert sraTree.run[0]['nbases'] == 152092036
        assert sraTree.run[0]['nreads'] == 1
        assert sraTree.run[0]['read_count_r1'] == 2001211.0
        assert sraTree.run[0]['read_len_r1'] == 76.0
        assert sraTree.run[0]['tax_analysis']['tax_counts']['Neoptera']['parent'] == 'Pterygota'
        assert sraTree.run[0]['tax_analysis']['tax_counts']['Neoptera']['self_count'] == 95
        assert sraTree.run[0]['tax_analysis']['tax_counts']['Neoptera']['total_count'] == 360608
        assert sraTree.run[0]['tax_analysis']['tax_counts']['Neoptera']['tax_id'] == '33340'
        assert sraTree.run[0]['tax_analysis']['tax_counts']['Neoptera']['rank'] == 'subclass'


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
        assert sraTree.organization['organization_type'] == 'center'
        assert sraTree.organization['abbreviation'] == 'GEO'
        assert sraTree.organization['name'] == 'NCBI'
        assert sraTree.organization['email'] == 'geo-group@ncbi.nlm.nih.gov'
        assert sraTree.organization['first_name'] == 'Geo'
        assert sraTree.organization['last_name'] == 'Curators'

    def test_parse_submission(self, sraTree):
        assert sraTree.submission['submission_id'] == 'SRA502713'
        assert sraTree.submission['submitter_id'][0]['id'] == 'GEO: GSE92305'

    def test_parse_study(self, sraTree):
        assert sraTree.study['study_id'] == 'SRP094985'
        assert sraTree.study['BioProject'] == 'PRJNA357158'
        assert sraTree.study['GEO'] == 'GSE92305'
        assert sraTree.study['title'].strip().split(' ')[0] == 'Characterization'
        assert sraTree.study['title'].strip().split(' ')[-1] == 'gametogenesis'
        assert sraTree.study['study_type'] == 'Transcriptome Analysis'
        assert sraTree.study['abstract'].strip().split(' ')[0] == 'Total'
        assert sraTree.study['abstract'].strip().split(' ')[-1] == 'testes'
        assert sraTree.study['center_project_name'] == 'GSE92305'

    def test_parse_experiment(self, sraTree):
        assert sraTree.experiment['experiment_id'] == 'SRX2416970'
        assert sraTree.experiment['submitter_id'][0]['id'] == 'GSM2425536'
        assert sraTree.experiment['title'].strip().split(' ')[1] == 'Drosophila'
        assert sraTree.experiment['title'].strip().split(' ')[-1] == 'RNA-Seq'
        assert sraTree.experiment['library_strategy'] == 'RNA-Seq'
        assert sraTree.experiment['library_source'] == 'TRANSCRIPTOMIC'
        assert sraTree.experiment['library_selection'] == 'cDNA'
        assert sraTree.experiment['library_layout'] == 'PAIRED'
        assert sraTree.experiment['library_construction_protocol'].strip().split(' ')[0] == 'Total'
        assert sraTree.experiment['library_construction_protocol'].strip().split(' ')[-1] == 'Center.'
        assert sraTree.experiment['platform'] == 'ILLUMINA'
        assert sraTree.experiment['instrument_model'] == 'Illumina HiSeq 2500'
        assert sraTree.experiment['GEO_Dataset'] == '302425536'
        assert sraTree.experiment['attributes']['geo accession'] == 'GSM2425536'

    def test_parse_sample(self, sraTree):
        assert sraTree.sample['sample_id'] == 'SRS1854358'
        assert sraTree.sample['BioSample'] == 'SAMN06134387'
        assert sraTree.sample['GEO'] == 'GSM2425536'
        assert sraTree.sample['title'] == 'Drosophila melanogaster w1118 testis total RNA replicate 1'
        assert sraTree.sample['taxon_id'] == '7227'
        assert sraTree.sample['scientific_name'] == 'Drosophila melanogaster'
        assert sraTree.sample['BioProject'] == '357158'
        assert sraTree.sample['attributes']['source_name'] == 'dissected testes'
        assert sraTree.sample['attributes']['tissue'] == 'testes'
        assert sraTree.sample['attributes']['genotype'] == 'W1118'

    def test_parse_pool(self, sraTree):
        assert sraTree.samples[0] == 'SRS1854358'

    def test_run(self, sraTree):
        assert sraTree.run[0]['run_id'] == 'SRR5100239'
        assert sraTree.run[0]['samples'][0] == 'SRS1854358'
        assert sraTree.run[0]['experiment_id'] == 'SRX2416970'
        assert sraTree.run[0]['nspots'] == 25818690
        assert sraTree.run[0]['nbases'] == 2581869000
        assert sraTree.run[0]['nreads'] == 2
        assert sraTree.run[0]['read_count_r1'] == 25818690.0
        assert sraTree.run[0]['read_len_r1'] == 50.0
        assert sraTree.run[0]['read_count_r2'] == 25818690.0
        assert sraTree.run[0]['read_len_r2'] == 50.0
