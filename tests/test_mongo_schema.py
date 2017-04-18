"""Tests for mongo_schema.py"""
from xml.etree import ElementTree
import pytest
from textwrap import dedent
from sramongo.sra import SraExperiment
from sramongo.biosample import BioSampleParse
from sramongo.bioproject import BioProjectParse
from sramongo.pubmed import PubmedParse
from sramongo.mongo_schema import URLLink, Xref, XrefLink, \
    EntrezLink, DDBJLink, ENALink, Submission, Organization, \
    Study, Sample, Experiment, Run, BioSample, Ncbi

from mongoengine import connect


def test_URLLink_str_complete():
    url = URLLink(label='test', url='http://test.com')
    assert url.label == 'test'
    assert url.url == 'http://test.com'


def test_URLLink_str_nolabel():
    url = URLLink(url='http://test.com')
    assert url.label is None
    assert url.url == 'http://test.com'


def test_URLLink_str_nourl():
    url = URLLink(label='test')
    assert url.label == 'test'
    assert url.url is None


def test_Xref_str_complete():
    xref = Xref(db='test', id='1030')
    assert xref.db == 'test'
    assert xref.id == '1030'


def test_XrefLink_str_complete():
    xref = XrefLink(db='test', id='1030', label='test label')
    assert xref.db == 'test'
    assert xref.id == '1030'
    assert xref.label == 'test label'


def test_EntrezLink_str_complete():
    xref = EntrezLink(db='test', id='1030', label='test label',
                      query='Test[orgn]')
    assert xref.db == 'test'
    assert xref.id == '1030'
    assert xref.label == 'test label'
    assert xref.query == 'Test[orgn]'


def test_DDBJLink_str_complete():
    xref = DDBJLink(db='test', id='1030', label='test label',
                    url='http://test.com')
    assert xref.db == 'test'
    assert xref.id == '1030'
    assert xref.label == 'test label'
    assert xref.url == 'http://test.com'


def test_ENALink_str_complete():
    xref = ENALink(db='test', id='1030', label='test label',
                   url='http://test.com')
    assert xref.db == 'test'
    assert xref.id == '1030'
    assert xref.label == 'test label'
    assert xref.url == 'http://test.com'


def test_Submission_str_complete():
    submission = Submission(submission_id='SRA12345', broker='GEO')
    submission.external_id.append(Xref(**{'db': 'test', 'id': '1030'}))
    submission.external_id.append(Xref(**{'db': 'test2', 'id': '1032'}))
    submission.external_id.append(Xref(**{'db': 'test3', 'id': '1033'}))
    submission.submitter_id.append(Xref(**{'db': 'test', 'id': 'test submitter'}))

    assert submission.submission_id == 'SRA12345'
    assert submission.broker == 'GEO'
    assert submission.external_id[0].db == 'test'
    assert submission.external_id[0].id == '1030'
    assert submission.external_id[1].db == 'test2'
    assert submission.external_id[1].id == '1032'
    assert submission.external_id[2].db == 'test3'
    assert submission.external_id[2].id == '1033'
    assert submission.submitter_id[0].db == 'test'
    assert submission.submitter_id[0].id == 'test submitter'
    assert submission.uuid == []


def test_Organization_str_complete():
    organization = Organization(
        organization_type='center', abbreviation='GEO', name='NCBI',
        email='geo@nih.gov', first_name='GEO', last_name='Curators')

    assert organization.organization_type == 'center'
    assert organization.abbreviation == 'GEO'
    assert organization.name == 'NCBI'
    assert organization.email == 'geo@nih.gov'
    assert organization.first_name == 'GEO'
    assert organization.last_name == 'Curators'

def test_Study_str_partial():
    abstract = ('A very long abstract with more than 80 characters per line. '
                'A very long abstract with more than 80 characters per line. '
                'A very long abstract with more than 80 characters per line. '
                'A very long abstract with more than 80 characters per line.')

    study = Study(
        study_id='SRA12345', GEO='GSE12345', BioProject='PRJN1234',
        pubmed='123456', title='Test title', study_type='Test type',
        abstract=abstract, center_project_name='Test center')

    study.external_id.append(Xref(**{'db': 'test', 'id': '1030'}))
    study.external_id.append(Xref(**{'db': 'test2', 'id': '1032'}))
    study.external_id.append(Xref(**{'db': 'test3', 'id': '1033'}))
    study.submitter_id.append(Xref(**{'db': 'test', 'id': 'test submitter'}))
    study.url_links.append(URLLink(**{'label': 'test', 'url': 'http://test.com'}))

    assert study.study_id == 'SRA12345'
    assert study.GEO == 'GSE12345'
    assert study.BioProject == 'PRJN1234'
    assert study.pubmed == '123456'
    assert study.title == 'Test title'
    assert study.study_type == 'Test type'
    assert study.abstract == abstract
    assert study.center_project_name == 'Test center'

    assert study.external_id[0].db == 'test'
    assert study.external_id[0].id == '1030'
    assert study.external_id[1].db == 'test2'
    assert study.external_id[1].id == '1032'
    assert study.external_id[2].db == 'test3'
    assert study.external_id[2].id == '1033'
    assert study.submitter_id[0].db == 'test'
    assert study.submitter_id[0].id == 'test submitter'
    assert study.url_links[0].label == 'test'
    assert study.url_links[0].url == 'http://test.com'
    assert study.entrez_links == []


class TestSRR3001915:
    @pytest.fixture(scope='class')
    def sraExperiment(self):
        """ Element tree of single sra experiment. """
        fname = 'tests/data/sra_SRR3001915.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return SraExperiment(root.find('EXPERIMENT_PACKAGE'))

    @pytest.fixture(scope='class')
    def bioSample(self):
        fname = 'tests/data/biosample_SAMN02981965.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return BioSampleParse(root.find('BioSample'))

    @pytest.fixture(scope='class')
    def bioProject(self):
        fname = 'tests/data/bioproject_PRJNA258012.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return BioProjectParse(root.find('DocumentSummary'))

    @pytest.fixture(scope='class')
    def Pubmed(self):
        fname = 'tests/data/pubmed_26732976.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return PubmedParse(root.find('PubmedArticle/MedlineCitation'))

    @pytest.fixture(scope='class')
    def Ncbi(self, sraExperiment, bioSample, bioProject, Pubmed):
        return Ncbi(srx=sraExperiment.srx, sra=sraExperiment.sra,
                    biosample=[bioSample.biosample,], bioproject=bioProject.bioproject,
                    pubmed=[Pubmed.pubmed,])

    def test_organization(self, Ncbi):
        assert Ncbi.sra.organization.name == 'NCBI'
        assert Ncbi.sra.organization.organization_type == 'center'
        assert Ncbi.sra.organization.abbreviation == 'GEO'
        assert Ncbi.sra.organization.name == 'NCBI'
        assert Ncbi.sra.organization.email == 'geo-group@ncbi.nlm.nih.gov'
        assert Ncbi.sra.organization.first_name == 'Geo'
        assert Ncbi.sra.organization.last_name == 'Curators'

    def test_submission(self, Ncbi):
        assert Ncbi.sra.submission.submission_id == 'SRA178685'
        assert Ncbi.sra.submission.submitter_id[0]['id'] == 'GEO: GSE60314'

    def test_study(self, Ncbi):
        assert Ncbi.sra.study.study_id == 'SRP045429'
        assert Ncbi.sra.study.BioProject == 'PRJNA258012'
        assert Ncbi.sra.study.GEO == 'GSE60314'
        assert Ncbi.sra.study.title.strip().split(' ')[0] == 'mRNA'
        assert Ncbi.sra.study.title.strip().split(' ')[-1] == 'environments'
        assert Ncbi.sra.study.study_type == 'Transcriptome Analysis'
        assert Ncbi.sra.study.abstract.strip().split(' ')[0] == 'Our'
        assert Ncbi.sra.study.abstract.strip().split(' ')[-1] == 'level.'
        assert Ncbi.sra.study.center_project_name == 'GSE60314'
        assert Ncbi.sra.study.pubmed == '26732976'

    def test_experiment(self, Ncbi):
        assert Ncbi.sra.experiment.experiment_id == 'SRX1483046'
        assert Ncbi.sra.experiment.submitter_id[0]['id'] == 'GSM1973128'
        assert Ncbi.sra.experiment.title.strip() == 'GSM1973128: 563M322; Drosophila melanogaster; RNA-Seq'
        assert Ncbi.sra.experiment.library_strategy == 'RNA-Seq'
        assert Ncbi.sra.experiment.library_source == 'TRANSCRIPTOMIC'
        assert Ncbi.sra.experiment.library_selection == 'cDNA'
        assert Ncbi.sra.experiment.library_layout == 'SINGLE'
        assert Ncbi.sra.experiment.library_construction_protocol.strip().split(' ')[0] == 'Single'
        assert Ncbi.sra.experiment.library_construction_protocol.strip().split(' ')[-1] == 'GEO_biosample_summary.xls.'
        assert Ncbi.sra.experiment.platform == 'ILLUMINA'
        assert Ncbi.sra.experiment.instrument_model == 'Illumina HiSeq 2000'
        assert Ncbi.sra.experiment.GEO_Dataset == '301973128'
        assert Ncbi.sra.experiment.attributes[0]['name'] == 'GEO Accession'
        assert Ncbi.sra.experiment.attributes[0]['value'] == 'GSM1973128'

    def test_parse_sample(self, Ncbi):
        assert Ncbi.sra.sample.sample_id == 'SRS679015'
        assert Ncbi.sra.sample.BioSample == 'SAMN02981965'
        assert Ncbi.sra.sample.GEO == 'GSM1471477'
        assert Ncbi.sra.sample.title == 'DGRP563 M_E3_2_L3'
        assert Ncbi.sra.sample.taxon_id == '7227'
        assert Ncbi.sra.sample.scientific_name == 'Drosophila melanogaster'
        assert Ncbi.sra.sample.BioProject == '258012'

        attrs = {
                'source_name': 'Whole body',
                'strain': 'DGRP-563',
                'developmental stage': 'Adult',
                'Sex': 'male',
                'tissue': 'Whole body',
                }

        for attribute in Ncbi.sra.sample.attributes:
            name = attribute['name']
            assert attrs[name] == attribute['value']

    def test_parse_pool(self, Ncbi):
        assert Ncbi.sra.pool[0] == 'SRS679015'

    def test_run(self, Ncbi):
        assert Ncbi.sra.run[0]['run_id'] == 'SRR3001915'
        assert Ncbi.sra.run[0]['samples'][0] == 'SRS679015'
        assert Ncbi.sra.run[0]['experiment_id'] == 'SRX1483046'
        assert Ncbi.sra.run[0]['nspots'] == 2001211
        assert Ncbi.sra.run[0]['nbases'] == 152092036
        assert Ncbi.sra.run[0]['nreads'] == 1
        assert Ncbi.sra.run[0]['read_count_r1'] == 2001211.0
        assert Ncbi.sra.run[0]['read_len_r1'] == 76.0

        assert Ncbi.sra.run[0]['tax_analysis']['tax_counts']['subclass'][0]['parent'] == 'Pterygota'
        assert Ncbi.sra.run[0]['tax_analysis']['tax_counts']['subclass'][0]['self_count'] == 95
        assert Ncbi.sra.run[0]['tax_analysis']['tax_counts']['subclass'][0]['total_count'] == 360608
        assert Ncbi.sra.run[0]['tax_analysis']['tax_counts']['subclass'][0]['tax_id'] == '33340'
        assert Ncbi.sra.run[0]['tax_analysis']['tax_counts']['subclass'][0]['name'] == 'Neoptera'

    def test_biosample(self, Ncbi):
        assert Ncbi.biosample[0].biosample_accn == 'SAMN02981965'
        assert Ncbi.biosample[0].sample_id == 'SRS679015'
        assert Ncbi.biosample[0].GEO == 'GSM1471477'
        assert Ncbi.biosample[0].title == 'DGRP563 M_E3_2_L3'
        assert Ncbi.biosample[0].tax_id == '7227'
        assert Ncbi.biosample[0].tax_name == 'Drosophila melanogaster'
        assert Ncbi.biosample[0].organism_name == 'Drosophila melanogaster'
        assert Ncbi.biosample[0].institute == 'Developmental Genomics, LCDB, NIDDK, NIH'
        assert Ncbi.biosample[0].access == 'public'
        assert Ncbi.biosample[0].publication_date == '2015-12-22'
        assert Ncbi.biosample[0].last_update == '2015-12-22'
        assert Ncbi.biosample[0].submission_date == '2014-08-11'
        assert Ncbi.biosample[0].contacts[0]['email'] == 'briano@helix.nih.gov'
        assert Ncbi.biosample[0].contacts[0]['first_name'] == 'Brian'
        assert Ncbi.biosample[0].contacts[0]['last_name'] == 'Oliver'
        assert Ncbi.biosample[0].models[0] == 'Generic'
        attr = {
                'source_name': 'Whole body',
                'strain': 'DGRP-563',
                'dev_stage': 'Adult',
                'sex': 'male',
                'tissue': 'Whole body',
                }
        for attribute in Ncbi.biosample[0].attributes:
            assert attribute['value'] == attr[attribute['name']]

    def test_bioproject(self, Ncbi):
        assert Ncbi.bioproject['bioproject_accn'] == 'PRJNA258012'
        assert Ncbi.bioproject['name'].strip().split(' ')[0] == 'mRNA'
        assert Ncbi.bioproject['name'].strip().split(' ')[-1] == 'environments'
        assert Ncbi.bioproject['title'].strip().split(' ')[0] == 'mRNA'
        assert Ncbi.bioproject['title'].strip().split(' ')[-1] == 'environments'
        assert Ncbi.bioproject['description'].strip().split(' ')[0] == 'Our'
        assert Ncbi.bioproject['description'].strip().split(' ')[-1] == 'level.'
        assert Ncbi.bioproject['publication'] == '26732976'
        assert Ncbi.bioproject['submission_id'] == 'SUB1475494'
        assert Ncbi.bioproject['submission_date'] == '2014-08-11'
        assert Ncbi.bioproject['last_update'] == '2016-04-20'

    def test_pubmed(self, Ncbi):
        assert Ncbi.pubmed[0].pubmed_id == '26732976'
        assert Ncbi.pubmed[0].title.strip().split(' ')[0] == 'Comparison'
        assert Ncbi.pubmed[0].title.strip().split(' ')[-1] == 'melanogaster.'
        assert Ncbi.pubmed[0].date_created == '2016-01-06'
        assert Ncbi.pubmed[0].date_completed == '2016-09-28'
        assert Ncbi.pubmed[0].date_revised == '2016-10-19'
        assert Ncbi.pubmed[0].citation == '1471-2164 17 BMC Genomics 2016'
        assert Ncbi.pubmed[0].abstract.strip().split('\n')[1].split(' ')[0] == 'We'
        assert Ncbi.pubmed[0].abstract.strip().split('\n')[2].split(' ')[1] == 'best'
        assert Ncbi.pubmed[0].authors[0]['last_name'] == 'Lin'
