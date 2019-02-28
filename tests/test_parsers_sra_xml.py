from sramongo import parsers_sra_xml, models


def test_add_library_layout(sra_xml_root):
    root = sra_xml_root
    sra = models.SraDocument()
    parsers_sra_xml.add_library_layout(root, sra)
    assert sra.library_layout == 'SINGLE'


def test_add_platform_information(sra_xml_root):
    root = sra_xml_root
    sra = models.SraDocument()
    parsers_sra_xml.add_platform_information(root, sra)
    assert sra.platform == 'ILLUMINA'
    assert sra.instrument_model == 'Illumina Genome Analyzer IIx'


def test_parse_sra_study(sra_xml_root):
    root = sra_xml_root
    study = parsers_sra_xml.parse_sra_study(root)
    assert study.accn == 'SRP056660'
    assert study.center_name == 'GEO'


def test_add_study_external_links(sra_xml_root2):
    root = sra_xml_root2
    study = models.Study()
    parsers_sra_xml.add_study_external_links(root, study)
    assert study.bioproject == 'PRJNA258012'
    assert study.geo == 'GSE60314'


def test_add_study_pubmed(sra_xml_root):
    root = sra_xml_root
    study = models.Study()
    parsers_sra_xml.add_study_pubmed(root, study)
    assert study.pubmed[0] == '26335107'


def test_parse_sra_organization(sra_xml_root2):
    root = sra_xml_root2
    organization = parsers_sra_xml.parse_sra_organization(root)
    assert organization.organization_type == 'center'


def test_add_sample_attributes(sra_xml_root2):
    root = sra_xml_root2
    sample = models.Sample()
    parsers_sra_xml.add_sample_attributes(root, sample)
    assert sample.attributes[0]['name'] == 'source_name'
    assert sample.attributes[0]['value'] == 'Whole body'


def test_add_sample_external_links(sra_xml_root2):
    root = sra_xml_root2
    sample = models.Sample()
    parsers_sra_xml.add_sample_external_links(root, sample)
    assert sample.biosample == 'SAMN02981965'
    assert sample.geo == 'GSM1471477'


def test_parse_sra_run_single_end(sra_xml_root):
    root = sra_xml_root
    runs = parsers_sra_xml.parse_sra_run(root)
    run = runs[0]
    assert len(runs) == 1
    assert run.srr == 'SRR1945105'
    assert run.nspots == 7_018_100
    assert run.nbases == 533_375_600
    assert run.nreads == 1
    assert run.read_count_r1 == 7_018_100
    assert run.read_len_r1 == 76


def test_parse_sra_run_paired_end(sra_xml_root_PE):
    root = sra_xml_root_PE
    runs = parsers_sra_xml.parse_sra_run(root)
    run = runs[0]
    assert run.nreads == 2
    assert run.read_count_r1 == 12_314_272
    assert run.read_len_r1 == 99
    assert run.read_count_r2 == 12_314_272
    assert run.read_len_r2 == 99


def test_parse_sra_experiment(sra_xml_root):
    root = sra_xml_root
    sra = parsers_sra_xml.parse_sra_experiment(root)
    assert sra.srx == 'SRX971855'
    assert sra.title == 'GSM1646282: AH_dsxNullM_3; Drosophila melanogaster; RNA-Seq'
    assert sra.design == None
    assert sra.library_name == ''
    assert sra.library_strategy == 'RNA-Seq'
    assert sra.library_source == 'TRANSCRIPTOMIC'
    assert sra.library_selection == 'cDNA'
    assert sra.library_construction_protocol[:5] == 'Total'
    assert sra.library_layout == 'SINGLE'
    assert sra.platform == 'ILLUMINA'
    assert sra.instrument_model == 'Illumina Genome Analyzer IIx'
    assert sra.study.accn == 'SRP056660'
    assert sra.organization.organization_type == 'center'
    assert sra.runs[0].srr == 'SRR1945105'
    assert sra.sample.attributes[0]['name'] == 'source_name'
    assert sra.sample.attributes[0]['value'] == 'Adult Head Tissue'
