"""Parse parts of the SRA XML into JSON."""
from typing import List
from xml.etree import cElementTree as ElementTree

from sramongo.services.entrez import EfetchPackage
from sramongo.xml_helpers import get_xml_text, get_xml_attribute
from .models import SraDocument, Study, Organization, Sample, Run


def parse_sra_experiment(root):
    sra = SraDocument()
    # Experiment information
    sra.srx = get_xml_text(root, 'EXPERIMENT/IDENTIFIERS/PRIMARY_ID')
    sra.title = get_xml_text(root, 'EXPERIMENT/TITLE')
    sra.design = get_xml_text(root, 'EXPERIMENT/DESIGN/DESIGN_DESCRIPTION')
    sra.library_name = get_xml_text(root, 'EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME')
    sra.library_strategy = get_xml_text(root, 'EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY')
    sra.library_source = get_xml_text(root, 'EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE')
    sra.library_selection = get_xml_text(root, 'EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION')
    sra.library_construction_protocol = get_xml_text(root, 'EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL')
    add_library_layout(root, sra)
    add_platform_information(root, sra)

    # Embedded Documents
    sra.study = parse_sra_study(root)
    sra.organization = parse_sra_organization(root)
    sra.sample = parse_sra_sample(root)
    sra.runs = parse_sra_run(root)
    return sra


def parse_sra_study(root):
    study = Study()
    study.accn = get_xml_text(root, 'STUDY/IDENTIFIERS/PRIMARY_ID')
    study.title = get_xml_text(root, 'STUDY/DESCRIPTOR/STUDY_TITLE')
    study.abstract = get_xml_text(root, 'STUDY/DESCRIPTOR/STUDY_ABSTRACT')
    study.center_name = get_xml_attribute(root, 'STUDY', 'center_name')
    study.center_project_name = get_xml_text(root, 'STUDY/DESCRIPTOR/CENTER_PROJECT_NAME')

    description = get_xml_text(root, 'STUDY/DESCRIPTOR/STUDY_DESCRIPTION')
    if description is not None and description != study.abstract:
        study.description = get_xml_text(root, '')
    return study


def parse_sra_organization(root):
    organization = Organization()
    organization.organization_type = get_xml_attribute(root, 'Organization', 'type')
    organization.abbreviation = get_xml_attribute(root, 'Organization/Name', 'abbr')
    organization.name = get_xml_text(root, 'Organization/Name')
    organization.email = get_xml_attribute(root, 'Organization/Contact', 'email')
    organization.first_name = get_xml_text(root, 'Organization/Contact/NAME/FIRST')
    organization.last_name = get_xml_text(root, 'Organization/Contact/NAME/FIRST')
    return organization


def parse_sra_sample(root):
    sample = Sample()
    sample.accn = get_xml_text(root, 'SAMPLE/IDENTIFIERS/PRIMARY_ID')
    sample.title = get_xml_text(root, 'SAMPLE/TITLE')
    sample.taxon_id = get_xml_text(root, 'SAMPLE/SAMPLE_NAME/TAXON_ID')
    sample.scientific_name = get_xml_text(root, 'SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME')
    sample.common_name = get_xml_text(root, 'SAMPLE/SAMPLE_NAME/COMMON_NAME')
    add_sample_attributes(root, sample)
    add_sample_external_links(root, sample)
    return sample


def parse_sra_run(root):
    runs = []
    for run in root.findall('RUN_SET/RUN'):
        sra_run = Run()
        sra_run.srr = get_xml_text(run, 'IDENTIFIERS/PRIMARY_ID')
        sra_run.nspots = int(run.attrib['total_spots'])
        sra_run.nbases = int(run.attrib['total_bases'])
        sra_run.nreads = int(get_xml_attribute(run, 'Statistics', 'nreads'))

        for read in run.findall('Statistics/Read'):
            idx = read.attrib['index']
            count = read.attrib['count']
            length = read.attrib['average']

            if idx == '0':
                sra_run.read_count_r1 = int(count)
                sra_run.read_len_r1 = int(length)
            elif idx == '1':
                sra_run.read_count_r2 = int(count)
                sra_run.read_len_r2 = int(length)
        runs.append(sra_run)
    return runs


def add_library_layout(root, sra):
    if root.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/SINGLE') is not None:
        sra.library_layout = 'SINGLE'
    else:
        sra.library_layout = 'PARIED'
        sra.library_layout_length = int(
            get_xml_attribute(root, 'EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED', 'NOMINAL_LENGTH'))
        sra.library_layout_sdev = float(
            get_xml_attribute(root, 'EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED', 'NOMINAL_SDEV'))
    return sra


def add_platform_information(root, sra):
    sra.platform = root.find('EXPERIMENT/PLATFORM').getchildren()[0].tag
    sra.instrument_model = get_xml_text(root, f'EXPERIMENT/PLATFORM/{sra.platform}/INSTRUMENT_MODEL')
    return sra


def add_sample_attributes(root, sample):
    for attribute in root.findall('SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE'):
        tag = attribute.find('TAG').text
        value = attribute.find('VALUE').text
        sample.attributes.append({'name': tag, 'value': value})


def add_sample_external_links(root, sample):
    for external_id in root.findall('SAMPLE/IDENTIFIERS/EXTERNAL_ID'):
        namespace = external_id.attrib.get('namespace')
        if namespace == 'BioSample':
            sample.biosample = external_id.text
        elif namespace == 'GEO':
            sample.geo = external_id.text

# TODO add parser somewhere to run table.
# # NOTE: Additional Fields not in the SRA XML but in summary table
# run.release_date = DateTimeField()
# run.load_date = DateTimeField()
# run.size_MB = IntField()

def parse_sra_experiment_set(xml: str) -> List[EfetchPackage]:
    root = ElementTree.fromstring(xml)
    for experiment in root.findall('EXPERIMENT_PACKAGE'):
        srx = experiment.find('EXPERIMENT').attrib['accession']
        experiment_xml = ElementTree.tostring(experiment).decode()
        yield EfetchPackage(srx, experiment_xml)