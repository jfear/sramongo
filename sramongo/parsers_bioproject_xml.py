from typing import List
from xml.etree import cElementTree as ElementTree

from .services.entrez import EfetchPackage, EsummaryResult
from .models import BioProject
from .xml_helpers import get_xml_text, get_xml_attribute, xml_to_root
from .utils import make_number, date_parse


def parse_bioproject(root):
    bioproject = BioProject()
    bioproject.accn = get_xml_attribute(root, 'Project/ProjectID/ArchiveID', 'accession')
    bioproject.id = make_number(get_xml_attribute(root, 'Project/ProjectID/ArchiveID', 'id'), int)
    bioproject.name = get_xml_text(root, 'Project/ProjectDescr/Name')
    bioproject.title = get_xml_text(root, 'Project/ProjectDescr/Title')
    bioproject.description = get_xml_text(root, 'Project/ProjectDescr/Description')
    bioproject.last_update = date_parse(get_xml_attribute(root, 'Submission', 'last_update'))
    bioproject.submission_date = date_parse(get_xml_attribute(root, 'Submission', 'submitted'))
    return bioproject


def parse_bioproject_efetch_result(xml: str) -> List[EfetchPackage]:
    root = ElementTree.fromstring(xml)
    for record in root.findall('DocumentSummary'):
        accn = record.find('Project/ProjectID/ArchiveID').attrib['accession']
        record_xml = ElementTree.tostring(record).decode()
        yield EfetchPackage(accn, record_xml)


def parse_bioproject_esummary_result(xml: str) -> List[EsummaryResult]:
    root = xml_to_root(xml)
    for doc in root.findall('*DocumentSummary'):
        uid = doc.find('Project_Id').text
        accn = doc.find("Project_Acc").text
        create_date = date_parse(doc.find("Registration_Date").text)
        yield EsummaryResult(uid, accn, create_date, '')
