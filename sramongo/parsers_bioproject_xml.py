from .models import BioProject
from .xml_helpers import get_xml_text, get_xml_attribute


def parse_bioproject(root):
    bioproject = BioProject()
    bioproject.accn = get_xml_attribute(root, 'DocumentSummary/Project/ProjectID/ArchiveID', 'accession')
    bioproject.id = int(get_xml_attribute(root, 'DocumentSummary/Project/ProjectID/ArchiveID', 'id'))
    bioproject.name = get_xml_text(root, 'DocumentSummary/Project/ProjectDescr/Name')
    bioproject.title = get_xml_text(root, 'DocumentSummary/Project/ProjectDescr/Title')
    bioproject.description = get_xml_text(root, 'DocumentSummary/Project/ProjectDescr/Description')
    # TODO make dates using datetime
    bioproject.last_update = get_xml_attribute(root, 'DocumentSummary/Submission', 'last_update')
    bioproject.submission_date = get_xml_attribute(root, 'DocumentSummary/Submission', 'submitted')
    return bioproject
