from .models import Pubmed
from .xml_helpers import get_xml_text


def parse_pubmed(root):
    pubmed = Pubmed()
    pubmed.accn = ''
    pubmed.title = ''
    pubmed.abstract = ''
    pubmed.authors = ''
    pubmed.citation = ''
    pubmed.date_created = ''
    pubmed.date_completed = ''
    pubmed.date_revised = ''

    return pubmed
