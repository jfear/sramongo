from typing import List
from xml.etree import cElementTree as ElementTree

from sramongo.services.entrez import EfetchPackage, EsummaryResult
from .models import Pubmed
from .utils import make_number, date_parse
from .xml_helpers import get_xml_text, xml_to_root


def parse_pubmed(root):
    pubmed = Pubmed()
    pubmed.accn = make_number(get_xml_text(root, 'MedlineCitation/PMID'), int)
    pubmed.authors = get_authors(root)
    pubmed.title = get_xml_text(root, 'MedlineCitation/Article/ArticleTitle')
    pubmed.abstract = get_abstract(root)
    pubmed.journal = get_xml_text(root, 'MedlineCitation/Article/Journal/ISOAbbreviation')

    pubmed.volume = get_xml_text(root, 'MedlineCitation/Article/Journal/JournalIssue/Volume')
    pubmed.page = get_xml_text(root, 'MedlineCitation/Article/Pagination/MedlinePgn')
    pubmed.year = get_xml_text(root, 'MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')

    date_created = get_date(root, 'MedlineCitation/DateCreated')
    if date_created:
        pubmed.date_created = date_created

    date_completed = get_date(root, 'MedlineCitation/DateCompleted')
    if date_completed:
        pubmed.date_completed = date_completed

    date_revised = get_date(root, 'MedlineCitation/DateRevised')
    if date_revised:
        pubmed.date_revised = date_revised

    pubmed.citation = create_citation(pubmed)
    return pubmed


def get_abstract(root):
    text = [
        paragraph.text
        for paragraph in root.findall('MedlineCitation/Article/Abstract/AbstractText')
    ]

    if len(text) == 0 or text[0] is None:
        return ''

    return '\n'.join(text)


def get_authors(root):
    authors = []
    for author in root.findall('MedlineCitation/Article/AuthorList/Author'):
        authors.append(
            {
                'first_name': get_xml_text(author, 'ForeName'),
                'last_name': get_xml_text(author, 'LastName'),
                'initials': get_xml_text(author, 'Initials'),
                'affiliation': get_xml_text(author, 'AffiliationInfo/Affiliation'),
            }
        )
    return authors


def get_date(root, path):
    curr_path = root.find(path)
    if curr_path is None:
        return ''

    year = curr_path.find('Year').text
    month = curr_path.find('Month').text
    day = curr_path.find('Day').text
    return date_parse(f'{year}-{month}-{day}')


def create_citation(pubmed):
    return (f'{pubmed.authors[0]["last_name"]}, {pubmed.authors[0]["first_name"]}, '
            f'et al. {pubmed.journal} {pubmed.volume} ({pubmed.year}): {pubmed.page}.')


def parse_pubmed_efetch_result(xml: str) -> List[EfetchPackage]:
    root = ElementTree.fromstring(xml)
    for record in root.findall('PubmedArticle'):
        accn = record.find('MedlineCitation/PMID').text
        record_xml = ElementTree.tostring(record).decode()
        yield EfetchPackage(accn, record_xml)


def parse_pubmed_esummary_result(xml: str) -> List[EsummaryResult]:
    root = xml_to_root(xml)
    for doc in root.findall('DocumentSummary'):
        uid = doc.find('Id').text
        accn = f'PMID{uid}'
        create_date = date_parse(doc.find("PubDate").text)
        yield EsummaryResult(uid, accn, create_date, '')
