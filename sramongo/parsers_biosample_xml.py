from typing import List
from xml.etree import cElementTree as ElementTree

from sramongo.services.entrez import EfetchPackage, EsummaryResult
from .models import BioSample
from .utils import make_number, date_parse
from .xml_helpers import get_xml_text, xml_to_root


def parse_biosample(root):
    biosample = BioSample()
    attributes = root.attrib
    biosample.accn = attributes['accession']
    biosample.id = make_number(attributes['id'], int)
    biosample.last_update = date_parse(attributes['last_update'])
    biosample.submission_date = date_parse(attributes['submission_date'])

    biosample.title = get_xml_text(root, 'Description/Title')
    biosample.description = get_biosample_description(root)
    biosample.attributes = get_biosample_attributes(root)
    biosample.contacts = get_biosample_contacts(root)

    return biosample


def get_biosample_attributes(root):
    attributes = []
    for attribute in root.findall('Attributes/Attribute'):
        name = attribute.attrib.get('harmonized_name', attribute.attrib.get('attribute_name', ''))
        value = attribute.text
        attributes.append({'name': name, 'value': value})
    return attributes


def get_biosample_description(root):
    text = [
        paragrpah.text
        for paragrpah in root.findall('Description/Comment/Paragraph')
    ]
    return '\n'.join(text)


def get_biosample_contacts(root):
    contacts = []
    for contact in root.findall('Owner/Contacts/Contact'):
        _contact = {}
        if contact.find('Name/First'):
            _contact['first_name'] = contact.find('Name/First').text

        if contact.find('Name/Last'):
            _contact['last_name'] = contact.find('Name/Last').text

        if contact.attrib.get('email', ''):
            _contact['email'] = contact.attrib['email']

        if _contact:
            contacts.append(_contact)

    return contacts


def parse_biosample_efetch_result(xml: str) -> List[EfetchPackage]:
    root = ElementTree.fromstring(xml)
    for record in root.findall('BioSample'):
        accn = record.attrib['accession']
        record_xml = ElementTree.tostring(record).decode()
        yield EfetchPackage(accn, record_xml)


def parse_biosample_esummary_result(xml: str) -> List[EsummaryResult]:
    root = xml_to_root(xml)
    for doc in root.findall('DocumentSummary'):
        uid = doc.find('Id').text
        accn = doc.find("Accession").text
        create_date = date_parse(doc.find("PublicationDate").text)
        update_date = date_parse(doc.find("ModificationDate").text)
        yield EsummaryResult(uid, accn, create_date, update_date)
