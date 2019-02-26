from typing import List
from xml.etree import cElementTree as ElementTree

from sramongo.services.entrez import EfetchPackage
from .models import BioSample
from .xml_helpers import get_xml_attribute, get_xml_text


def parse_biosample(root):
    biosample = BioSample()

    biosample.accn = get_xml_attribute(root, 'BioSample', 'accession')
    biosample.id = int(get_xml_attribute(root, 'BioSample', 'id'))
    biosample.title = get_xml_text(root, 'BioSample/Description/Title')
    biosample.description = get_biosample_description(root)
    biosample.last_update = get_xml_attribute(root, 'BioSample', 'last_update')
    biosample.submission_date = get_xml_attribute(root, 'BioSample', 'submission_date')
    biosample.attributes = get_biosample_attributes(root)
    biosample.contacts = get_biosample_contacts(root)
    return biosample


def get_biosample_attributes(root):
    attributes = []
    for attribute in root.findall('BioSample/Attributes/Attribute'):
        attributes.append(
            {
                'name': attribute.attrib['harmonized_name'],
                'value': attribute.text
            }
        )
    return attributes


def get_biosample_description(root):
    text = [
        paragrpah.text
        for paragrpah in root.findall('BioSample/Description/Comment/Paragraph')
    ]
    return '\n'.join(text)


def get_biosample_contacts(root):
    contacts = []
    for contact in root.findall('BioSample/Owner/Contacts/Contact'):
        contacts.append(
            {
                'first_name': contact.find('Name/First').text,
                'last_name': contact.find('Name/Last').text,
                'email': contact.attrib['email']
            }
        )
    return contacts


def parse_biosample_set(xml: str) -> List[EfetchPackage]:
    root = ElementTree.fromstring(xml)
    for record in root.findall('BioSampleSet'):
        accn = record.find('BioSample').attrib['accession']
        record_xml = ElementTree.tostring(record).decode()
        yield EfetchPackage(accn, record_xml)