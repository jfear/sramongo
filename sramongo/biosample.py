"""Module to handle BioSample XML metadat."""
import re
from sramongo.xml_helpers import valid_path, parse_tree_from_dict


class BioSample(object):
    def __init__(self, node):
        """Represents the biosample XML.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            A Row element from a runinfo xml as an ElemenTree element.

        """
        self.__dict__.update(self._parse_ids(node.find('Ids')))
        self.contacts = self._parse_contacts(node.find('Owner/Contacts'))
        self.models = self._parse_models(node.find('Models'))
        self.attributes = self._parse_attributes(node.find('Attributes'))

        locs = {'title': ('Description/Title', 'text'),
                'tax_id': ('Description/Organism', 'taxonomy_id'),
                'tax_name': ('Description/Organism', 'taxonomy_name'),
                'organism_name': ('Description/Organism/OrganismName', 'text'),
                'institute': ('Owner/Name', 'text'),
                'access': ('.', 'access'),
                'publication_date': ('.', 'publication_date'),
                'last_update': ('.', 'last_update'),
                'submission_date': ('.', 'submission_date'),
                }
        self.__dict__.update(parse_tree_from_dict(node, locs))

        self._clean_dates()

    @valid_path
    def _parse_ids(self, node):
        d = dict()
        for id in node:
            d[id.get('db')] = id.text
        return d

    @valid_path
    def _parse_contacts(self, node):
        contacts = []
        for contact in node:
            d = {
                'email': contact.find('.').get('email'),
                'first_name': contact.find('Name/First').text,
                'last_name': contact.find('Name/Last').text,
                }
            contacts.append(d)
        return contacts

    @valid_path
    def _parse_models(self, node):
        models = []
        for model in node:
            models.append(model.text)
        return models

    @valid_path
    def _parse_attributes(self, node):
        d = {}
        for attribute in node:
            try:
                d[attribute.get('harmonized_name')] = attribute.text
            except:
                d[attribute.get('attribute_name')] = attribute.text
        return d

    def _clean_dates(self):
        """Cleans up biosample dates."""
        regex = r'(\d+\-\d+\-\d+)T.*'

        self.publication_date = re.match(
                regex, self.publication_date).groups()[0]

        self.last_update = re.match(
                regex, self.last_update).groups()[0]

        self.submission_date = re.match(
                regex, self.submission_date).groups()[0]
