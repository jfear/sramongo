"""Module to handle BioSample XML metadat."""
import re
from sramongo.xml_helpers import valid_path, parse_tree_from_dict


class BioSampleParse(object):
    def __init__(self, node):
        """Represents the biosample XML.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            A Row element from a runinfo xml as an ElemenTree element.

        """
        self.biosample = {}
        self.biosample.update(self._parse_ids(node.find('Ids')))
        self.biosample['contacts'] = self._parse_contacts(node.find('Owner/Contacts'))
        self.biosample['models'] = self._parse_models(node.find('Models'))
        self.biosample['attributes'] = self._parse_attributes(node.find('Attributes'))

        locs = {'title': ('Description/Title', 'text'),
                'tax_id': ('Description/Organism', 'taxonomy_id'),
                'tax_name': ('Description/Organism', 'taxonomy_name'),
                'organism_name': ('Description/Organism/OrganismName', 'text'),
                'description': ('Description/Comment/Paragraph', 'text'),
                'institute': ('Owner/Name', 'text'),
                'access': ('.', 'access'),
                'biosample_accn': ('.', 'accession'),
                'biosample_id': ('.', 'id'),
                'publication_date': ('.', 'publication_date'),
                'last_update': ('.', 'last_update'),
                'submission_date': ('.', 'submission_date'),
                }
        self.biosample.update(parse_tree_from_dict(node, locs))

        # Clean up
        self._drop_empty()
        self._clean_dates()

    def _drop_empty(self):
        """Scans through a dictionary and removes empy lists or None elements."""

        drop_keys = []
        for k, v in self.biosample.items():
            try:
                if (v is None) or (len(v) == 0):
                    drop_keys.append(k)
            except TypeError as err:
                # Some types (int, float) will cause an exception, if they are
                # of this type I don't really care then they are not blank so
                # skip.
                pass

        for k in drop_keys:
            del self.biosample[k]

    @valid_path
    def _parse_ids(self, node):
        d = dict()
        dbs_of_interest = ['GEO', 'SRA', 'BioSample']
        for _id in node:
            db = _id.get('db')
            primary = _id.get('is_primary')

            # Only keep useful dbs
            if (db is not None) and (db in dbs_of_interest):
                # If Primary BioSample then make biosample_id
                if (db == 'BioSample') and (primary is not None):
                    d['biosample_primary'] = _id.text
                elif (db == 'BioSample'):
                    d['biosample_secondary'] = _id.text
                elif (db == 'SRA'):
                    d['sample_id'] = _id.text
                else:
                    d[db] = _id.text

        # Try to ensure there is a biosample_id
        if ('biosample_primary' not in d) and ('biosample_secondary' in d):
            d['biosample_primary'] = d['biosample_secondary']
            del d['biosample_secondary']

        return d

    @valid_path(rettype=list)
    def _parse_contacts(self, node):
        contacts = []
        for contact in node:
            locs = {
                'email': ('.', 'email'),
                'first_name': ('Name/First', 'text'),
                'last_name': ('Name/Last', 'text'),
                }
            contacts.append(parse_tree_from_dict(contact, locs))
        return contacts

    @valid_path(rettype=list)
    def _parse_models(self, node):
        models = []
        for model in node:
            models.append(model.text)
        return models

    @valid_path(rettype=list)
    def _parse_attributes(self, node):
        attributes = []
        for attribute in node:
            name = attribute.get('harmonized_name')
            if name is None:
                name = attribute.get('attribute_name')
            attributes.append({'name': name,
                      'value': attribute.text})
        return attributes

    def _clean_dates(self):
        """Cleans up biosample dates."""
        regex = r'(\d+\-\d+\-\d+)T.*'

        self.biosample['publication_date'] = re.match(
            regex, self.biosample['publication_date']).groups()[0]

        self.biosample['last_update'] = re.match(
            regex, self.biosample['last_update']).groups()[0]

        self.biosample['submission_date'] = re.match(
            regex, self.biosample['submission_date']).groups()[0]
