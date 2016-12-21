#!/usr/bin/env python
""" Module to handle SRA XML metadata. """

from xml.etree import ElementTree
from collections import defaultdict

from mongoengine import connect

from sramongo import db_schema as dbs

def valid_path(func):
    def new_func(*args, **kwargs):
        # If the current path is present
        if args[0]:
            return func(*args, **kwargs)
        else:
            print('Not valid path.')
            return {}

    return new_func


class AmbiguousElementException(Exception):
    pass


class AmbiguousElementException(Exception):
    pass


class SraExperiment(object):
    """ Class to parse SRA experiments and add them to mongoDB.

    Parameters
    ----------
    node: xml.etree.ElementTree.ElementTree.Element
        Experiment level node as an ElemenTree element.

    db_name: str
        Name of the database to use. NOTE: monoDB server needs to be running.

    """
    def __init__(self, node, db_name='sraDB'):
        # connect to database
        #self.db = MongoClient(db_name)

        # parse different sections for SRA xml
        self.organization = self.parse_organization(node.find('Organization'))
        self.submission = self.parse_submission(node.find('SUBMISSION'))
        self.study = self.parse_study(node.find('STUDY'))
        self.experiment = self.parse_experiment(node.find('EXPERIMENT'))
        #self.sample = self.parse_experiment(node.find('SAMPLE'))
        #self.run = self.parse_run(node.find('RUN_SET'))
        #self.pool = self.parse_experiment(node.find('Pool'))


    def parse_submission(self, node):
        d = dict()
        d.update(self._parse_ids(node.find('IDENTIFIERS'), 'submission'))
        return d

    def parse_organization(self, node):
        d = dict()
        locs = {
                'type': ('.', 'type'),
                'abbreviation': ('Name', 'abbr'),
                'name': ('Name', 'text'),
                'email': ('Contact', 'email'),
                'first_name': ('Contact/Name/First', 'text'),
                'last_name': ('Contact/Name/Last', 'text'),
                }
        d.update(self._parse_relationships(node, locs))

        return d

    def parse_study(self, node):
        d = dict()
        d.update(self._parse_ids(node.find('IDENTIFIERS'), 'study'))
        d.update(self._parse_links(node.find('STUDY_LINKS'), 'study'))
        locs = {
                'study_title': ('DESCRIPTOR/STUDY_TITLE', 'text'),
                'study_abstract': ('DESCRIPTOR/STUDY_ABSTRACT', 'text'),
                'center_project_name': ('DESCRIPTOR/CENTER_PROJECT_NAME', 'text')
                }
        d.update(self._parse_relationships(node, locs))
        d.update(self._parse_study_type(node.find('DESCRIPTOR/STUDY_TYPE')))

        return d

    def parse_experiment(self, node):
        d = dict()

        d.update(self._parse_ids(node.find('IDENTIFIERS'), 'experiment'))
        d.update(self._parse_links(node.find('EXPERIMENT_LINKS'), 'experiment'))
        d.update(self._parse_attributes(node.find('EXPERIMENT_ATTRIBUTES'), 'experiment'))

        locs = {
                'title': ('TITLE', 'text'),
                'study': ('STUDY_REF/IDENTIFIERS/PRIMARY_ID', 'text'),
                'design': ('DESIGN/DESIGN_DESCRIPTION', 'text'),
                'library_name': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME', 'text'),
                'library_strategy': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY', 'text'),
                'library_source': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE', 'text'),
                'library_selection': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION', 'text'),
                'library_layout': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT', 'child', 'tag'),
                'library_layout_orientation': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED', 'ORIENTATION'),
                'library_layout_length': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED', 'NOMINAL_LENGTH'),
                'library_layout_sdev': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED', 'NOMINAL_SDEV'),
                'pooling_stategy': ('DESIGN/LIBRARY_DESCRIPTOR/POOLING_STRATEGY', 'text'),
                'library_construction_protocol': ('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL', 'text'),
                'platform': ('PLATFORM', 'child', 'tag'),
                'instrument_model': ('PLATFORM/*/INSTRUMENT_MODEL', 'text'),
                }

        d.update(self._parse_relationships(node, locs))

        return d

    def parse_run(self, node):
        d = dict()
        return d

    def parse_sample(self, node):
        d = dict()

        return d

    @valid_path
    def _parse_study_type(self, node):
        """ Processes study types.

        node: xml.etree.ElementTree.ElementTree.element
            'STUDY/DESCRIPTOR/STUDY_TYPE'
        """
        d = dict()
        if node.get('existing_study_type'):
            d['study_type'] = node.get('existing_study_type')
        elif node.get('new_study_type'):
            d['study_type'] = node.get('new_study_type')

        return d

    @valid_path
    def _parse_relationships(self, node, locs):
        """ Processes key locations.

        node: xml.etree.ElementTree.ElementTree.element
            Current node.
        locs: dict
            A dictionary mapping key to a tuple. The tuple can either be 2 or 3
            elements long. The first element maps to the location in the
            current node. The second element given a processing hint. Possible
            values are:

                * 'text': assumes the wanted is the text element of the path.
                * 'child': assumes that the child of the given path is wanted.
                * str: Any other string will be treated as an attribute lookup of the path.

            If 'child' is given, then a third element needs to be given
            indicating the type of processing. Possible values are:

                * 'text': assumes the wanted is the text element of the path.
                * 'tag': assumes the wanted is the class tag of the path.
                * str: Any other string will be treated as an attribute lookup of the path.

        """
        d = dict()
        for n, l in locs.items():
            try:
                if l[1] == 'text':
                    d[n] = node.find(l[0]).text
                elif l[1] == 'child':
                    child = node.find(l[0]).getchildren()

                    if len(child) > 1:
                        raise AmbiguousElementException('There are too many elements')
                    elif l[2] == 'text':
                        d[n] = child[0].text
                    elif l[2] == 'tag':
                        d[n] = child[0].tag
                else:
                    d[n] = node.find(l[0]).get(l[1])
            except:
                pass

        return d

    @valid_path
    def _parse_ids(self, node, namespace):
        """ Helper function to parse IDENTIFIERS section.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is 'IDENTIFIERS'.
        namespace: str
            Name to use for the primary ID.

        """
        d = defaultdict(list)

        for _id in node:
            if _id.tag == 'PRIMARY_ID':
                d[namespace + '_id'] = _id.text

            elif _id.tag == 'EXTERNAL_ID':
                d['external_id'].append({'db': _id.get('namespace'), 'id': _id.text})

            elif _id.tag == 'SECONDARY_ID':
                d['secondary_id'].append({'db': _id.get('namespace'), 'id': _id.text})

            elif _id.tag == 'SUBMITTER_ID':
                d['submitter_id'].append({'db': _id.get('namespace'), 'id': _id.text})

        return d

    @valid_path
    def _parse_link_parts(self, node):
        """ Parse the different parts of a link.

        In SRA a link can have a {label, db, id, query, url}. This function
        iterates over the different attributes and pulls out the bits of
        information.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is '{URL,XREF,ENTRE,DDBJ,ENA}_LINK'.

        """
        d = dict()
        for child in node:
            d[child.tag.lower()] = child.text
        return d

    @valid_path
    def _parse_links(self, node, namespace):
        """ Helper script to parse various types of links in the SRA.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is '*_LINKS'.
        namespace: str
            Name to use for the primary ID.

        """
        d = defaultdict(list)

        for url in node.findall('*/URL_LINK'):
            d[namespace + '_url_link'].append(self._parse_link_parts(url))
        for xref in node.findall('*/XREF_LINK'):
            d[namespace + '_xref_link'].append(self._parse_link_parts(xref))
        for entrez in node.findall('*/ENTREZ_LINK'):
            d[namespace + '_entrez_link'].append(self._parse_link_parts(entrez))
        for ddbj in node.findall('*/DDBJ_LINK'):
            d[namespace + '_ddbj_link'].append(self._parse_link_parts(ddbj))
        for ena in node.findall('*/ENA_LINK'):
            d[namespace + '_ena_link'].append(self._parse_link_parts(ena))

        return d

    @valid_path
    def _parse_attributes(self, node, namespace):
        """ Helper script to parse various types of links in the SRA.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is 'EXPERIMENT_ATTRIBUTES'.
        namespace: str
            Name to use for the primary ID.

        """
        n = namespace + '_attribute'
        d = {n: {}}

        for attribute in node.getchildren():
            d[n][attribute.find('TAG').text] = attribute.find('VALUE').text

        return d


def parse_sra_xml(xml):
    tree = ElementTree.parse(xml)
    root = tree.getroot()

    for package in root:
        sraExp = SraExperiment(package)


if __name__ == '__main__':
    pass
