"""Module to handle BioProject XML metadat."""
import re
from sramongo.xml_helpers import valid_path, parse_tree_from_dict


class BioProject(object):
    def __init__(self, node):
        """Represents the biosample XML.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            A Row element from a runinfo xml as an ElemenTree element.

        """
        locs = {
                'BioProject': ('Project/ProjectID/ArchiveID', 'accession'),
                'name': ('Project/ProjectDescr/Name', 'text'),
                'title': ('Project/ProjectDescr/Title', 'text'),
                'description': ('Project/ProjectDescr/Description', 'text'),
                'publication': ('Project/ProjectDescr/Publication', 'id'),
                'submission_id': ('Submission', 'submission_id'),
                'last_update': ('Submission', 'last_update'),
                'submission_date': ('Submission', 'submitted'),
                }

        self.__dict__.update(parse_tree_from_dict(node, locs))
        self.xref = self._parse_links(node.find('ExternalLink'))

    @valid_path
    def _parse_links(self, node):
        links = []
        for xref in node:
            d = {}
            d[xref.get('db')] = xref.find('ID').text
            links.append(d)
