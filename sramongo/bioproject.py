"""Module to handle BioProject XML metadat."""
import re
from sramongo.xml_helpers import valid_path, parse_tree_from_dict


class BioProjectParse(object):
    def __init__(self, node):
        """Represents the bioproject XML.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            A Row element from a runinfo xml as an ElemenTree element.

        """
        locs = {
                'bioproject_accn': ('Project/ProjectID/ArchiveID', 'accession'),
                'bioproject_id': ('Project/ProjectID/ArchiveID', 'id'),
                'name': ('Project/ProjectDescr/Name', 'text'),
                'title': ('Project/ProjectDescr/Title', 'text'),
                'description': ('Project/ProjectDescr/Description', 'text'),
                'publication': ('Project/ProjectDescr/Publication', 'id'),
                'publication_date': ('Project/ProjectDescr/Publication', 'date'),
                'submission_id': ('Submission', 'submission_id'),
                'last_update': ('Submission', 'last_update'),
                'submission_date': ('Submission', 'submitted'),
                }

        self.bioproject = {}
        self.bioproject.update(parse_tree_from_dict(node, locs))
        self.bioproject['external_id'] = self._parse_links(node.find('ExternalLink'))

        # Clean up
        self._drop_empty()
        self._clean_dates()

    def _drop_empty(self):
        """Scans through a dictionary and removes empy lists or None elements."""

        drop_keys = []
        for k, v in self.bioproject.items():
            try:
                if (v is None) or (len(v) == 0):
                    drop_keys.append(k)
            except TypeError as err:
                # Some types (int, float) will cause an exception, if they are
                # of this type I don't really care then they are not blank so
                # skip.
                pass

        for k in drop_keys:
            del self.bioproject[k]

    def _clean_dates(self):
        """Cleans up bioproject dates.

        Sometimes dates are stored as 'YYYY-MM-DDTHH:MM:SSZ', strip out the
        time information if it is there.
        """
        regex = r'(\d+\-\d+\-\d+)T.*'

        try:
            self.bioproject['publication_date'] = re.match(
                regex, self.bioproject['publication_date']).groups()[0]
        except (AttributeError, KeyError):
            pass

        try:
            self.bioproject['last_update'] = re.match(
                regex, self.bioproject['last_update']).groups()[0]
        except (AttributeError, KeyError):
            pass

        try:
            self.bioproject['submission_date'] = re.match(
                regex, self.bioproject['submission_date']).groups()[0]
        except (AttributeError, KeyError):
            pass

    @valid_path(rettype=list)
    def _parse_links(self, node):
        links = []
        for xref in node:
            d = {}
            d[xref.get('db')] = xref.find('ID').text
            links.append(d)
