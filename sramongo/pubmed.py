"""Module to parse pubmed XML."""
from sramongo.xml_helpers import valid_path, parse_tree_from_dict

class Pubmed(object):
    def __init__(self, node):
        locs = {
                'pubmed': ('PMID', 'text'),
                'title': ('Article/ArticleTitle', 'text'),
                }
        self.__dict__.update(parse_tree_from_dict(node, locs))
        self.date_created = self._parse_dates(node.find('DateCreated'))
        self.date_completed = self._parse_dates(node.find('DateCompleted'))
        self.date_revised = self._parse_dates(node.find('DateRevised'))
        self.citation = self._parse_citation(node.find('Article/Journal'))
        self.abstract = self._parse_abstract(node.find('Article/Abstract'))
        self.authors = self._parse_authors(node.find('Article/AuthorList'))

    @valid_path
    def _parse_dates(self, node):
        year = node.find('Year').text
        month = node.find('Month').text
        day = node.find('Day').text
        return '-'.join([year, month, day])

    @valid_path
    def _parse_citation(self, node):
        issue = node.find('ISSN').text
        volume = node.find('JournalIssue/Volume').text
        year = node.find('JournalIssue/PubDate/Year').text
        journal = node.find('ISOAbbreviation').text
        return '{v} {j} {y}'.format(i=issue, y=year, j=journal, v=volume)

    @valid_path
    def _parse_abstract(self, node):
        parts = []
        for part in node:
            parts.append(part.text)
        return '\n'.join(parts)

    @valid_path
    def _parse_authors(self, node):
        authors =[]
        for author in node:
            d = {
                'first_name': author.find('ForeName').text,
                'last_name': author.find('LastName').text,
                'affiliation': author.find('AffiliationInfo/Affiliation').text,
                }
            authors.append(d)
        return authors
