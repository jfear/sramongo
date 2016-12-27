#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.pubmed import Pubmed


class TestPMID26732976:
    @pytest.fixture(scope='class')
    def pub_etree(self):
        """ Element tree of single pubmed id. """
        fname = 'data/pubmed_26732976.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('PubmedArticle/MedlineCitation')

    @pytest.fixture(autouse=True, scope='class')
    def pubTree(self, pub_etree):
        return Pubmed(pub_etree)

    def test_parse(self, pubTree):
        assert pubTree.pubmed == '26732976'
        assert pubTree.title.strip().split(' ')[0] == 'Comparison'
        assert pubTree.title.strip().split(' ')[-1] == 'melanogaster.'
        assert pubTree.date_created == '2016-01-06'
        assert pubTree.date_completed == '2016-09-28'
        assert pubTree.date_revised == '2016-10-19'
        assert pubTree.citation == '17 BMC Genomics 2016'
        assert pubTree.abstract.strip().split('\n')[1].split(' ')[0] == 'We'
        assert pubTree.abstract.strip().split('\n')[2].split(' ')[1] == 'best'
        assert pubTree.authors[0]['last_name'] == 'Lin'
