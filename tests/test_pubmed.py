#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.pubmed import PubmedParse


class TestPMID26732976:
    @pytest.fixture(scope='class')
    def pub_etree(self):
        """ Element tree of single pubmed id. """
        fname = 'tests/data/pubmed_26732976.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('PubmedArticle/MedlineCitation')

    @pytest.fixture(autouse=True, scope='class')
    def pubTree(self, pub_etree):
        return PubmedParse(pub_etree)

    def test_parse(self, pubTree):
        assert pubTree.pubmed['pubmed_id'] == '26732976'
        assert pubTree.pubmed['title'].strip().split(' ')[0] == 'Comparison'
        assert pubTree.pubmed['title'].strip().split(' ')[-1] == 'melanogaster.'
        assert pubTree.pubmed['date_created'] == '2016-01-06'
        assert pubTree.pubmed['date_completed'] == '2016-09-28'
        assert pubTree.pubmed['date_revised'] == '2016-10-19'
        assert pubTree.pubmed['citation'] == '1471-2164 17 BMC Genomics 2016'
        assert pubTree.pubmed['abstract'].strip().split('\n')[1].split(' ')[0] == 'We'
        assert pubTree.pubmed['abstract'].strip().split('\n')[2].split(' ')[1] == 'best'
        assert pubTree.pubmed['authors'][0]['last_name'] == 'Lin'
