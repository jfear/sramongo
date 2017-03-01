#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.bioproject import BioProjectParse


class TestPRJNA258012:
    @pytest.fixture(scope='class')
    def bio_etree(self):
        """ Element tree of single biosample. """
        fname = 'tests/data/bioproject_PRJNA258012.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('DocumentSummary')

    @pytest.fixture(autouse=True, scope='class')
    def bioTree(self, bio_etree):
        return BioProjectParse(bio_etree)

    def test_parse(self, bioTree):
        assert bioTree.bioproject['bioproject_accn'] == 'PRJNA258012'
        assert bioTree.bioproject['bioproject_id'] == '258012'
        assert bioTree.bioproject['name'].strip().split(' ')[0] == 'mRNA'
        assert bioTree.bioproject['name'].strip().split(' ')[-1] == 'environments'
        assert bioTree.bioproject['title'].strip().split(' ')[0] == 'mRNA'
        assert bioTree.bioproject['title'].strip().split(' ')[-1] == 'environments'
        assert bioTree.bioproject['description'].strip().split(' ')[0] == 'Our'
        assert bioTree.bioproject['description'].strip().split(' ')[-1] == 'level.'
        assert bioTree.bioproject['publication'] == '26732976'
        assert bioTree.bioproject['submission_id'] == 'SUB1475494'
        assert bioTree.bioproject['submission_date'] == '2014-08-11'
        assert bioTree.bioproject['last_update'] == '2016-04-20'
