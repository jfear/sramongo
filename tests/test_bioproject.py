#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.bioproject import BioProject


class TestPRJNA258012:
    @pytest.fixture(scope='class')
    def bio_etree(self):
        """ Element tree of single biosample. """
        fname = 'data/bioproject_PRJNA258012.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('DocumentSummary')

    @pytest.fixture(autouse=True, scope='class')
    def bioTree(self, bio_etree):
        return BioProject(bio_etree)

    def test_parse(self, bioTree):
        assert bioTree.BioProject == 'PRJNA258012'
        assert bioTree.name.strip().split(' ')[0] == 'mRNA'
        assert bioTree.name.strip().split(' ')[-1] == 'environments'
        assert bioTree.title.strip().split(' ')[0] == 'mRNA'
        assert bioTree.title.strip().split(' ')[-1] == 'environments'
        assert bioTree.description.strip().split(' ')[0] == 'Our'
        assert bioTree.description.strip().split(' ')[-1] == 'level.'
        assert bioTree.publication == '26732976'
        assert bioTree.submission_id == 'SUB1475494'
        assert bioTree.submission_date == '2014-08-11'
        assert bioTree.last_update == '2016-04-20'
