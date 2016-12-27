#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.biosample import BioSample


class TestSAMN02981065:
    @pytest.fixture(scope='class')
    def bio_etree(self):
        """ Element tree of single biosample. """
        fname = 'data/biosample_SAMN02981965.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('BioSample')

    @pytest.fixture(autouse=True, scope='class')
    def bioTree(self, bio_etree):
        return BioSample(bio_etree)

    def test_parse(self, bioTree):
        assert bioTree.BioSample == 'SAMN02981965'
        assert bioTree.SRA == 'SRS679015'
        assert bioTree.GEO == 'GSM1471477'
        assert bioTree.title == 'DGRP563 M_E3_2_L3'
        assert bioTree.tax_id == '7227'
        assert bioTree.tax_name == 'Drosophila melanogaster'
        assert bioTree.organism_name == 'Drosophila melanogaster'
        assert bioTree.institute == 'Developmental Genomics, LCDB, NIDDK, NIH'
        assert bioTree.access == 'public'
        assert bioTree.publication_date == '2015-12-22'
        assert bioTree.last_update == '2015-12-22'
        assert bioTree.submission_date == '2014-08-11'
        assert bioTree.contacts[0]['email'] == 'briano@helix.nih.gov'
        assert bioTree.contacts[0]['first_name'] == 'Brian'
        assert bioTree.contacts[0]['last_name'] == 'Oliver'
        assert bioTree.models[0] == 'Generic'
        assert bioTree.attributes['source_name'] == 'Whole body'
        assert bioTree.attributes['strain'] == 'DGRP-563'
        assert bioTree.attributes['dev_stage'] == 'Adult'
        assert bioTree.attributes['sex'] == 'male'
        assert bioTree.attributes['tissue'] == 'Whole body'
