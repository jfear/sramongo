#!/usr/bin/env python
from xml.etree import ElementTree

import pytest

from sramongo.biosample import BioSampleParse


class TestSAMN02981065:
    @pytest.fixture(scope='class')
    def bio_etree(self):
        """ Element tree of single biosample. """
        fname = 'tests/data/biosample_SAMN02981965.xml'
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        return root.find('BioSample')

    @pytest.fixture(autouse=True, scope='class')
    def bioTree(self, bio_etree):
        return BioSampleParse(bio_etree)

    def test_parse(self, bioTree):
        assert bioTree.biosample['biosample_accn'] == 'SAMN02981965'
        assert bioTree.biosample['biosample_primary'] == 'SAMN02981965'
        assert bioTree.biosample['biosample_id'] == '2981965'
        assert bioTree.biosample['sample_id'] == 'SRS679015'
        assert bioTree.biosample['GEO'] == 'GSM1471477'
        assert bioTree.biosample['title'] == 'DGRP563 M_E3_2_L3'
        assert bioTree.biosample['tax_id'] == '7227'
        assert bioTree.biosample['tax_name'] == 'Drosophila melanogaster'
        assert bioTree.biosample['organism_name'] == 'Drosophila melanogaster'
        assert bioTree.biosample['institute'] == 'Developmental Genomics, LCDB, NIDDK, NIH'
        assert bioTree.biosample['access'] == 'public'
        assert bioTree.biosample['publication_date'] == '2015-12-22'
        assert bioTree.biosample['last_update'] == '2015-12-22'
        assert bioTree.biosample['submission_date'] == '2014-08-11'
        assert bioTree.biosample['contacts'][0]['email'] == 'briano@helix.nih.gov'
        assert bioTree.biosample['contacts'][0]['first_name'] == 'Brian'
        assert bioTree.biosample['contacts'][0]['last_name'] == 'Oliver'
        assert bioTree.biosample['models'][0] == 'Generic'
        attr = {
                'source_name': 'Whole body',
                'strain': 'DGRP-563',
                'dev_stage': 'Adult',
                'sex': 'male',
                'tissue': 'Whole body',
                }
        for attribute in bioTree.biosample['attributes']:
            assert attribute['value'] == attr[attribute['name']]
