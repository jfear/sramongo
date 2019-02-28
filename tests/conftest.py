#!/usr/bin/env python
""" Fixtures for running pytest. """
from pathlib import Path

import pytest

from sramongo.xml_helpers import xml_to_root


@pytest.fixture(scope='session')
def sra_xml_root():
    fname = Path(__file__).parent / 'data/SRX971855.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def sra_xml_root2():
    fname = Path(__file__).parent / 'data/SRR3001915.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def sra_xml_root_PE():
    fname = Path(__file__).parent / 'data/ERR1662611.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def pubmed_xml():
    fname = Path(__file__).parent / 'data/pubmed_26732976.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def bioproject_xml():
    fname = Path(__file__).parent / 'data/bioproject_PRJNA258012.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def biosample_xml():
    fname = Path(__file__).parent / 'data/biosample_SAMN02981965.xml'
    return xml_to_root(fname)
