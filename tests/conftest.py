#!/usr/bin/env python
""" Fixtures for running pytest. """
from xml.etree import ElementTree

import pytest

from sramongo.mongo import start_mongo
from sramongo.sra import SraExperiment
from sramongo.biosample import BioSampleParse


@pytest.fixture(scope='session')
def mongo_folders(tmpdir_factory):
    """ Create TMPDIR for database and logs. """
    data = tmpdir_factory.mktemp('db')
    log = tmpdir_factory.mktemp('log')
    return str(data), str(log)


@pytest.fixture(scope='session')
def mongoDB(mongo_folders):
    """ Start the database server.

    Server will be killed after the fixture is no longer neaded.
    """
    mongoDB = start_mongo(dbDir=mongo_folders[0], logDir=mongo_folders[1])
    yield mongoDB
    print('Shutting down mongoDB.')
    mongoDB.kill()


@pytest.fixture(scope='session')
def sraExperiment():
    fname = 'tests/data/sra_SRR3001915.xml'
    tree = ElementTree.parse(fname)
    root = tree.getroot()
    return SraExperiment(root.find('EXPERIMENT_PACKAGE'))


@pytest.fixture(scope='session')
def bioSample():
    fname = 'tests/data/biosample_SAMN02981965.xml'
    tree = ElementTree.parse(fname)
    root = tree.getroot()
    return BioSampleParse(root.find('BioSample'))
