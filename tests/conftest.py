#!/usr/bin/env python
""" Fixtures for running pytest. """

import pytest

from sramongo.mongo import start_mongo, stop_mongo
from sramongo.xml_helpers import xml_to_root


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
    stop_mongo(dbDir=mongo_folders[0])


@pytest.fixture(scope='session')
def sra_xml_root():
    fname = 'data/SRX971855.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def sra_xml_root2():
    fname = 'data/SRR3001915.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def sra_xml_root_PE():
    fname = 'data/ERR1662611.xml'
    return xml_to_root(fname)


@pytest.fixture(scope='session')
def pubmed_xml():
    fname = 'data/pubmed_26732976.xml'
    return xml_to_root(fname)
