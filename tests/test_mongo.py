#!/usr/bin/env python
""" Fixtures for running pytest """

from time import sleep
from pymongo import MongoClient
from sramongo.mongo import start_mongo, stop_mongo, MongoDB


def test_start_mongo(mongo_folders):
    """ Test starting the database server. """
    mongoDB = start_mongo(dbDir=mongo_folders[0], logDir=mongo_folders[1])
    assert isinstance(mongoDB.pid, int)
    stop_mongo(dbDir=mongo_folders[0])
    sleep(3)
    assert mongoDB.poll() is not None


def test_start_mongo_command_line_args(mongo_folders):
    """ Test starting the database server. """
    mongoDB = start_mongo(dbDir=mongo_folders[0], logDir=mongo_folders[1],
                          command_line_args='--verbose')
    assert isinstance(mongoDB.pid, int)
    assert mongoDB.poll() is None
    stop_mongo(dbDir=mongo_folders[0])


def test_start_context_manager(mongo_folders):
    with MongoDB(dbDir=mongo_folders[0], logDir=mongo_folders[1]) as mongoDB:
        assert isinstance(mongoDB.pid, int)
        assert mongoDB.poll() is None
    sleep(3)
    assert mongoDB.poll() is not None


def test_mongo_connect(mongoDB):
    """ Test basic functionality of mongo and pymongo. """
    client = MongoClient()
    assert sorted(client.database_names()) == sorted(['admin', 'local'])

    db = client.test_database
    collection = db.test_collection
    collection.insert_one({'name': 'test'})

    assert 'test_database' in client.database_names()
    assert db.collection_names() == ['test_collection']
    assert collection.find_one({}, {'name': 1, '_id': 0}) == {'name': 'test'}

    client.close()
