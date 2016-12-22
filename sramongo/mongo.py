#!/usr/bin/python
""" Helper functions for using mongoDB. """

import os
from time import sleep
from subprocess import Popen

class MongoDBFolderError(Exception):
    pass


def start_mongo(dbDir=None, logDir=None, port=27017, command_line_args=''):
    """ Start the mongod server.

    Use a linux subprocess to start the mongoDB server.

    Parameters
    ----------
    dbDir: str
        The path to the mongoDB database direcotry. If not given will try to
        use `MONGODB_DATA_DIR` environmental variable.
    logDir: str
        The path to the mongoDB log directory. If not given will try to use
        `MONGODB_LOG_DIR` environmental variable.
    port: int
        The port to run the server on.
    command_line_args: str
        Additional mongod command line arguments. Passed directly to the mongod
        command.

    """
    if dbDir is None:
        try:
            dbDir = os.environ('MONGODB_DATA_DIR')
        except KeyError:
            raise MongoDBFolderError('No database direcotry given. '
                                     'Either set the MONGODB_DATA_DIR environmental variable '
                                     'or provide a dbDir.')
    if not os.path.exists(dbDir):
        raise MongoDBFolderError('Database direcotry must exists.')

    if logDir is None:
        try:
            logDir = os.environ('MONGODB_LOG_DIR')
        except KeyError:
            raise MongoDBFolderError('No log direcotry given. '
                                     'Either set the MONGODB_LOG_DIR environmental variable '
                                     'or provide a logDir.')
    if not os.path.exists(logDir):
        raise MongoDBFolderError('Log direcotry must exists.')

    log = os.path.join(logDir, 'mongoDB.log')

    if command_line_args:
        cmd = ['mongod', '--dbpath', dbDir, '--logpath', log, '--port', str(port), command_line_args]
    else:
        cmd = ['mongod', '--dbpath', dbDir, '--logpath', log, '--port', str(port)]

    process =  Popen(cmd)
    sleep(10)
    return process
