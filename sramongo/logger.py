"""
Configures a global logger
"""
from sys import stderr, stdout
import logging


class InfoFilter(logging.Filter):
    """ Simple filter to only output debug and info to stdout.

    http://stackoverflow.com/questions/16061641/python-logging-split-between-stdout-and-stderr

    """
    def filter(self, rec):
        return rec.levelno in (logging.DEBUG, logging.INFO)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Error handler to stderr
h1 = logging.StreamHandler(stream=stderr)
h1.setLevel(logging.WARNING)
h1.setFormatter(formatter)

# INFO and DEBUG handler to stdout
h2 = logging.StreamHandler(stream=stdout)
h2.setLevel(logging.DEBUG)
h2.setFormatter(formatter)
h2.addFilter(InfoFilter())

# Add Handles to logger
logger.addHandler(h1)
logger.addHandler(h2)
