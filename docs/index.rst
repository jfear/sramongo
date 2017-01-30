sramongo User Documentation
====================================

sramongo is a python library and command line tool that queries NCBI's sequence
read archive(`SRA <https://www.ncbi.nlm.nih.gov/sra>`_) and dumps all relevant
information into a mongo database. `Mongo <https://www.mongodb.com/>`_ is a
popular document based database, which stores information as `key:value` pairs;
unlike a typical relational database (e.g., SQL) which stores data as rows in
multiple related tables. One major advantage of a document based database to a
relation database is that there is no need a defined schema; meaning attributes
can be arbitrarily added or removed and don't need to be the same across
records. What this means is that you can run sra2mongo to query and populate
your database, then freely modify or add new fields/documents as part of a
processing pipeline.

Roughly speak sramongo is made up of 3 parts. The first is a parser for SRA XML,
the second is an object relational mapper to allow easy interface with mongo,
and the third is a command line utility which uses Biopython and Entrez
utilities to query the SRA and download the resulting XML.

.. warning::
    Please use this tool responsibly, querying the SRA and dumping large amounts
    of data can be taxing on their system and may result in blacklisting of your
    IP address.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickStart
   sramongoODM
   flags
   sraConstants



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
