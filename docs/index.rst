sra2mongo User Documentation
====================================

.. _`SRA <https://www.ncbi.nlm.nih.gov/sra>`

sra2mongo is a python library and command line tool that queries the `SRA`_
and dumps all relevant information into a mongo database. Mongo is a
document based database, meaning that it stores information in`key:value` pairs.
One major advantage to this over a SQL solution is that mongo does not need a
defined schema; meaning attributes can be arbitrarily added or removed and don't
need to be the same across records. What this means is that you can run
sra2mongo to query and populate your database, then incorporate and modify the
database as part of processing pipelines.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   flags



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
