.. _dbFields:

=================
sramongo mappings
=================

The database created by ``sra2mongo`` consists of 7 different documents:

.. contents:: :local:


Each document can be thought of as a separate notebook containing various bits
of information. Here I describe what information is stored by default.


Study
-----

.. autoclass:: sramongo.mongo_schema.Study


Sample
------

.. autoclass:: sramongo.mongo_schema.Sample

Experiment
----------

.. autoclass:: sramongo.mongo_schema.Experiment

Run
---

.. autoclass:: sramongo.mongo_schema.Run

BioSample
---------

.. autoclass:: sramongo.mongo_schema.BioSample

BioProject
----------

.. autoclass:: sramongo.mongo_schema.BioProject

Pubmed
------

.. autoclass:: sramongo.mongo_schema.Pubmed

