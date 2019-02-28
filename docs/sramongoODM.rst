.. _dbFields:

=================
sramongo mappings
=================

The database created by ``sra2mongo`` consists of a single document that is
organized hierarchically:

.. contents:: :local:



This can be thought of as a giant JSON or python dict which various levels can
be accessed by indexing through (e.g., ``ncbi.sra.run.run_id``). MongoDB has a
very nice querying system which allows easy searching through the document.

.. note::
    One downside of storing all of this information as a single document is that
    mongoDB has a max document size of 16 MB. This is more than enough for
    storing metadata and text, but if you start adding data tables you may hit
    this limit.

ncbi
====

This is the top level document. Information from each database is stored under
its name. As I add data normalization steps I intend to aggregate data from the
different databases and store them up in this top level document.


sra
---

This stores all from the Sra. There are also a couple of summary fields that are
stored at this level. Each section of the SRA record are represented as
subdocuments.

.. autoclass:: sramongo.models.SraDocument

organization
++++++++++++

.. autoclass:: sramongo.models.Organization

submission
++++++++++

.. autoclass:: sramongo.models.Submission

study
+++++

.. autoclass:: sramongo.models.Study


run
+++

.. autoclass:: sramongo.models.Run

sample
++++++

.. autoclass:: sramongo.models.Sample


biosample
---------

Information from the BioSample database is stored here.

.. autoclass:: sramongo.models.BioSample


bioproject
----------

Information from the BioProject database is stored here.

.. autoclass:: sramongo.models.BioProject


pubmed
------

Information from the Pubmed is stored here.

.. autoclass:: sramongo.models.Pubmed

