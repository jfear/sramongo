Quick Start
===========

Installation
------------

MongoDB
+++++++

sramongo requires a working version of MongoDB community server >=3.4.1. Please
download the appropriate version
`here <https://www.mongodb.com/download-center#community>`__. There are also a
number of online hosts that will run a mongoDB server for free (small databases)
or relatively cheaply.

sramongo
++++++++

The preferred way to install sramongo is using the conda package manager, which
can be downloaded `here <https://conda.io/miniconda.html>`__. To install simply
run:

.. code:: bash

    conda install -c jfear sramongo

sra2mongo Usage
---------------

``sra2mongo`` is the command line tool provided by sramongo. To get a full set of
options run ``sra2mongo -h``. A simple query would look like:

.. code:: bash

   sra2mongo \
       --email john.smith@example.com \
       --dbDir $HOME/db \
       --logDir $HOME/logs \
       --db sra \
       --query '"Drosophila melanogaster[orgn]"'

The ``\`` allows for breaking the command on multiple lines. This command will
query the SRA for ``"Drosophila melanogaster"[orgn]``, download the XML for all
of the runs, and parse the XML into a database named 'sra' that is stored in the
user's home folder. If the mongo database was running then the script will
connect to the instance, if it is not running the script will startup a mongodb
instance.

A (see :ref:`dbFields`) for a list of database fields.

.. note::
    The query string is passed directly to SRA, so any query options such as
    [orgn], [pid], or [author] will work. Also queries can include boolean
    operators (i.e., AND, OR).

Querying the Database
---------------------

.. todo::
    Add section about querying the database using mongoengine. Until then follow
    `mongoengines docs <http://docs.mongoengine.org/guide/querying.html>`__

Extending Database
------------------

sramongo uses what is called an object document mapping or ODM (called
mongoengine_). An ODM allows mapping of attributes from a python class to the
database documents. While mongo is schemaless, the use of an ODM allows the
addition of validation and control of what fields are included in the database.
In order to add additional fields to the database you can either directly
interact with the database using the mongo API (i.e., pymongo) or you can
subclass the ODMs and add fields to them.

For example, let's imagine that you have a workflow that runs fastq screen on
each run and you want to store a list of potential contaminants in the database.
We would start by subclassing the ``Ncbi`` ODM and then adding a field to store
these
contaminants. My base ODM class uses mongoengine_ so you must use the
`mongoengine syntax <http://docs.mongoengine.org/guide/defining-documents.html>`__.

.. code:: python

   from mongoengine import ListField, StringField
   from sramongo import mongo_schema

   class myNcbi(mongo_schema.Ncbi):
       """My custom Run ODM."""
       fastq_screen = ListField(StrinField(), default=list)

.. note::
    When using mongoengine it stores a hidden variable in each document
    describing the mongoengine class used to create the document. For example,
    if you used sra2mongo to build the database the ``Ncbi`` document would have
    ``_cls = Ncbi``. **If you subclass Ncbi you must change this value to add
    your subclass name.** Given our example above:

    .. code:: python

        from mongoengine import connect

        client = connect('sra')
        client.sra.ncbi.update_many({}, {'$set': {'_cls': 'Ncbi.myNcbi'}})

.. _mongoengine: http://mongoengine.org

.. _pymongo: https://api.mongodb.com/python/current/
