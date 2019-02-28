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
sramongo can be installed using pip:

.. code:: bash

    pip install git+https://github.com/jfear/sramongo


sra2mongo Usage
---------------

``sra2mongo`` is the command line tool provided by sramongo. To get a full set of
options run ``sra2mongo -h``. A simple query would look like:

.. code:: bash

   sra2mongo \
       --email john.smith@example.com \
       --query '"Drosophila melanogaster"[orgn]'

The ``\`` allows for breaking the command on multiple lines. This command will
query the SRA for ``"Drosophila melanogaster"[orgn]``, download the XML for all
of the runs, and parse the XML into a database named 'sramongo'.

A (see :ref:`dbFields`) for a list of database fields.

.. note::
    The query string is passed directly to SRA, so any query options such as
    [orgn], [pid], or [author] will work. Also queries can include boolean
    operators (i.e., AND, OR).

Querying the Database
---------------------

.. todo::
    Add section about querying the database using mongoengine and pymongo. Until
    then follow
    `mongoengines docs <http://docs.mongoengine.org/guide/querying.html>`__

.. _mongoengine: http://mongoengine.org

.. _pymongo: https://api.mongodb.com/python/current/
