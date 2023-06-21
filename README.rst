============================================================================
*pypath:* A Python module for molecular signaling prior knowledge processing
============================================================================

|Demo|

OmniPath
========

Are you interested in OmniPath data? Check out our R package OmnipathR_,
the most popular and most versatile access point to OmniPath, a database
built from more than 150 original resources. If you use Python and don't
need to build the database yourself, try our `Python client`_. Read more
about the `web service here`_.

.. _OmnipathR: https://r.omnipathdb.org
.. _`Python client`: https://github.com/saezlab/omnipath
.. _`web service here`: https://pypath.omnipathdb.org/webservice.html

Do you need pypath?
===================

Pypath is the database builder of OmniPath. For most people the data
distributed in OmniPath is satisfying (see above), they don't really need
pypath. Typically you need pypath to:

* Build a custom or very fresh version of the OmniPath database(s)
* Use one of the utilities such as ID translation, homology translation, etc.
  (see the `utils module`_)
* Access the raw or preprocessed data directly from the original resources
  (see the `inputs module`_)

.. _`utils module`: https://github.com/saezlab/pypath/tree/master/pypath/utils
.. _`inputs module`: https://github.com/saezlab/pypath/tree/master/pypath/inputs

Installation
============

**From PyPI:**

.. code:: bash

    pip install pypath-omnipath

**From Git:**

.. code:: bash

    pip install git+https://github.com/saezlab/pypath.git

Docs
====

Read the `reference documentation`_ or check out the tutorials_. The most
comprehensive guide to *pypath* is `The Pypath Book`_.

.. _`reference documentation`: https://pypath.omnipathdb.org/
.. _tutorials: https://workflows.omnipathdb.org/
.. _`The Pypath Book`: https://pypath.omnipathdb.org/notebooks/manual.html

Get help
========

Should you have a question or experiencing an issue, please write us by
the `Github issues`_ page.

Features
========

**pypath** is a Python module for processing molecular biology data resources,
combining them into databases and providing a versatile interface in Python
as well as exporting the data for access through other platforms such as
R_, `web service`_, Cytoscape_ and BEL (Biological Expression Language).

.. _R: https://r.omnipathdb.org/
.. _`web service`: https://omnipathdb.org/
.. _Cytoscape: https://apps.cytoscape.org/apps/omnipath

**pypath** provides access to more than 100 resources! It builds 5 major
combined databases and within these we can distinguish different datasets.
The 5 major databases are interactions (molecular interaction network or
pathways), enzyme-substrate relationships, protein complexes, molecular
annotations (functional roles, localizations, and more) and inter-cellular
communication roles.

**pypath** consists of a number of submodules and each of them again contains
a number of submodules. Overall **pypath** consists of around 100 modules.
The most important higher level submodules:

* *pypath.core:* contains the database classes e.g. network, complex,
  annotations, etc
* *pypath.inputs:* contains the resource specific methods which directly
  downlad and preprocess data from the original sources
* *pypath.omnipath:* higher level applications, e.g. a database manager, a
  web server
* *pypath.utils:* stand alone useful utilities, e.g. identifier translator,
  Gene Ontology processor, BioPax processor, etc

Integrated databases
--------------------

In the beginning the primary aim of ``pypath`` was to build networks from
multiple sources using an igraph object as the fundament of the integrated
data structure. From version 0.7 and 0.8 this design principle started to
change. Today ``pypath`` builds a number of different databases, exposes them
by a rich API and each of them can be converted to ``pandas.DataFrame``.
The modules and classes responsible for the integrated databases are located
in ``pypath.core``. The five main databases are the followings:

* *network* - ``core.network``
* *enzyme-substrate* - ``core.enz_sub``
* *complexes* - ``core.complex``
* *annotations* - ``core.annot``
* *intercell* - ``core.intercell``

Some of the databases have different variants (e.g. PPI and transcriptional
network) and all can be customized by many parameters.

Database management
-------------------

The databases above can be loaded by calling the appropriate classes.
However building the databases require time and memory so we want to avoid
building them more often than necessary or keeping more than one copies
in the memory. Some of the modules listed above have a method ``get_db``
which ensures only one instance of the database is loaded. But there is a
more full featured database management system available in **pypath**,
this is the **pypath.omnipath** module. This module is able to build the
databases, automatically saves them to ``pickle`` files and loads them from
there in subsequent sessions. **pypath** comes with a number of database
definitions and users can add more. The ``pickle`` files are located by
default in the ``~/.pypath/pickles/`` directory. With the ``omnipath``
module it's easy to get an instance of a database. For example to get the
`omnipath` PPI network dataset:

.. code:: python

    from pypath import omnipath
    op = omnipath.db.get_db('omnipath')

**Important:** Building the databases for the first time requires the
download of several MB or GB of data from the original resources. This
normally takes long time and is prone of errors (e.g. truncated or empty
downloads due to interrupted HTTP connection). In this case you should check
the log to find the path of the problematic cache file, check the contents
of this file to find out the reason and possibly delete the file to ensure
another download attempt when you call the database build again. Sometimes
the original resources change their content or go offline. If you encounter
such case please open an issue at https://github.com/saezlab/pypath/issues
so we can fix it in ``pypath``. Once all the necessary contents are
downloaded and stored in the cache, the database builds are much faster,
but still can take minutes.

Further modules in pypath
-------------------------

Apart from the databases, **pypath** has many submodules with standalone
functionality which can be used in other modules and scripts. Below we
present a few of these.

ID conversion
-------------

The ID conversion module ``utils.mapping`` translates between a large variety
of gene, protein, miRNA and small molecule ID types. It has the feature to
translate secondary UniProt ACs to primaries, and Trembl ACs to SwissProt,
using primary Gene Symbols to find the connections. This module automatically
loads and stores the necessary conversion tables. Many tables
are predefined, such as all the IDs in **UniProt mapping service,** while
users are able to load any table from **file** using the classes provided
in the module ``input_formats``. An example how to translate identifiers:

.. code:: python

    from pypath.utils import mapping
    mapping.map_name('P00533', 'uniprot', 'genesymbol')
    # {'EGFR'}


Homology translation
--------------------

The ``pypath.utils.homology`` module is able to find the orthologs of genes
between two organisms. It uses data both from NCBI HomoloGene, Ensembl and
UniProt. This module is really simple to use:

.. code:: python

    from pypath.utils import homology
    homology.translate('P00533', 10090) # translating the human EGFR to mouse
    # ['Q01279'] # it returns the mouse Egfr UniProt AC

It is able to handle any ID type supported by ``pypath.utils.mapping``.
Alternatively, you can access a complete dictionary of orthologous genes,
or translate columns in a pandas data frame.

FAQ
===

**Does it run on my old Python?**

Most likely it doesn't. The oldest supported version, currently 3.9, is
defined in our `pyproject.toml`_.

.. _`pyproject.toml`: https://github.com/saezlab/pypath/blob/master/pyproject.toml

**Is there something similar in R?**

`OmniPath's R client`_, besides accessing data from OmniPath, provides many
similar services as pypath: `ID translation`_, `homology translation`_,
`taxonomy support`_, `GO support`_, and many more.

.. _`OmniPath's R client`: https://r.omnipathdb.org
.. _`ID translation`: https://r.omnipathdb.org/reference/translate_ids.html
.. _`homology translation`: https://r.omnipathdb.org/reference/homologene_uniprot_orthology.html
.. _`taxonomy support`: https://r.omnipathdb.org/reference/ncbi_taxid.html
.. _`GO support`: https://r.omnipathdb.org/reference/go_annot_download.html

`Questions about OmniPath`_

.. _`Questions about OmniPath`: https://omnipathdb.org/#faq

Contact
=======

We prefer to keep all communication within the `Github issues`_. About private
or sensitive matters feel free to contact us by omnipathdb@gmail.com.

.. _`Github issues`: https://github.com/saezlab/pypath/issues

Impressum
=========

The development of ``pypath`` is coordinated by `Dénes Türei`_ in the
`Saez Lab`_, with the contribution of developers and scientists from
other groups:

* Erva Ulusoy, Melih Darcan, Ömer Kaan Vural, Tennur Kılıç, Elif Çevrim,
  Bünyamin Şen, Atabey Ünlü and Mert Ergün in the
  `HU Biological Data Science Lab (PI: Tunca Doğan)`_ created many new input
  modules in `pypath`;
* Leila Gul, Dezső Módos, Márton Ölbei and Tamás Korcsmáros in the
  `Korcsmaros Lab`_ contributed to the overall design of OmniPath, the
  design and implementation of the intercellular communication database,
  and with various case studies and tutorials;
* Michael Klein from the group of `Fabian Theis`_ developed the
  `Python client`_ for the OmniPath web service;
* Charles Tapley Hoyt and Daniel Domingo-Fernández added the BEL export
  module.
* From the `Saez Lab`_, Olga Ivanova introduced the resource manager in
  `pypath`, Sophia Müller-Dott added the CollecTRI gene regulatory network,
  while Nicolàs Palacio, Sebastian Lobentanzer and Ahmet Rifaioglu
  have done various maintenance and refactoring works. Aurelien Dugourd and
  Christina Schmidt helped with the design of the metabolomics related
  datasets and services.
* The `R package`_ and the `Cytoscape app`_ are developed and maintained by
  Francesco Ceccarelli, Attila Gábor, Alberto Valdeolivas, Dénes Türei and
  Nicolàs Palacio;
* The first logo of OmniPath has been designed by Jakob Wirbel (Saez Lab),
  the current logo by Dénes Türei, while the cover graphics for Nature Methods
  is the work of Spencer Phillips from EMBL-EBI.

.. _`Saez Lab`: https://saezlab.org/
.. _`HU Biological Data Science Lab (PI: Tunca Doğan)`: https://yunus.hacettepe.edu.tr/~tuncadogan/
.. _`Dénes Türei`: https://denes.omnipathdb.org/
.. _`R package`: https://r.omnipathdb.org
.. _`Cytoscape app`: https://apps.cytoscape.org/apps/omnipath
.. _`Fabian Theis`: https://www.helmholtz-munich.de/en/icb/research-groups/theis-lab/
.. _`Korcsmaros Lab`: https://korcsmaroslab.org/

History and releases
====================

See here_ a bird eye view of pypath's development history. For more details
about recent developments see the `Github releases`_.

.. _here: https://pypath.omnipathdb.org/releasehistory.html
.. _`Github releases`: https://github.com/saezlab/pypath/releases

.. |Demo| image:: https://raw.githubusercontent.com/saezlab/pypath/master/docs/source/_static/img/pypath-demo.webp
