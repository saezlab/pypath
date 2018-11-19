.. pypath documentation master file, created by
   sphinx-quickstart on Fri Oct 19 15:13:54 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

############################################################
PyPath: *A Python module for molecular interaction networks*
############################################################


**pypath** is a Python package built around igraph to work with molecular
network representations e.g. protein, miRNA and drug compound interaction
networks.

:note: ``pypath`` supports both Python 2.7 and Python 3.6+. In the beginning,
    pypath has been developed only for Python 2.7. Then the code have been
    adjusted to Py3 however we can not guarantee no incompatibilities
    remained. If you find any method does not work please submit an issue on
    github. For few years I develop and test ``pypath`` in Python 3. Therefore
    this is the better supported Python variant.

:documentation: http://saezlab.github.io/pypath
:issues: https://github.com/saezlab/pypath/issues

.. toctree::
    :maxdepth: 1
    :caption: Contents:

    installation
    main
    webservice
    changelog


Features
========

The primary aim of **pypath** is to build up networks from multiple sources on
one igraph object. **pypath** handles ambiguous ID conversion, reads custom
edge and node attributes from text files and **MySQL**.

Submodules perform various features, e.g. graph visualization, working with
rug compound data, searching drug targets and compounds in **ChEMBL**.

ID conversion
-------------

The ID conversion module ``mapping`` can be used independently. It has the
feature to translate secondary UniProt IDs to primaries, and Trembl IDs to
SwissProt, using primary Gene Symbols to find the connections. This module
automatically loads and stores the necessary conversion tables. Many tables
are predefined, such as all the IDs in **UniProt mapping service,** while
users are able to load any table from **file** or **MySQL,** using the classes
provided in the module ``input_formats``.

Pathways
--------

**pypath** includes data and predefined format descriptions for more than 25
high quality, literature curated databases. The inut formats are defined in
the ``data_formats`` module. For some resources data downloaded on the fly,
where it is not possible, data is redistributed with the module. Descriptions
and comprehensive information about the resources is available in the
``descriptions`` module.

Structural features
-------------------

One of the modules called ``intera`` provides many classes for representing
structures and mechanisms behind protein interactions. These are ``Residue``
(optionally mutated), ``Motif``, ``Ptm``, ``Domain``, ``DomainMotif``,
``DomainDomain`` and ``Interface``. All these classes have ``__eq__()``
methods to test equality between instances, and also ``__contains__()``
methods to look up easily if a residue is within a short motif or protein
domain, or is the target residue of a PTM.

Sequences
---------

The module ``seq`` contains a simple class for quick lookup any residue or
segment in **UniProt** protein sequences while being aware of isoforms.

Tissue expression
-----------------

For 3 protein expression databases there are functions and modules for
downloading and combining the expression data with the network. These are the
Human Protein Atlas, the ProteomicsDB and GIANT. The ``giant`` and
``proteomicsdb`` modules can be used also as stand alone Python clients for
these resources.

Functional annotations
----------------------

**GSEA** and **Gene Ontology** are two approaches for annotating genes and
gene products, and enrichment analysis technics aims to use these annotations
to highlight the biological functions a given set of genes is related to. Here
the ``enrich`` module gives abstract classes to calculate enrichment
statistics, while the ``go`` and the ``gsea`` modules give access to GO and
GSEA data, and make it easy to count enrichment statistics for sets of genes.

Drug compounds
--------------

**UniChem** submodule provides an interface to effectively query the UniChem
service, use connectivity search with custom settings, and translate SMILEs to
ChEMBL IDs with ChEMBL web service.

**ChEMBL** submodule queries directly your own ChEMBL MySQL instance, has the
features to search targets and compounds from custom assay types and
relationship types, to get activity values, binding domains, and action types.
You need to download the ChEMBL MySQL dump, and load into your own server.

Technical
---------

**MySQL** submodule helps to manage MySQL connections and track queries. It is
able to run queries parallely to optimize CPU and memory usage on the server,
handling queues, and serve the result by server side or client side storage.
The ``chembl`` and potentially the ``mapping`` modules rely on this ``mysql``
module.

The most important function in module ``dataio`` is a very flexible **download
manager** built around ``curl``. The function ``dataio.curl()`` accepts
numerous arguments, tries to deal in a smart way with local **cache,**
authentication, redirects, uncompression, character encodings, FTP and HTTP
transactions, and many other stuff. Cache can grow to several GBs, and takes
place in ``./cache`` by default. Please be aware of this, and use for example
symlinks in case of using multiple working directories.

A simple **webservice** comes with this module: the ``server`` module based on
``twisted.web.server`` opens a custom port and serves plain text tables over
HTTP with REST style querying.


OmniPath in R
=============

You can download the data from the webservice and load into R. Look
`here <https://github.com/saezlab/pypath/tree/master/r_import>`_ for an
example.


.. only:: html

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
