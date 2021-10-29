********
Features
********

    *Warning:*
    The sections below are outdated, will be updated soon

In the beginning the primary aim of **pypath** was to build networks from
multiple sources using an igraph object as the fundament of the integrated
data structure. From version 0.7 and 0.8 this design principle started to
change. Today **pypath** builds a number of different databases each having
**pandas.DataFrame** as a final format. Each of these integrates a specific
kind of data from various databases (e.g. protein complexes, interactions,
enzyme-PTM relationships, etc). **pypath** has many submodules with standalone
functionality which can be used in other modules and scripts. For example
the ID conversion module **pypath.mapping**.

Submodules perform various features, e.g. graph visualization, working with
rug compound data, searching drug targets and compounds in **ChEMBL**.

ID conversion
-------------

The ID conversion module ``utils.mapping`` can be used independently. It has
the feature to translate secondary UniProt IDs to primaries, and Trembl IDs to
SwissProt, using primary Gene Symbols to find the connections. This module
automatically loads and stores the necessary conversion tables. Many tables
are predefined, such as all the IDs in **UniProt mapping service,** while
users are able to load any table from **file** or **MySQL,** using the classes
provided in the module ``input_formats``.

Pathways
--------

**pypath** includes data and predefined format descriptions for more than 25
high quality, literature curated databases. The input formats are defined in
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

For three protein expression databases there are functions and modules for
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

The module ``pypath.curl`` provides a very flexible **download manager**
built on top of ``pycurl``. The classes ``pypath.curl.Curl()`` and
``pypath.curl.FileOpener`` accept numerous arguments to deal in a smart
way with local **cache**, authentication, redirects, uncompression, character
encodings, FTP and HTTP transactions, and many other stuff. Cache can grow to
several GBs, and takes place in ``~/.pypath/cache`` by default. If you
experience issues using ``pypath`` these are most often related to failed
downloads which often result nonsense cache contents. To debug such issues
you can see the cache file names and cache usage in the log, and you can use
the context managers in ``pypath.curl`` to show, delete or bypass the cache
for some particular method calls (``pypath.curl.cache_print_on()``,
``pypath.curl.cache_delete_on()`` and ``pypath.curl.cache_off()``.
You can always set up an alternative cache directory for the entire session
using the ``pypath.settings`` module.

The ``pypath.session`` and ``pypath.log`` modules take care of setting up
session level parameters and logging. Each session has a random 5 character
identifier e.g. ``y5jzx``. The default log file in this case is
``pypath_log/pypath-y5jzx.log``. The log messages are flushed every 2 seconds
by default. You can always change these things using the ``settings`` module.
In this module you can get and set the values of various parameters using
the ``pypath.settings.setup()`` and the ``pypath.settings.get()`` methods.

A simple **webservice** comes with this module: the ``server`` module based on
``twisted.web.server`` opens a custom port and serves plain text tables over
HTTP with REST style querying.
