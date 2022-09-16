============================================================================
*pypath:* A Python module for molecular signaling prior knowledge processing
============================================================================

OmniPath
========

Are you interested in OmniPath data? Check out our R package OmnipathR_,
the most popular and most versatile access point to OmniPath, a database
built from more than 150 original resources. If you use Python and don't
need to build the database yourself, try our `Python client`_.

.. _OmnipathR: https://github.com/saezlab/OmnipathR
.. _`Python client`: https://github.com/saezlab/omnipath

    **Important:** New module structure and new network API (January 2020)

    Around the end of December we added a new network API to ``pypath`` which
    is not based on ``igraph`` any more and provides a modular and versatile
    access interface to the network data (since version ``0.9``). In January
    we reorganized the submodules in ``pypath`` in order to create a clear
    structure (since version ``0.10``). These are important milestones
    towards version ``1.0`` and we hope they will make ``pypath`` more
    convenient to use for everyone. By 18 February we merged these changes
    to the master branch however the *pypath guide* is still to be updated.
    Apologies for this inconvenience and please don't hesitate to ask
    questions by opening an issue on github. The old ``igraph`` based network
    class is still available in the ``pypath.legacy`` module.

:Py2/3: The oldest suitable Python is version 3.9, as defined in
        ``pyproject.toml``. Until about April 2022 we still kept most of
        the module compatible with Python 2 and earlier Python 3 versions.
        Since then we support only the recent Pythons.

:documentation: https://saezlab.github.io/pypath
:issues: https://github.com/saezlab/pypath/issues
:contact: omnipathdb@gmail.com
:developers: ``pypath`` is developed in the Saez Lab (https://saezlab.org) by
  Dénes Türei, Sebastian Lobentanzer and Ahmet Rifaioglu, and Erva Ulusoy
  and Tennur Kılıç in Volkan Atalay's group
  (https://blog.metu.edu.tr/vatalay/). Olga Ivanova and Nicolàs Palacio also
  contributed in the past. The R package and the Cytoscape app are developed
  and maintained by Francesco Ceccarelli, Attila Gábor, Alberto Valdeolivas,
  Dénes Türei and Nicolàs Palacio. The `Python client`_ for the OmniPath web
  service has been developed and is maintained by Michael Klein in the group
  of Fabian Theis.

.. _`Python client`: https://github.com/saezlab/omnipath

**pypath** is a Python module for processing molecular biology data resources,
combining them into databases and providing a versatile interface in Python
as well as exporting the data for access through other platforms such as
the R (the OmnipathR R/Bioconductor package;
https://github.com/saezlab/OmnipathR), web service (at
https://omnipathdb.org), Cytoscape (the OmniPath Cytoscape app;
https://apps.cytoscape.org/apps/omnipath) and BEL
(Biological Expression Language).

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


Webservice
==========

**New webservice** from 14 June 2018: the queries slightly changed, have been
largely extended. See the examples below.

The webservice implements a very simple REST style API, you can make requests
by the HTTP protocol (browser, wget, curl or whatever). After defining the
query type and optionally a set of molecular entities (proteins) you can
add further GET parameters encoded in the URL.

Query types
-----------

The webservice currently recognizes 7 types of queries: ``interactions``,
``enz_sub``, ``annotations``, ``complexes``, ``intercell``, ``queries`` and
``info``.
The query types ``resources``, ``network`` and ``about`` have not been
implemented yet in the new webservice.

Interaction datasets
--------------------

The instance of the ``pypath`` webserver running at the domain
https://omnipathdb.org/, serves not only the OmniPath data but also other
datasets. Each of them has a short name what you can use in the queries
(e.g. ``&datasets=omnipath,pathwayextra``).

* ``omnipath``: the OmniPath data as defined in the paper, an arbitrary
  optimum between coverage and quality
* ``pathwayextra``: activity flow interactions without literature reference
* ``kinaseextra``: enzyme-substrate interactions without literature reference
* ``ligrecextra``: ligand-receptor interactions without literature reference
* ``dorothea``: transcription factor (TF)-target interactions from DoRothEA
* ``tf_target``: transcription factor (TF)-target interactions from other
  sources
* ``mirnatarget``: miRNA-mRNA and TF-miRNA interactions

TF-target interactions from DoRothEA, a large collection additional
enzyme-substrate interactions, and literature curated miRNA-mRNA interacions
combined from 4 databases.

Mouse and rat
-------------

Except the miRNA interactions all interactions are available for human, mouse
and rat. The rodent data has been translated from human using the NCBI
Homologene database. Many human proteins do not have known homolog in rodents
hence rodent datasets are smaller than their human counterparts. Note, if you
work with mouse omics data you might do better to translate your dataset to
human (for example using the ``pypath.homology`` module) and use human
interaction data.


Examples
--------

A request without any parameter provides the main webpage:

    https://omnipathdb.org

The ``info`` returns a HTML page with comprehensive information about the
resources. The list here should be and will be updated as currently OmniPath
includes much more databases:

    https://omnipathdb.org/info

Molecular interaction network
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``interactions`` query accepts some parameters and returns interactions in
tabular format. This example returns all interactions of EGFR (P00533), with
sources and references listed.

    https://omnipathdb.org/interactions/?partners=P00533&fields=sources,references

By default only the OmniPath dataset used, to include any other dataset you
have to set additional parameters. For example to query the transcriptional regulators of EGFR:

    https://omnipathdb.org/interactions/?targets=EGFR&types=transcriptional

The DoRothEA database assigns confidence levels to the interactions. You
might want to select only the highest confidence, *A* category:

    https://omnipathdb.org/interactions/?targets=EGFR&types=transcriptional&dorothea_levels=A

Show the transcriptional targets of Smad2 homology translated to rat including
the confidence levels from TF Regulons:

    https://omnipathdb.org/interactions/?genesymbols=1&fields=type,ncbi_tax_id,dorothea_level&organisms=10116&sources=Smad2&types=transcriptional

Query interactions from PhosphoNetworks which is part of the *kinaseextra*
dataset:

    https://omnipathdb.org/interactions/?genesymbols=1&fields=sources&databases=PhosphoNetworks&datasets=kinaseextra

Get the interactions from Signor, SPIKE and SignaLink3:

    https://omnipathdb.org/interactions/?genesymbols=1&fields=sources,references&databases=Signor,SPIKE,SignaLink3

All interactions of MAP1LC3B:

    https://omnipathdb.org/interactions/?genesymbols=1&partners=MAP1LC3B

By default ``partners`` queries the interaction where either the source or the
arget is among the partners. If you set the ``source_target`` parameter to
``AND`` both the source and the target must be in the queried set:

    https://omnipathdb.org/interactions/?genesymbols=1&fields=sources,references&sources=ATG3,ATG7,ATG4B,SQSTM1&targets=MAP1LC3B,MAP1LC3A,MAP1LC3C,Q9H0R8,GABARAP,GABARAPL2&source_target=AND

As you see above you can use UniProt IDs and Gene Symbols in the queries and
also mix them. Get the miRNA regulating NOTCH1:

    https://omnipathdb.org/interactions/?genesymbols=1&fields=sources,references&datasets=mirnatarget&targets=NOTCH1

Note: with the exception of mandatory fields and genesymbols, the columns
appear exactly in the order you provided in your query.

Enzyme-substrate interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another query type available is ``ptms`` which provides enzyme-substrate
interactions. It is very similar to the ``interactions``:

    https://omnipathdb.org/enz_sub?genesymbols=1&fields=sources,references,isoforms&enzymes=FYN

Is there any ubiquitination reaction?

    https://omnipathdb.org/ens_sub?genesymbols=1&fields=sources,references&types=ubiquitination

And acetylation in mouse?

    https://omnipathdb.org/ptms?genesymbols=1&fields=sources,references&types=acetylation&organisms=10090

Rat interactions, both directly from rat and homology translated from human,
from the PhosphoSite database:

    https://omnipathdb.org/enz_sub?genesymbols=1&fields=sources,references&organisms=10116&databases=PhosphoSite,PhosphoSite_noref


Molecular complexes
^^^^^^^^^^^^^^^^^^^

The ``complexes`` query provides a comprehensive database of more than 22,000
protein complexes. For example, to query all complexes from CORUM and PDB
containing MTOR (P42345):

    https://omnipathdb.org/complexes?proteins=P42345&databases=CORUM,PDB


Annotations
^^^^^^^^^^^

The ``annotations`` query provides a large variety of data about proteins,
complexes and in the future other kinds of molecules. For example an
annotation can tell if a protein is a kinase, or if it is expressed in the
hearth muscle. These data come from dozens of databases and each kind of
annotation record contains different fields. Because of this here we have
a ``record_id`` field which is unique within the records of each database.
Each row contains one key value pair and you need to use the ``record_id``
to connect the related key-value pairs. You can easily do this with ``tidyr``
and ``dplyr`` in R or ``pandas`` in Python. An example to query the pathway
annotations from SignaLink:

    https://omnipathdb.org/annotations?databases=SignaLink_pathway

Or the tissue expression of BMP7 from Human Protein Atlas:

    https://omnipathdb.org/annotations?databases=HPA_tissue&proteins=BMP7


Roles in inter-cellular communication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another query type is the ``intercell``, providing information about the
roles in inter-cellular signaling. E.g. if a protein is a ligand, a receptor,
an extracellular matrix (ECM) component, etc. The proteins and protein
complexes are classified into categories. The categories are defined by a
number of attributes:

* `aspect`: funtional (e.g. ion channel) or locational (e.g. plasma
  membrane transmembrane).
* `scope`: generic (e.g. ligand) or specific (e.g. interleukin)
* `source`: resource specific (from one resource) or composite (combined
  from more resources)
* `causality`: transmitter (delivering signal from the expressing cell)
  or receiver (receiving signal into the expressing cell) or both
* `topology`: major localization categories derived from the locational
  categories: plasma membrane transmembrane or peripheral or secreted

The `intercell` database defines 25 functional and 10 locational generic,
composite categories. The number of specific categories is above 1,000.

You can use all these attributes in your queries, see the exact keys and
values at https://omnipathdb.org/queries/intercell

Some example queries:

    https://omnipathdb.org/intercell?proteins=EGFR,ULK1,ATG4A,BMP8B

All the resource specific functional classes for one protein:

    https://omnipathdb.org/intercell?source=resource_specific&aspect=functional&proteins=P00533

A list of all ECM proteins:

    https://omnipathdb.org/intercell?categories=ecm


Exploring possible parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes the names and values of the query parameters are not intuitive,
even though in many cases the server accepts multiple alternatives. To see
the possible parameters with all possible values you can use the ``queries``
query type. The server checks the parameter names and values exactly against
these rules and if any of them don't match you will get an error message
instead of reply. To see the parameters for the ``interactions`` query:

    https://omnipathdb.org/queries/interactions


Can I use OmniPath in R?
========================

You can download the data from the webservice and load into R. Thanks to
our colleague Attila Gabor we have a dedicated package for this:

    https://github.com/saezlab/OmnipathR


Installation
============

**Warning:** ``pip install pypath`` installs another package, you find
``pypath`` in PyPI under the name ``pypath-omnipath``:

.. code:: bash

    pip install pypath-omnipath

Linux
-----

In almost any up-to-date Linux distribution the dependencies of **pypath** are
built-in, or provided by the distributors. You can simply install **pypath**
by **pip** (see below).
If any non mandatory dependency is still missing, you can install them the
usual way by *pip* or your package manager.

igraph C library, cairo and pycairo
-----------------------------------

For the legacy network class or the ``igraph`` conversion from the current
network class *python-igraph* must be installed.
*python(2)-igraph* is a Python interface to use the igraph C library. The
C library must be installed. The same goes for *cairo*, *py(2)cairo* and
*graphviz*.


Directly from git
-----------------

.. code:: bash

    pip install git+https://github.com/saezlab/pypath.git

With pip
--------

Download the package from /dist, and install with pip:

.. code:: bash

    pip install pypath-x.y.z.tar.gz

Build source distribution
-------------------------

Clone the git repo, and run setup.py:

.. code:: bash

    python setup.py sdist

Mac OS X
--------

Recently the installation on Mac should not be more complicated than on Linux:
you can simply install by **pip** (see above).

When ``igraph`` was a mandatory dependency and it didn't provide wheels
the OS X installation was not straightforward primarily because cairo needs to
be compiled from source. If you want igraph and cairo we provide two scripts
in ``scripts``: the **mac-install-brew.sh** installs everything with HomeBrew and
**mac-install-conda.sh** installs from Anaconda distribution. With these
scripts, installation of igraph, cairo and graphviz goes smoothly most of the
time and options are available to omit the last two. To know more, see
the description in the script header. There is a third script
**mac-install-source.sh** which compiles everything from source and presumes
only Python 2.7 and Xcode installed. We do not recommend this as it is time
consuming and troubleshooting requires expertise.

Troubleshooting
^^^^^^^^^^^^^^^

* ``no module named ...`` when you try to load a module in Python. Did
  the installation of the module run without error? Try to run again the specific
  part from the mac install shell script to see if any error comes up. Is the
  path where the module has been installed in your ``$PYTHONPATH``? Try ``echo
  $PYTHONPATH`` to see the current paths. Add your local install directories if
  those are not there, e.g.
  ``export PYTHONPATH="/Users/me/local/python2.7/site-packages:$PYTHONPATH"``.
  If it works afterwards, don't forget to append these export path statements to
  your ``~/.bash_profile``, so these will be set every time you launch a new
  shell.

* ``pkgconfig`` not found. Check if the ``$PKG_CONFIG_PATH`` variable is
  set correctly, and pointing on a directory where pkgconfig really can be
  found.

* Error while trying to install py(2)cairo by pip. py(2)cairo could not be
  installed by pip, but only by waf. Please set the ``$PKG_CONFIG_PATH`` before.
  See **mac-install-source.sh** on how to install with waf.

* Error at pygraphviz build: ``graphviz/cgraph.h file not found``. This is
  because the directory of graphviz detected wrong by pkgconfig. See
  **mac-install-source.sh** how to set include dirs and library dirs by
  ``--global-option`` parameters.

* Can not install bioservices, because installation of jurko-suds fails. Ok,
  this fails because pip is not able to install the recent version of
  setuptools, because a very old version present in the system path. The
  development version of jurko-suds does not require setuptools, so you can
  install it directly from git as it is done in **mac-install-source.sh**.

* In **Anaconda**, *pypath* can be imported, but the modules and classes are
  missing. Apparently Anaconda has some built-in stuff called *pypath*. This
  has nothing to do with this module. Please be aware that Anaconda installs a
  completely separated Python distribution, and does not detect modules in the
  main Python installation. You need to install all modules within Anaconda's
  directory. **mac-install-conda.sh** does exactly this. If you still
  experience issues, please contact us.

* ``error: openssl/ssl.h: No such file or directory``: In order to install
  the Python modules ``pyopenssl`` and its dependency ``cryptography`` on
  some systems the development headers of OpenSSL need to be available. This
  is not the case if you can install ``pyopenssl`` from a wheel. If you get
  an error about a missing libssl header, just install the appropriate
  packages, in Debian based distros these are ``libssl-dev`` and
  ``libffi-dev``, in Red Hat based distros ``openssl-devel`` and
  ``libffi-devel``. In Mac OS X install ``openssl`` by ``homebrew``.

Microsoft Windows
-----------------

Not many people have used *pypath* on Microsoft computers so far. Please share
your experiences and contact us if you encounter any issue. We appreciate
your feedback, and it would be nice to have better support for other computer
systems.

With Anaconda
^^^^^^^^^^^^^

The same workflow like you see in ``mac-install-conda.sh`` should work for
Anaconda on Windows. The only problem you certainly will encounter is that not
all the channels have packages for all platforms. If certain channel provides
no package for Windows, or for your Python version, you just need to find an
other one. For this, do a search:

.. code:: bash

    anaconda search -t conda <package name>

For example, if you search for *pycairo*, you will find out that *vgauther*
provides it for osx-64, but only for Python 3.4, while *richlewis* provides
also for Python 3.5. And for win-64 platform, there is the channel of
*KristanAmstrong*. Go along all the commands in ``mac-install-conda.sh``, and
modify the channel if necessary, until all packages install successfully.

With other Python distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here the basic principles are the same as everywhere: first try to install all
external dependencies, after *pip* install should work. On Windows certain
packages can not be installed by compiled from source by *pip*, instead the
easiest to install them precompiled. These are in our case *fisher, lxml,
numpy (mkl version), pycairo, igraph, pygraphviz, scipy and statsmodels*. The
precompiled packages are available `here <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`__.
We tested the setup with Python 3.4.3 and Python 2.7.11. The former should just
work fine, while with the latter we have issues to be resolved.

Known issues
^^^^^^^^^^^^

* *"No module fabric available."* -- or *pysftp* missing: this is not
  important, only certain data download methods rely on these modules, but
  likely you won't call those at all.
* Progress indicator floods terminal: sorry about that, will be fixed soon.
* Encoding related exceptions in Python2: these might occur at some points in
  the module, please send the traceback if you encounter one, and we will fix
  as soon as possible.
* For Mac OS X (v >= 10.11 El Capitan) import of pypath fails with error:
  "libcurl link-time ssl backend (openssl) is different from compile-time ssl
  backend (none/other)". To fix it, you may need to reinstall pycurl library
  using special flags. More information and steps can be found
  `here <https://cscheng.info/2018/01/26/installing-pycurl-on-macos-high-sierra.html>`_.

*Special thanks to Jorge Ferreira for testing pypath on Windows!*

Release History
===============

Main improvements in the past releases:

0.1.0
-----

* First release of PyPath, for initial testing.

0.2.0
-----

* Lots of small improvements in almost every module
* Networks can be read from local files, remote files, lists or provided by any function
* Almost all redistributed data have been removed, every source downloaded from the original provider.

0.3.0
-----

* First version with partial Python 3 support.

0.4.0
-----

* **pyreact** module with **BioPaxReader** and **PyReact** classes added
* Process description databases, BioPax and PathwayCommons SIF conversion rules are supported
* Format definitions for 6 process description databases included.

0.5.0
-----

* Many classes have been added to the **plot** module
* All figures and tables in the manuscript can be generated automatically
* This is supported by a new module, **analysis**, which implements a generic workflow in its **Workflow** class.

0.5.32
------

* `chembl`, `unichem`, `mysql` and `mysql_connect` modules made Python3 compatible

0.6.31
------

* Orthology translation of network
* Homologene UniProt dict to translate between different organisms UniProt-to-UniProt
* Orthology translation of PTMs
* Better processing of PhosphoSite regulatory sites

0.7.0
-----

* TF-target, miRNA-mRNA and TF-miRNA interactions from many databases

0.7.74
------

* New web server based on `pandas` data frames
* New module `export` for generating data frames of interactions or enzyme-substrate interactions
* New module `websrvtab` for exporting data frames for the web server
* TF-target interactions from DoRothEA

0.7.93
------

* New `dataio` methods for Gene Ontology

0.7.110
-------

* Many new docstrings


0.8
---

* New module `complex`: a comprehensive database of complexes
* New module `annot`: database of protein annotations (function, location)
* New module `intercell`: special methods for data integration focusing on intercellular communication
* New module `bel`: BEL integration
* Module `go` and all the connected `dataio` methods have been rewritten offering a workaround for
  data access despite GO's terrible web services and providing much more versatile query methods
* Removed MySQL support (e.g. loading mapping tables from MySQL)
* Modules `mapping`, `reflists`, `complex`, `ptm`, `annot`, `go` became services:
  these modules build databases and provide query methods, sometimes they even automatically
  delete data to free memory
* New interaction category in `data_formats`: `ligand_receptor`
* Improved logging and control over verbosity
* Better control over parameters by the `settings` module
* Many methods in `dataio` have been improved or fixed, docs and code style largely improved
* Started to add tests especially for methods in `dataio`

0.9
---
* The network database is not dependent any more on `python-igraph` hence it
  has been removed from the mandatory dependencies
* New API for the network, interactions, evidences, molecular entities

0.10.0
------
* New module structure: modules grouped into `core`, `inputs`, `internals`,
  `legacy`, `omnipath`, `resources`, `share` and `utils` submodules.

0.11.0
------
* Redesign of the intercell (intercellular communication roles) database

Upcoming
--------

* New, more flexible network reader class
* Full support for multi-species molecular interaction networks
  (e.g. pathogene-host)
* Better support for not protein only molecular interaction networks
  (metabolites, drug compounds, RNA)

Features
========

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
of gene, protein and miRNA ID types. It has the feature to translate
secondary UniProt ACs to primaries, and Trembl ACs to SwissProt, using
primary Gene Symbols to find the connections. This module automatically
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

The ``pypath.utils.homology`` module is able to find the homologues of genes
between two organisms. It uses data from NCBI HomoloGene, and soon we will
extend it to use Ensembl or UniProt as alternatives. This module is really
simple to use:

.. code:: python

    from pypath.utils import homology
    homology.translate('P00533', 10090) # translating the human EGFR to mouse
    # ['Q01279'] # it returns the mouse Egfr UniProt AC
