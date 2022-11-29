############################################################################
*pypath:* A Python module for molecular signaling prior knowledge processing
############################################################################

    **Important:** New module structure and new network API

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

:Py2/3: Although we still keep the compatibility with Python 2, we don't
        test ``pypath`` in this environment and very few people uses it
        already. We highly recommend to use ``pypath`` in Python 3.6+.

:documentation: http://saezlab.github.io/pypath
:issues: https://github.com/saezlab/pypath/issues
:contact: omnipathdb@gmail.com
:developers: ``pypath`` is developed in the Saez Lab (http://saezlab.org) by
  Olga Ivanova, Nicolàs Palacio and Dénes Türei; the R package and the
  Cytoscape app are developed and maintained by Francesco Ceccarelli, Attila
  Gábor, Alberto Valdeolivas and Nicolàs Palacio.

.. toctree::
    :maxdepth: 5
    :caption: Contents:

    installation
    webservice
    releasehistory

**pypath** is a Python module for processing molecular biology data resources,
combining them into databases and providing a versatile interface in Python
as well as exporting the data for access through other platforms such as
the R (the OmnipathR R/Bioconductor package), web service (at
http://omnipathdb.org), Cytoscape (the OmniPath Cytoscape app) and BEL
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


*********
Reference
*********

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:

   pypath
