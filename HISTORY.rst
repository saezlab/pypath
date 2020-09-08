Release History
------------------
This is a summary of the changelog.

0.1.0:
^^^^^^^^^^^
* First release of PyPath, for initial testing.

0.2.0:
^^^^^^^^^^^
* Lots of small improvements in almost every module
* Networks can be read from local files, remote files, lists or provided
  by any function
* Almost all redistributed data have been removed, every source downloaded
  from the original provider.

0.3.0:
^^^^^^^^^^^
* First version with partial Python 3 support.

0.4.0:
^^^^^^^^^^^
* **pyreact** module with **BioPaxReader** and **PyReact** classes
* Process description databases, BioPax and PathwayCommons SIF conversion
  rules are supported
* Format definitions for 6 process description databases included.

0.5.0:
^^^^^^^^^^^
* Many new classes in the **plot** module
* All figures and tables in the manuscript can be generated automatically
* This is supported by a new module, **analysis**, which implements a
  generic workflow in its **Workflow** class.

0.5.32:
^^^^^^^^^^^
* `chembl`, `unichem`, `mysql` and `mysql_connect` modules are Python3
  compatible

0.6.31:
^^^^^^^^^^^
* Orthology translation of network
* Homologene UniProt dict to translate between different organisms
  UniProt-to-UniProt
* Orthology translation of PTMs
* Better processing of PhosphoSite regulatory sites

0.7.0
^^^^^^^^^^^
* TF-target, miRNA-mRNA and TF-miRNA interactions from many databases

0.7.74
^^^^^^^^^^^
* New web server based on `pandas` data frames
* New module `export` for generating data frames of interactions or
  enzyme-substrate interactions
* New module `websrvtab` for exporting data frames for the web server
* TF-target interactions from DoRothEA

0.7.93
^^^^^^^^^^^
* New `dataio` methods for Gene Ontology

0.7.110
^^^^^^^^^^^
* Many new docstrings

0.8
^^^^^^^^^^^
* New module `complex`: a comprehensive database of complexes
* New module `annot`: database of protein annotations (function, location)
* New module `intercell`: special methods for data integration focusing on
  intercellular communication
* New module `bel`: BEL integration
* Module `go` and all the connected `dataio` methods have been rewritten
  offering a workaround for data access despite GO's terrible web services
  and providing much more versatile query methods
* Removed MySQL support (e.g. loading mapping tables from MySQL)
* Modules `mapping`, `reflists`, `complex`, `ptm`, `annot`, `go` became
  services: these modules build databases and provide query methods,
  sometimes they even automatically delete data to free memory
* New interaction category in `data_formats`: `ligand_receptor`
* Improved logging and control over verbosity
* Better control over paremeters by the `settings` module
* Many methods in `dataio` have been improved or fixed, docs and code style
  largely improved
* Started to add tests especially for methods in `dataio`

0.9
^^^^^^^^^^^
* The network database is not dependent any more on `python-igraph` hence it
  has been removed from the mandatory dependencies
* New API for the network, interactions, evidences, molecular entities

0.10.0
^^^^^^^^^^^
* A complete reorganization of the module structure: submodules sorted into a
  few major groups: core, inputs, internals, omnipath, share, utils

0.11.0
^^^^^^^^^^^
* A complete redesign of the intercellular communication database
* License support

Upcoming:
^^^^^^^^^^^
* New, more flexible network reader class
* Full support for multi-species molecular interaction networks (e.g. pathogene-host)
* Better support for not protein only molecular interaction networks (metabolites, drug compounds, RNA)