Release History
------------------
This is a summary of the changelog.

0.1.0:
+++++++++++
* First release of PyPath, for initial testing.

0.2.0:
+++++++++++
* Lots of small improvements in almost every module
* Networks can be read from local files, remote files, lists or provided by any function
* Almost all redistributed data have been removed, every source downloaded from the original provider.

0.3.0:
+++++++++++
* First version whith partial Python 3 support.

0.4.0:
+++++++++++
* **pyreact** module with **BioPaxReader** and **PyReact** classes added
* Process description databases, BioPax and PathwayCommons SIF conversion rules are supported
* Format definitions for 6 process description databases included.

0.5.0:
+++++++++++
* Many classes have been added to the **plot** module
* All figures and tables in the manuscript can be generated automatically
* This is supported by a new module, **analysis**, which implements a generic workflow in its **Workflow** class.

0.5.32:
+++++++++++
* `chembl`, `unichem`, `mysql` and `mysql_connect` modules made Python3 compatible

0.6.31:
+++++++++++
* Orthology translation of network
* Homologene UniProt dict to translate between different organisms UniProt-to-UniProt
* Orthology translation of PTMs
* Better processing of PhosphoSite regulatory sites

Upcoming:
+++++++++++
* New, more flexible network reader class
* Full support for multi-species molecular interaction networks (e.g. pathogene-host)
* Better support for not protein only molecular interaction networks (metabolites, drug compounds, RNA)
* Silent mode: a way to suppress messages and progress bars
