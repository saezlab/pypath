***************
Release history
***************

Main improvements in the past releases:


0.1.0
=====

* First release of pypath, for initial testing.


0.2.0
=====

* Lots of small improvements in almost every module
* Networks can be read from local files, remote files, lists or provided by
  any function
* Almost all redistributed data have been removed, every source downloaded
  from the original provider.


0.3.0
=====

* First version with partial Python 3 support.


0.4.0
=====

* **pyreact** module with **BioPaxReader** and **PyReact** classes added
* Process description databases, BioPax and PathwayCommons SIF conversion
  rules are supported
* Format definitions for 6 process description databases included.


0.5.0
=====

* Many classes have been added to the **plot** module
* All figures and tables in the manuscript can be generated automatically
* This is supported by a new module, **analysis**, which implements a generic
workflow in its **Workflow** class.


0.7.74
======

* **homology** module: finds the homologs of proteins using the NCBI
Homologene database and the homologs of PTM sites using UniProt sequences
and PhosphoSitePlus homology table
* **ptm** module: fully integrated way of processing enzyme-substrate
interactions from many databases and their translation by homology to other
species
* **export** module: creates ``pandas.DataFrame`` or exports the network into
tabular file
* New webservice
* TF Regulons database included and provides much more comprehensive
transcriptional regulation resources, including literature curated, in silico
predicted, ChIP-Seq and expression pattern based approaches
* Many network resources added, including miRNA-mRNA and TF-miRNA interactions


Upcoming
========

* New, more flexible network reader class
* Full support for multi-species molecular interaction networks
(e.g. pathogene-host)
* Better support for not protein only molecular interaction networks
(metabolites, drug compounds, RNA)
* ChEMBL webservice interface, interface for PubChem and eventually
forDrugBank
* Silent mode: a way to suppress messages and progress bars
