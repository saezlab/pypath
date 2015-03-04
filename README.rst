Bioigraph
#########


:note: Bioigraph supported only in Python 2.7, and no version 3.x is available currently.

:contributions: denes@ebi.ac.uk
:issues: https://github.com/saezlab/bioigraph/issues

**Bioigraph** is a Python package built around igraphthat to work with molecular network representations e.g. PPI, miRNA, drug compound interaction networks.

Installation
============

With pip
--------

Download the package from /dist, and install with pip:

.. code:: bash
    
    pip2 install bioigraph-x.y.z.tar.gz

Build source distribution
-------------------------

Clone the git repo, and run build.sh:

.. code:: bash
    
    cd path/
    ./build.sh

Features
========

The primary aim of **Bioigraph** is to build up networks from multiple sources on one igraph object. **Bioigraph** handles ambiguous ID conversion, reads custom edge and node attributes from text files and **MySQL**.

**Bioigraph** includes data and format definition from 19 high quality, literature curated databases. Descriptions and comprehensive information about the resources is available in the package. 

Submodules serves the handling of graph visualization, working with drug compound data, searching drug targets and compounds in **ChEMBL**. 

ID conversion module (bioigaph.mapping) can be used independently. It has the feature to translate secondary UniProt IDs to primaries, and Trembl IDs to Swissprot, using primary gene symbols to find the connections. 

**UniChem** submodule provides an interface to effectively query the UniChem service, use connectivity search with custom settings, and translate SMILEs to ChEMBL IDs with ChEMBL web service.

**ChEMBL** submodule queries directly your own ChEMBL MySQL instance, has the features to search targets and compounds from custom assay types and relationship types, to get activity values, binding domains, and action types. You need to download the ChEMBL MySQL dump, and load into your own server.

**MySQL** submodule helps to manage MySQL connections and track queries.

**ModuLand** module for running the ModuLand method family
