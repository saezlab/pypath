.. pypath documentation master file, created by
   sphinx-quickstart2 on Mon Mar  2 12:43:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pypath's documentation!
=====================================

pypath is a Python module for cellular signaling pathways analysis. It is accompanied and developed together with a high confidence, literature curated signaling network, OmniPath_.

.. _OmniPath: http://omnipathdb.org/

Source code
+++++++++++

pypath is free software, licensed under GPLv3_.
The code is available at github_.
Recent package for PIP is available at http://pypath.omnipathdb.org/releases/latest/, archives are available at http://pypath.omnipathdb.org/releases/archive/.

.. _GPLv3: http://www.gnu.org/licenses/gpl.html
.. _github: https://github.com/saezlab/pypath

Main features
+++++++++++++

- Undirected or directed networks.
- Easy, often seemless protein ID conversion.
- Efficient handling of annotations, especially sources, literature references, directions, effect signs (stimulation/inhibition) and enzyme-substrate interactions.
- Ready integration of dozens of bioinformatics resources, most of them being downloaded on the fly, directly from the original source.
- Caching: after downloading data, a local copy saved, so faster and offline run is possible afterwards.
- Based on igraph_: plenty of graph methods are available with excellent computational performance.
- Partial support for non-human species.
- Partial support for other molecular species than proteins.  

.. _igraph: http://igraph.org/

.. toctree:
   :maxdepth: 3

Indices and tables
==================

   * :ref:`genindex`
   * :ref:`modindex`

Examples
========

Example 1: building the network
+++++++++++++++++++++++++++++++

Import the module and create an instance:

.. code-block:: python

    import pypath
    pa = pypath.PyPath()

The tables for ID conversions, and the network of the selected manually curated pathway databases (OmniPath_) can be initialized this way:

.. _OmniPath: http://omnipathdb.org/

.. code-block:: python

    pa.init_network()

Resources without references for each interaction are available separately. For example, to build a comprehensive phosphonetwork, you might want to load interactions from high-throughput data in PhosphoSite, kinase-substrate relationships from PhosphoNetworks and MIMP. This is the way to do it:

.. code-block:: python

    from pypath.data_formats import ptm_misc
    pa.load_resources(lst = {'mimp': ptm_misc['mimp']})
    pa.load_resources(lst = {'pnetworks': ptm_misc['pnetworks']})
    pa.load_resources(lst = {'psite_noref': ptm_misc['psite_noref']})

In case you have your own interaction data, and you wish to merge this into the network, you first need to define its format, and then call the load_resource() function. In this simple example, the proteins are noted by Gene Symbols, their names are in columns #0 and #1, while from column #2 an edge attribute named 'score' will be read. The file is tab separated, and the file name is 'mylist.sif'. (This won't work unless you really have a file like this.)

.. code-block:: python

    mylist = pypath.input_formats.ReadSettings(name = "mylist", 
        separator = "\t", nameColA = 0, nameColB = 1,
        nameTypeA = "genesymbol", nameTypeB = "genesymbol",
        isDirected = False, inFile = 'mylist.sif', 
        extraEdgeAttrs = {
            'score': 2
        })
    pa.load_resource(mylist)


Although proteins are identified by UniProt IDs, it is useful to have standard gene symbols as labels, to make things more human understandable. 

.. code-block:: python

    pa.genesymbol_labels()

Now the igraph object is in pa.graph, vertex and edge attributes are accessible the usual way:

.. code-block:: python

    pa.graph.vs[0]['name']
    pa.graph.vs[0]['label']

Let's see the default attributes of edges:

.. code-block:: python

    pa.graph.es[0]['sources']
    pa.graph.es[0]['references']
    pa.graph.es[0]['dirs']

The graph is undirected by default, and 'Direction' objects can be found in <edge>['dirs'] to describe the directions and signs of the edge:

.. code-block:: python

    print pa.graph.es[111]['dirs']
    pa.graph.es[111].is_directed()
    pa.graph.es[111]['dirs'].is_inhibition()
    pa.graph.es[111]['dirs'].which_dirs()
    pa.graph.es[111]['dirs'].get_dir(pa.graph.es[111]['dirs'].which_dirs()[0])

To convert the igraph object to a directed graph, use this function:

.. code-block:: python

    pa.get_directed()
    pa.dgraph.ecount()
    pa.dgraph.vcount()

To replace the graph with the giant component of the directed graph:

.. code-block:: python

    pa.graph = pa.get_giant(graph = pa.dgraph)

Get compounds for all the proteins in the network:

.. code-block:: python

    pa.compounds_from_chembl((None, 'chembl_ebi'), pchembl = True)

This works only if you have a MySQL server hosting an instance of ChEMBL:

.. code-block:: python

    pa.compounds_from_chembl((None, 'chembl'), pchembl = True)

You can find the list of compounds and the requested data (now only the pchembl values) in vertex attributes:

.. code-block:: python

    pa.graph.vs[2]['compounds_chembl']
    pa.graph.vs[2]['compounds_data']
    pa.graph.vs[2]['compounds_data'][0]['pchembl']

Loading PTMs. This function loads PTMs from all the available resources. It takes long time first, but once it saves the downloaded data under `./cache` directory, it is much faster to run again.

.. code-block:: python

    pa.load_ptms()

Individual PTM resources can be loaded the following way:

.. code-block:: python

    pa.load_phospho_dmi(source = 'phosphoELM')

In the latter case, this function is needed to merge the identical PTMs loaded from multiple resources:

.. code-block:: python

    pa.uniq_ptms()

PTMs are stored in objects. This example shows how to access the type of modification, and the number and name of the modified residue.

.. code-block:: python

    pa.graph.es[70]['ptm'][0].ptm.typ
    pa.graph.es[70]['ptm'][0].ptm.protein
    pa.graph.es[70]['ptm'][0].ptm.residue.name
    pa.graph.es[70]['ptm'][0].ptm.residue.number

Example 2: using the Mapper class for translating IDs
+++++++++++++++++++++++++++++++++++++++++++++++++++++

The mapping submodule of pypath can be used for ID conversion. Here is a basic example:

.. code-block:: python

    from pypath import mapping
    m = mapping.Mapper(9606, mysql_conf = (None, 'mapping_gelati'))
    m.load_mappings()
    result = {}
    gene_list = ['EGFR', 'AKT1', 'GABARAPL1', 'TP53']
    for g in gene_list:
        result[g] = m.map_name(g, 'genesymbol', 'uniprot')

Reference
=========

.. automodule:: pypath.pypath
   :members:

.. autoclass:: pypath.pypath.PyPath
   :members:

.. autoclass:: pypath.pypath.Direction
   :members:

.. automodule:: pypath.intera
   :members:

.. autoclass:: pypath.intera.Residue
   :members:

.. autoclass:: pypath.intera.Mutation
   :members:

.. autoclass:: pypath.intera.Motif
   :members:

.. autoclass:: pypath.intera.Ptm
   :members:

.. autoclass:: pypath.intera.Domain
   :members:

.. autoclass:: pypath.intera.DomainMotif
   :members:

.. autoclass:: pypath.intera.DomainDomain
   :members:

.. autoclass:: pypath.intera.Interface
   :members:

.. automodule:: pypath.mapping
   :members:

.. autoclass:: pypath.mapping.Mapper
   :members:

.. automodule:: pypath.moduland
   :members:

.. automodule:: pypath.dataio
   :members:

.. automodule:: pypath.chembl
   :members:

.. automodule:: pypath.unichem
   :members:

.. automodule:: pypath.seq
   :members:

.. automodule:: pypath.mysql
   :members:

.. automodule:: pypath.mysql_connect
   :members:
