.. bioigraph documentation master file, created by
   sphinx-quickstart2 on Mon Mar  2 12:43:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to bioigraph's documentation!
=====================================

Contents:

.. toctree::
   :maxdepth: 2

   index

Examples
========

Example 1: building the network
+++++++++++++++++++++++++++++++

Import the module and create an instance. MySQL settings are stored by default in mysql_config/defaults.mysql, now you don't need to modify it:

.. code-block:: python

    import bioigraph
    mysql_gelati = (None, 'mapping_gelati')
    mysql_chembl = (None, 'chembl_ebi')
    net = bioigraph.BioGraph(9606, mysql = mysql_gelati)

The tables for ID conversions, and the network of the selected manually curated pathway databases can be initialized this way:

.. code-block:: python

    net.init_network()

Resources without references for each interaction are available separately. For example, to build a comprehensive phosphonetwork, you might want to load interactions from high-throughput data in PhosphoSite, kinase-substrate relationships from PhosphoNetworks and MIMP. This is the way to do it:

.. code-block:: python

    from bioigraph.data_formats import good
    net.load_resources(lst = {'mimp': good['mimp']})
    net.load_resources(lst = {'pnetworks': good['pnetworks']})
    net.load_resources(lst = {'psite_noref': good['psite_noref']})

In case you have your own interaction data, and you wish to merge this into the network, you first need to define its format, and then call the load_resource() function. In this simple example, the proteins are noted by Gene Symbols, their names are in columns #0 and #1, while from column #2 an edge attribute named 'score' will be read. The file is tab separated, and the file name is 'mylist.sif'. (This won't work unless you really have a file like this.)

.. code-block:: python

    mylist = bioigraph.input_formats.ReadSettings(name = "mylist", 
        separator = "\t", nameColA = 0, nameColB = 1,
        nameTypeA = "genesymbol", nameTypeB = "genesymbol",
        isDirected = False, inFile = 'mylist.sif', 
        extraEdgeAttrs = {
            'score': 2
        })
    net.load_resource(mylist)


Although proteins are identified by UniProt IDs, it is useful to have standard gene symbols as labels, to make things more human understandable. 

.. code-block:: python

    net.genesymbol_labels()

Now the igraph object is in net.graph, vertex and edge attributes are accessible the usual way:

.. code-block:: python

    net.graph.vs[0]['name']
    net.graph.vs[0]['label']

Let's see the default attributes of edges:

.. code-block:: python

    net.graph.es[0]['sources']
    net.graph.es[0]['references']
    net.graph.es[0]['dirs']

The graph is undirected by default, and 'Direction' objects can be found in <edge>['dirs'] to describe the directions and signs of the edge:

.. code-block:: python

    print net.graph.es[111]['dirs']
    net.graph.es[111].is_directed()
    net.graph.es[111]['dirs'].is_inhibition()
    net.graph.es[111]['dirs'].which_dirs()
    net.graph.es[111]['dirs'].get_dir(net.graph.es[111]['dirs'].which_dirs()[0])

To convert the igraph object to a directed graph, use this function:

.. code-block:: python

    net.get_directed()
    net.dgraph.ecount()
    net.dgraph.vcount()

To replace the graph with the giant component of the directed graph:

.. code-block:: python

    net.graph = net.get_giant(graph = net.dgraph)

Get compounds for all the proteins in the network:

.. code-block:: python

    net.compounds_from_chembl((None, 'chembl_ebi'), pchembl = True)

This works only from EBI, as the ChEMBL MySQL instance running on an EBI virtual machine. To use it from home, start an SSH tunnel to EBI: ssh -o proxycommand="ssh -p 2222 gate.ebi.ac.uk proxy %h" login.ebi.ac.uk -L 3306:172.22.68.151:3306
And while the tunnel is open, run the command above with mysql settings pointing to localhost:

.. code-block:: python

    net.compounds_from_chembl((None, 'chembl'), pchembl = True)

You can find the list of compounds and the requested data (now only the pchembl values) in vertex attributes:

.. code-block:: python

    net.graph.vs[2]['compounds_chembl']
    net.graph.vs[2]['compounds_data']
    net.graph.vs[2]['compounds_data'][0]['pchembl']

Loading PTMs. This function loads PTMs from all the available resources. It takes long time first, but once it saves the downloaded data under ./cache directory, it is much faster to run again.

.. code-block:: python

    net.load_ptms()

Individual PTM resources can be loaded the following way:

.. code-block:: python

    net.load_phospho_dmi(source = 'phosphoELM')

In the latter case, this function is needed to merge the identical PTMs loaded from multiple resources:

.. code-block:: python

    net.uniq_ptms()

PTMs are stored in objects. This example shows how to access the type of modification, and the number and name of the modified residue.

.. code-block:: python

    net.graph.es[70]['ptm'][0].ptm.typ
    net.graph.es[70]['ptm'][0].ptm.protein
    net.graph.es[70]['ptm'][0].ptm.residue.name
    net.graph.es[70]['ptm'][0].ptm.residue.number

Example 2: using the Mapper class for translating IDs
+++++++++++++++++++++++++++++++++++++++++++++++++++++

The mapping submodule of bioigraph can be used for ID conversion. Here is a basic example:

.. code-block:: python

    from bioigraph import mapping
    m = mapping.Mapper(9606, mysql_conf = (None, 'mapping_gelati'))
    m.load_mappings()
    result = {}
    gene_list = ['EGFR', 'AKT1', 'GABARAPL1', 'TP53']
    for g in gene_list:
        result[g] = m.map_name(g, 'genesymbol', 'uniprot')

Reference
=========

.. automodule:: bioigraph.bioigraph
   :members:

.. autoclass:: bioigraph.bioigraph.BioGraph
   :members:

.. automodule:: bioigraph.intera
   :members:

.. autoclass:: bioigraph.intera.Residue
   :members:

.. autoclass:: bioigraph.intera.Mutation
   :members:

.. autoclass:: bioigraph.intera.Motif
   :members:

.. autoclass:: bioigraph.intera.Ptm
   :members:

.. autoclass:: bioigraph.intera.Domain
   :members:

.. autoclass:: bioigraph.intera.DomainMotif
   :members:

.. autoclass:: bioigraph.intera.DomainDomain
   :members:

.. autoclass:: bioigraph.intera.Interface
   :members:

.. automodule:: bioigraph.mapping
   :members:

.. autoclass:: bioigraph.mapping.Mapper
   :members:

.. automodule:: bioigraph.moduland
   :members:

.. automodule:: bioigraph.gdsc
   :members:

.. automodule:: bioigraph.dataio
   :members:

.. automodule:: bioigraph.chembl
   :members:

.. automodule:: bioigraph.unichem
   :members:

.. automodule:: bioigraph.seq
   :members:

.. automodule:: bioigraph.mysql
   :members:

.. automodule:: bioigraph.mysql_connect
   :members:

.. automodule:: bioigraph.seq
   :members: