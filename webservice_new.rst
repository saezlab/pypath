Webservice
==========

One instance of the pypath webservice runs at the domain http://omnipathdb.org/, serving not only the OmniPath data but other datasets: TF-target interactions from TF Regulons, a large collection additional enzyme-substrate interactions, and literature curated miRNA-mRNA interacions combined from 4 databases. The webservice implements a very simple REST style API, you can make requests by HTTP protocol (browser, wget, curl or whatever). The webservice currently recognizes 6 types of queries: ``interactions``, ``ptms``, ``resources``, ``network``, ``about`` and ``info``.

Mouse and rat
-------------
Except the miRNA interactions all interactions are available for human, mouse and rat. The rodent data has been translated from human using the NCBI Homologene database. Many human proteins have no known homolog in rodents hence rodent datasets are smaller than their human counterparts. Note, if you work with mouse omics data you might do better to translate your dataset to human (for example using the ``pypath.homology`` module) and use human interaction data.


A request without any parameter, gives some basic numbers about the actual loaded dataset:

    http://omnipathdb.org

The ``about`` tells the version number:

    http://omnipathdb.org/about

The ``network`` prints basic statistics about the whole network:
    
    http://omnipathdb.org/network

The ``resources`` returns the list of all resources with their size:
    
    http://omnipathdb.org/resources

The ``info`` returns a HTML page with comprehensive information about the resources:

    http://omnipathdb.org/info

The ``interactions`` accepts some parameters and returns interactions in tabular format. This example returns all interactions of EGFR (P00533), with sources and references listed.

    http://omnipathdb.org/interactions/?partners=P00533&fields=sources,references

By default only the OmniPath dataset used, to query the TF Regulons or add the extra enzyme-substrate interactions you need to set additional parameters. For example to query the transcriptional regulators of EGFR:

    http://omnipathdb.org/interactions/?targets=EGFR&types=TF

The TF Regulons database assigns confidence levels to the interactions. You might want to select only the highest confidence, *A* category:

    http://omnipathdb.org/interactions/?targets=EGFR&types=TF&tfregulons_levels=A

Show the transcriptional targets of Smad2 homology translated to rat including the confidence levels from TF Regulons:

    http://localhost:33333/interactions/?genesymbols=1&fields=type,ncbi_tax_id,tfregulons_level&organisms=10116&sources=Smad2&types=TF

Query interactions from PhosphoNetworks which is part of the *kinaseextra* dataset:

    http://localhost:33333/interactions/?genesymbols=1&fields=sources&databases=PhosphoNetworks&datasets=kinaseextra

Get the interactions from Signor, SPIKE and SignaLink3:

    http://localhost:33333/interactions/?genesymbols=1&fields=sources,references&databases=Signor,SPIKE,SignaLink3

All interactions of MAP1LC3B:

    http://localhost:33333/interactions/?genesymbols=1&partners=MAP1LC3B

By default `partners` queries the interaction where either the source or the target is among the partners. If you set the `source_target` parameter to `AND` both the source and the target must be in the queried set:
    
    http://localhost:33333/interactions/?genesymbols=1&fields=sources,references&sources=ATG3,ATG7,ATG4B,SQSTM1&targets=MAP1LC3B,MAP1LC3A,MAP1LC3C,Q9H0R8,GABARAP,GABARAPL2&source_target=AND

As you see above you can use UniProt IDs and Gene Symbols in the queries and also mix them.
Get the miRNA regulating NOTCH1:
    
    http://localhost:33333/interactions/?genesymbols=1&fields=sources,references&datasets=mirnatarget&targets=NOTCH1

Note: with the exception of mandatory fields and genesymbols, the columns appear exactly in the order you provided in your query.

Another query type available is ``ptms`` which provides enzyme-substrate interactions. It is very similar to the ``interactions``:

    http://localhost:33333/ptms?genesymbols=1&fields=sources,references,isoforms&enzymes=FYN

Is there any ubiquitination reaction?

    http://localhost:33333/ptms?genesymbols=1&fields=sources,references&types=ubiquitination

And acetylation in mouse?

    http://localhost:33333/ptms?genesymbols=1&fields=sources,references&types=acetylation&organisms=10090

Rat interactions and 
