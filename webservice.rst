Webservice
==========

I run and update time to time the pypath webservice on my virtual machine on EBI infrastructure, serving the OmniPath data with additional kinase-substrate interactions and PTMs. The current IP address of the VM is 172.22.71.26, but this might eventually change (given by the Systems, every time when I reboot the VM after upgrades). The webservice is set up to listen on port 33333. It serves data in REST style, by HTTP protocol (browser, wget, curl or anything can make requests). This host is accessible only on the EBI internal network. From outside you can use SSH tunnel of course. The webservice currently recognizes 6 types of queries: ``interactions``, ``ptms``, ``resources``, ``network``, ``about`` and ``info``. 

A request without any parameter, gives some basic numbers about the actual loaded dataset:

    http://omnipathdb.org

The ``about`` tells the version number:

    http://omnipathdb.org/about

The ``network`` prints basic statistics about the whole network:
    
    http://omnipathdb.org/network

The ``resources`` returns the list of all resources with their size:
    
    http://omnipathdb.org/network

The ``info`` returns a HTML page with comprehensive information about the resources:

    http://omnipathdb.org/info

The ``interactions`` accepts some parameters and returns interactions in tabular format. This example returns all interactions of EGFR (P00533), with sources and references listed:

    http://omnipathdb.org/interactions/P00533/?fields=sources&fields=references

The parameters can be omitted. More UniProts can be given separated by comma, and JSON format is available too (better for import to Python!):

    http://omnipathdb.org/interactions/P00533,O15117,Q96FE5?format=json

Another interface is ``ptms``, to list enzymes, substrates and PTMs. 

    http://omnipathdb.org/ptms/P00533?ptm_type=phosphorylation&fields=sources&fields=references

To list all interactions simply request:

    http://omnipathdb.org/interactions

To list all PTMs similarly:

    http://omnipathdb.org/ptms