from __future__ import annotations

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls

def gutmgene_human_raw() -> list[tuple]:
    """
    Retrieves raw data about gut microbiota genes in humans from the gut microbiota gene database.

    Returns:
    A list of named tuples containing information about gut microbiota genes in human.
    """
    HumanFields = collections.namedtuple(
        'gutMGeneHuman',
        [
            'gut_microbiota',
            'gut_microbiata_tax_id',
            'gut_microbiota_id',
            'classification',
            'gene_symbol',
            'gene_entrez_id',
            'alteration',
            'high_low_throughput',
        ],
        defaults = None
    )
    
    url = urls.urls['gutmgene']['url_human']
    c = curl.Curl(url, silent = False, large = True)
    
    result = set()
        
    for l in c.result:
        l = l.strip().split('\t')
        l = [None if not i else i for i in l]
        if l:
            result.add(HumanFields(*l))
            
    return list(result)


def gutmgene_mouse_raw() -> list[tuple]:
    """
    Retrieves raw data about gut microbiota genes in mice from the gut microbiota gene database.

    Returns:
    A list of named tuples containing information about gut microbiota genes in mouse.
    """

    MouseFields = collections.namedtuple(
        'gutMGeneMouse',
        [
            'gut_microbiota',
            'gut_microbiata_tax_id',
            'gut_microbiota_id',
            'classification',
            'gene_symbol',
            'gene_entrez_id',
            'alteration',
            'high_low_throughput',
        ],
        defaults = None
    )
    
    url = urls.urls['gutmgene']['url_mouse']
    c = curl.Curl(url, silent = False, large = True)
    
    result = set()
        
    for l in c.result:
        if l.startswith('"'):
            continue
        
        l = l.replace('"', '')
        l = l.strip().split('\t')
        l = [None if not i else i for i in l]
        if l:
            result.add(MouseFields(*l))
            
    return list(result)
            
        