from __future__ import annotations

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common

def drugcombdb_syner_antag_voting() -> list[tuple]:
    Fields = collections.namedtuple(
        'drugcombdb_syner_antag_voting',
        [
            'ID',
            'Drug1',
            'Drug2',
            'CellLine',
            'ZIP',
            'Bliss',
            'Loewe',
            'HSA',
            'ZIPclass',
            'Blissclass',
            'Loeweclass',
            'HSAclass',
            'synthetic',
            'classification',
        ],
        defaults = None
    )
    
    url = urls.urls['drugcombdb']['url_syner_antag_voting']
    c = curl.Curl(url, silent = False, large = True)
    
    result = set()
    
    for l in c.result:
        l = l.strip().split(',')
        l = [None if not i else i for i in l]
        if l and len(l) == len(Fields._fields):
            result.add(Fields(*l))
            
    return list(result)

def drugcombdb_drug_chemical_info() -> list[tuple]:
    Fields = collections.namedtuple(
        'drugcombdb_drug_chemical_info',
        [
            'DrugName',
            'cIDs',
            'drugNameOfficial',
            'molecularWeight',
            'smiles',
        ],
        defaults = None
    )
    
    url = urls.urls['drugcombdb']['url_drug_chemical_info']
    c = curl.Curl(url, silent = False, large = True, encoding = 'ISO-8859-1')
    
    result = set()
    
    for l in c.result:
        l = l.strip().split(',')
        l = [None if not i or i == '#N/A' else i for i in l]
        l = [None if not i else i for i in l]
        if l and len(l) == len(Fields._fields):
            result.add(Fields(*l))
            
    return list(result)


def drugcombdb_syndrugcomb_fda() -> list[tuple]:
    Fields = collections.namedtuple(
        'drugcombdb_syndrugcomb_fda',
        [
            'ID',
            'Drug1',
            'Drug2',
            'Machenism',
            'Source',
        ],
        defaults = None
    )
    
    url = urls.urls['drugcombdb']['url_syndrugcomb_fda']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='2drugs')
    
    result = set()
    
    for l in contents:
        l = [None if not i else i for i in l]
        if l and len(l) == len(Fields._fields):
            result.add(Fields(*l))
            
    return list(result)


def drugcombdb_syndrugcomb_textmining() -> list[tuple]:
    Fields = collections.namedtuple(
        'drugcombdb_syndrugcomb_textmining',
        [
            'Drug1',
            'Drug2',
            'Target',
            'Source',
        ],
        defaults = None
    )
    
    url = urls.urls['drugcombdb']['url_syndrugcomb_textmining']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='2drugs')
    
    result = set()
    
    for l in contents:
        l = l[:-1]
        l = [None if not i else i for i in l]
        if l:
            result.add(Fields(*l))
            
    return list(result)
    
    
def drugcombdb_syndrugcomb_external_synergism() -> list[tuple]:
    Fields = collections.namedtuple(
        'drugcombdb_syndrugcomb_external',
        [
            'Drug1',
            'Drug2',
            'PubmedID',
        ],
        defaults = None
    )
    
    url = urls.urls['drugcombdb']['url_syndrugcomb_external']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='ASDCD_synergism')
    
    result = set()
    
    for l in contents:
        l = l[:3]
        l = [None if not i else i for i in l]
        if l:
            result.add(Fields(*l))
            
    return list(result)


def drugcombdb_syndrugcomb_external_antagonism() -> list[tuple]:
    Fields = collections.namedtuple(
        'drugcombdb_syndrugcomb_external',
        [
            'Interaction_A',
            'Interaction_B',
        ],
        defaults = None
    )
    
    url = urls.urls['drugcombdb']['url_syndrugcomb_external']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='Drugbank_antagonism')
    
    result = set()
    
    for l in contents:
        if l[0] and l[1]:
            l = [None if not i else i for i in l]
            if l:
                result.add(Fields(*l))
            
    return list(result)