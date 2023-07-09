from __future__ import annotations

import collections
from typing import NamedTuple

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common


class DrugcombdbSynerAntagVoting(NamedTuple):
    ID: str
    Drug1: str
    Drug2: str
    CellLine: str
    ZIP: str
    Bliss: str
    Loewe: str
    HSA: str
    ZIPclass: str
    Blissclass: str
    Loeweclass: str
    HSAclass: str
    synthetic: str
    classification: str
    
    
class DrugcombdbDrugChemicalInfo(NamedTuple):
    DrugName: str
    cIDs: str
    drugNameOfficial: str
    molecularWeight: str
    smiles: str
    
    
class DrugcombdbSyndrugcombFda(NamedTuple):
    ID: str
    Drug1: str
    Drug2: str
    Machenism: str
    Source: str
    

class DrugcombdbSyndrugcombTextmining(NamedTuple):
    Drug1: str
    Drug2: str
    Target : str
    Source : str
    

class DrugcombdbExternalSynergy(NamedTuple):
    Drug1: str
    Drug2: str
    PubmedID : str
    

class DrugcombdbExternalAntagonism(NamedTuple):
    InteractionA : str
    InteractionB : str

def drugcombdb_syner_antag_voting() -> list[tuple]:
    url = urls.urls['drugcombdb']['url_syner_antag_voting']
    c = curl.Curl(url, silent = False, large = True)
    
    result = set()
    
    for l in c.result:
        l = l.strip().split(',')
        l = [None if not i else i for i in l]
        if l and len(l) == len(DrugcombdbSynerAntagVoting._fields):
            result.add(DrugcombdbSynerAntagVoting(*l))
            
    return list(result)

def drugcombdb_drug_chemical_info() -> list[tuple]:
    url = urls.urls['drugcombdb']['url_drug_chemical_info']
    c = curl.Curl(url, silent = False, large = True, encoding = 'ISO-8859-1')
    
    result = set()
    
    for l in c.result:
        l = l.strip().split(',')
        l = [None if not i or i == '#N/A' else i for i in l]
        l = [None if not i else i for i in l]
        if l and len(l) == len(DrugcombdbDrugChemicalInfo._fields):
            result.add(DrugcombdbDrugChemicalInfo(*l))
            
    return list(result)


def drugcombdb_syndrugcomb_fda() -> list[tuple]:
    url = urls.urls['drugcombdb']['url_syndrugcomb_fda']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='2drugs')
    
    result = set()
    
    for l in contents:
        l = [None if not i else i for i in l]
        if l and len(l) == len(DrugcombdbSyndrugcombFda._fields):
            result.add(DrugcombdbSyndrugcombFda(*l))
            
    return list(result)


def drugcombdb_syndrugcomb_textmining() -> list[tuple]:
    url = urls.urls['drugcombdb']['url_syndrugcomb_textmining']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='2drugs')
    
    result = set()
    
    for l in contents:
        l = l[:-1]
        l = [None if not i else i for i in l]
        if l:
            result.add(DrugcombdbSyndrugcombTextmining(*l))
            
    return list(result)
    
    
def drugcombdb_syndrugcomb_external_synergism() -> list[tuple]:
    url = urls.urls['drugcombdb']['url_syndrugcomb_external']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='ASDCD_synergism')
    
    result = set()
    
    for l in contents:
        l = l[:3]
        l = [None if not i else i for i in l]
        if l:
            result.add(DrugcombdbExternalSynergy(*l))
            
    return list(result)


def drugcombdb_syndrugcomb_external_antagonism() -> list[tuple]:
    url = urls.urls['drugcombdb']['url_syndrugcomb_external']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, sheet='Drugbank_antagonism')
    
    result = set()
    
    for l in contents:
        if l[0] and l[1]:
            l = [None if not i else i for i in l]
            if l:
                result.add(DrugcombdbExternalAntagonism(*l))
            
    return list(result)