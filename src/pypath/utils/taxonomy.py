#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Contains helper functions shared by different modules.
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import pypath.share.common as common
import pypath.share.session as session_mod

_logger = session_mod.Logger(name = 'taxonomy')
_log = _logger._log

# XXX: Shouldn't we keep all functions and variables separated
#      (together among them)?
taxids = {
    9606: 'human',
    10090: 'mouse',
    10116: 'rat',
    9031: 'chicken',
    9913: 'cow',
    9986: 'rabbit',
    9940: 'sheep',
    10141: 'guinea pig',
    10036: 'hamster',
    7227: 'fruit fly',
    9615: 'dog',
    9823: 'pig',
    8355: 'frog',
    9091: 'quail',
    9796: 'horse',
    9925: 'goat',
    89462: 'water buffalo',
    9598: 'monkey',
    9103: 'turkey',
    9685: 'cat',
    7604: 'starfish',
    7609: 'spiny starfish',
    1213717: 'torpedo',
    9669: 'ferret',
    8839: 'duck',
    9593: 'gorilla',
    7460: 'honeybee',
    8407: 'european common frog',
    9544: 'rhesus macaque',
}


taxa = common.swap_dict(taxids)


taxa_synonyms = {
    'bovine': 'cow',
}


phosphoelm_taxids = {
    9606: 'Homo sapiens',
    10090: 'Mus musculus',
    9913: 'Bos taurus',
    9986: 'Oryctolagus cuniculus',
    9615: 'Canis familiaris',
    10029: 'Cricetulus griseus',
    9267: 'Didelphis virginiana',
    9031: 'Gallus gallus',
    10036: 'Mesocricetus auratus',
    9940: 'Ovis aries',
    10116: 'Rattus norvegicus',
    9823: 'Sus scrofa',
    8355: 'Xenopus laevis',
    10141: 'Cavia porcellus',
    9796: 'Equus caballus',
    7227: 'Drosophila melanogaster',
    487: 'Neisseria meningitidis',
    562: 'Escherichia coli',
    5207: 'Cryptococcus neoformans',
    470: 'Acinetobacter baumannii',
    1280: 'Staphylococcus aureus',
    4939: 'Saccharomyces cerevisiae',
    34: 'Myxococcus xanthus',
    1392: 'Bacillus anthracis',
    210: 'Helicobacter pylori',
    6239: 'Caenorhabditis elegans',
}

dbptm_taxids = {
    9606: 'HUMAN',
    10090: 'MOUSE',
    7227: 'DROME',
    10116: 'RAT',
    559292: 'YEAST',
    284812: 'SCHPO',
    4081: 'SOLLC',
    3702: 'ARATH',
    9940: 'SHEEP',
    9913: 'BOVIN',
    9925: 'CAPHI',
    44689: 'DICDI',
    4577: 'MAIZE',
    9823: 'PIG',
    9615: 'CANLF',
    6239: 'CAEEL',
    8455: 'XENLA',
    83333: 'ECOLI',
    1891767: 'SV40',
}


mirbase_taxids = {
    9606: 'hsa',
    10090: 'mmu',
    10116: 'rno',
    7227: 'dme',
}


nonstandard_taxids = {
    'drosophila': 7227,
    'c.elegans': 6239,
    'xenopus': 8355,
    'Synechocystis_sp.': 1142,
}


dbptm_to_ncbi_tax_id = common.swap_dict(dbptm_taxids)
latin_name_to_ncbi_tax_id = common.swap_dict(phosphoelm_taxids)


def taxid_from_common_name(taxon_name):
    
    taxon_name = taxon_name.lower().strip()
    
    if not taxon_name or taxon_name in {'none', 'unknown'}:
        
        return None
    
    if taxon_name in taxa_synonyms:
        
        taxon_name = taxa_synonyms[taxon_name]
    
    if taxon_name in taxa:
        
        return taxa[taxon_name]


def taxid_from_latin_name(taxon_name):
    
    if taxon_name in latin_name_to_ncbi_tax_id:
        
        return latin_name_to_ncbi_tax_id[taxon_name]


def taxid_from_dbptm_taxon_name(taxon_name):
    
    if taxon_name in dbptm_to_ncbi_tax_id:
        
        return dbptm_to_ncbi_tax_id[taxon_name]


def taxid_from_nonstandard(taxon_name):
    
    if taxon_name in nonstandard_taxids:
        
        return nonstandard_taxids[taxon_name]


def ensure_ncbi_tax_id(taxon_id):
    
    if isinstance(taxon_id, int):
        
        return taxon_id
        
    else:
        
        if '(' in taxon_id:
            
            part0, part1 = taxon_id.split('(', maxsplit = 1)
            
            ncbi_tax_id = (
                ensure_ncbi_tax_id(part0) or
                ensure_ncbi_tax_id(part1.split(')', maxsplit = 1)[0])
            )
            
        elif hasattr(taxon_id, 'isdigit') and taxon_id.isdigit():
            
            nxbi_tax_id = int(taxon_id)
            
        else:
            
            ncbi_tax_id = (
                taxid_from_common_name(taxon_id) or
                taxid_from_latin_name(taxon_id) or
                taxid_from_dbptm_taxon_name(taxon_id) or
                taxid_from_nonstandard(taxon_id)
            )
        
        if not ncbi_tax_id:
            
            _log('Could not map to NCBI Taxonomy ID: `%s`.' % str(taxon_id))
        
        return ncbi_tax_id
