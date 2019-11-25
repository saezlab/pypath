#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Creates a complex meta-database accessible through either pypath
#  or Omnipath.
#
#  Copyright
#  2014-2019
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

from future.utils import iteritems

import pypath.session_mod as session_mod

_logger = session_mod.Logger(name = 'db_categories')
_log = _logger._log


categories = {
    'Vidal HI-III': 'i',
    'CancerCellMap': 'p',
    'InnateDB': 'i',
    'SPIKE': 'p',
    'LMPID': 'm',
    'DIP': 'i',
    'HPRD': 'i',
    'HPRD-phos': 'm',
    'PDZBase': 'p',
    'dbPTM': 'm',
    'MatrixDB': 'i',
    'DOMINO': 'm',
    'Signor': 'p',
    'Macrophage': 'p',
    'Adhesome': 'p',
    'NetPath': 'r',
    'ELM': 'm',
    'SignaLink2': 'p',
    'SignaLink3': 'p',
    'NRF2ome': 'p',
    'DEPOD': 'm',
    'phosphoELM': 'm',
    'MPPI': 'i',
    'Guide2Pharma': 'l',
    'Guide2Pharma_CP': 'l',
    'TRIP': 'p',
    'AlzPathway': 'r',
    'PhosphoSite': 'm',
    'CA1': 'p',
    'NCI-PID': 'r',
    'DeathDomain': 'p',
    'ARN': 'p',
    'BioGRID': 'i',
    'IntAct': 'i',
    'Reactome': 'r',
    'ACSN': 'r',
    'WikiPathways': 'r',
    'TRIP': 'p',
    'PANTHER': 'r',
    'ABS': 't',
    'MIMP': 'm',
    'PhosphoNetworks': 'm',
    'Li2012': 'm',
    'PhosphoPoint': 'm',
    'PhosphoSite_noref': 'm',
    'Ramilowski2015': 'l',
    'Kirouac2010': 'l',
    'HPMR': 'l',
    'CellPhoneDB': 'l',
    'Guide2Pharma': 'l',
    'GO_lig_rec': 'l',
    'guidetopharmacology.org': 'l',
    'UniProt': 'l',
    'InnateDB-All': 'i',
    'MINT': 'i',
    'HIPPIE': 'i',
    'Wang': 'p',
    'KEGG': 'p',
    'ProtMapper': 'm',
    'RLIMS-P': 'm',
    'REACH': 'm',
    # TF-target
    'ENCODE_distal': 't',
    'PAZAR': 't',
    'ENCODE_proximal': 't',
    'ORegAnno': 't',
    'HTRI': 't',
    'ARACNe-GTEx': 't',
    'DoRothEA_reviews': 't',
    'FANTOM4': 't',
    'HOCOMOCO_v11': 't',
    'HTRIdb': 't',
    'JASPAR_v2018': 't',
    'NFIRegulomeDB': 't',
    'ReMap': 't',
    'RegNetwork': 't',
    'TFactS': 't',
    'TFe': 't',
    'TRED': 't',
    'TRRD': 't',
    'TRRUST': 't',
    # miRNA-mRNA
    'miR2Disease': 'n',
    'miRDeathDB': 'n',
    'miRecords': 'n',
    'miRTarBase': 'n',
    # TF-miRNA
    'TransmiR': 'u',
    'ENCODE_tf-mirna': 'u',
}

p = set()
i = set()
r = set()
m = set()
t = set()
l = set()
n = set() # miRNA-target
u = set() # TF-mirna

for db, cats in iteritems(categories):
    
    for c in cats:
        
        locals()[c].add(db)

catnames = {
    'm': 'Enzyme-substrate',
    'p': 'Activity flow',
    'i': 'Undirected PPI',
    'r': 'Process description',
    't': 'Transcription',
    'l': 'Ligand-receptor',
    'n': 'miRNA-mRNA',
    'u': 'TF-miRNA',
    '':  'No category',
    None: 'No category',
}

catletters = dict(map(reversed, iteritems(catnames)))

pathway_resources = p
interaction_resources = i
ptm_resources = m
reaction_resources = r
transctiption_resources = t
ligand_receptor_resources = l


def get_categories(database):
    
    result = (
        {letter for letter in categories[database]}
            if database in categories else
        get_categories(
            '_'.join(
                reversed(tuple(reversed(database.split('_')))[:-1])
            )
        )
            if '_' in database else
        set()
    )
    
    if not result:
        
        _log(
            'Could not find database `%s` in any '
            'of the categories.' % database
        )
    
    return result


def get_category(database):
    
    db_categories = get_categories(database)
    
    return list(db_categories)[0] if db_categories else None
