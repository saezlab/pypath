#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from future.utils import iteritems

import pypath.share.session as session_mod

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
    'SIGNOR': 'p',
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
    'ENCODE-distal': 't',
    'PAZAR': 't',
    'ENCODE-proximal': 't',
    'ORegAnno': 't',
    'HTRI': 't',
    'ARACNe-GTEx': 't',
    'DoRothEA_reviews': 't',
    'FANTOM4': 't',
    'HOCOMOCO-v11': 't',
    'HTRIdb': 't',
    'JASPAR-v2018': 't',
    'NFIRegulomeDB': 't',
    'ReMap': 't',
    'RegNetwork': 't',
    'TFactS': 't',
    'TFe': 't',
    'TRED': 't',
    'TRRD': 't',
    'TRRUST': 't',
    'DoRothEA': 't',
    # miRNA-mRNA
    'miR2Disease': 'n',
    'miRDeathDB': 'n',
    'miRecords': 'n',
    'miRTarBase': 'n',
    'ncRDeathDB': 'nw',
    # TF-miRNA
    'TransmiR': 'u',
    'ENCODE_tf-mirna': 'u',
    # lncRNA-mRNA
    'LncRNADisease': 'w',
    'lncrnadb': 'w',
}

p = set()
i = set()
r = set()
m = set()
t = set()
l = set()
n = set() # miRNA-target
u = set() # TF-mirna
w = set() # lncRNA-target

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
    'w': 'lncRNA-mRNA',
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


def get_categories(database, names = False, top = True):
    
    result = (
        {letter for letter in categories[database]}
            if database in categories else
        (
            (
                get_categories(
                    '_'.join(
                        reversed(tuple(reversed(database.split('_')))[:-1])
                    ),
                    top = False,
                )
            ) or (
                get_categories(
                    '_'.join(database.split('_')[:-1]),
                    top = False,
                )
            )
        )
            if '_' in database else
        set()
    )
    
    if not result and top:
        
        _log(
            'Could not find database `%s` in any '
            'of the categories.' % database
        )
    
    if names:
        
        result = {catnames[cat] for cat in result}
    
    return result


def get_category(database):
    
    db_categories = get_categories(database)
    
    return list(db_categories)[0] if db_categories else None
