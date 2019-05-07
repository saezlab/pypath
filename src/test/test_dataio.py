#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#


import sys
import pytest

import pypath.dataio as dataio
import pypath.data_formats as data_formats
import pypath.settings as settings


exclude = {
    'reactome_interactions'
}
network_methods = []

for var in data_formats.__dir__():
    
    var = getattr(data_formats, var)
    
    if isinstance(var, dict):
        
        for k, v in var.items():
            
            if hasattr(v, 'inFile'):
                
                method_name = getattr(v, 'inFile')
                
                if (
                    hasattr(dataio, method_name) and
                    method_name not in exclude
                ):
                    
                    network_methods.append((method_name, v.inputArgs))


class TestDataio(object):
    
    @pytest.mark.parametrize('method, args', network_methods)
    def test_network_data_methods(self, method, args, cachedir):
        
        settings.setup(cachedir = cachedir)
        
        method = getattr(dataio, method)
        result = method(**args)
        
        assert len(result)
    
    
    def test_signor_interactions(self):
        
        s = list(dataio.signor_interactions())
        
        assert len(s) > 18000
        assert 'phosphorylation' in set(i[9] for i in s)
    
    
    def test_ptm_orthology(self):
        
        o = dataio.ptm_orthology()
        
        assert ('Q08775', 5, 'S', 125, 10090, 'phosphorylation') in o
    
    
    def test_get_oreganno(self):
        
        o = list(dataio.get_oreganno())
        
        assert ('NFAT5', 'TNF', '16027109') in o
    
    
    def test_get_pathwaycommons(self):
        
        p = dataio.get_pathwaycommons()
        
        assert ['CASK', 'interacts-with', 'PARK2', 'WikiPathways'] in p
    
    
    def test_go_annotations_goa(self):
        
        goa = dataio.go_annotations_goa()
        
        assert 'GO:0043235' in goa['C']['P00533']
    
    
    def test_go_ancestors_goose(self):
        
        g = dataio.go_ancestors_goose()
        
        assert 'GO:0003824' in g['GO:0004799']
    
    
    def test_go_terms_goose(self):
        
        g = dataio.go_terms_goose()
        
        assert g['P']['GO:0061524'] == 'central canal development'
    
    
    def test_go_annotations_goose(self):
        
        a = dataio.go_annotations_goose()
        
        assert a[0]['P']['GO:0060662'] == 'salivary gland cavitation'
    
    
    def test_get_tfcensus(self):
        
        t = dataio.get_tfcensus()
        
        assert 'FOXD4L6' in t['hgnc']
    
    
    def test_get_cellphonedb_interactions(self):
        
        example = (
            'Q8NHL6',
            'P17693',
            'CellPhoneDB',
            '',
            'ligand-receptor',
            'ligand',
            'receptor'
        )
        
        c = list(dataio.cellphonedb_interactions())
        
        assert example in c
    
    
    def test_hprd_htp(self):
        
        example = [
            'PLCG1',
            '01398',
            'NP_002651.2',
            'KIT',
            '01287',
            'NP_000213.1',
            'in vivo',
            '7536744'
        ]
        
        h = dataio.hprd_htp()
        
        assert example in h
    
    
    def test_biogrid_interactions(self):
        
        b = dataio.biogrid_interactions()
        
        assert ['SNF2', 'GCN4', '16982689'] in b
    
    
    def test_get_matrixdb(self):
        
        m = dataio.get_matrixdb()
        
        assert ['O60568', 'O60568', '20955792', 'molecular sieving'] in m
    
    
    def test_get_innatedb(self):
        
        i = dataio.get_innatedb()
        
        assert ('P27695', 'APEX1', 'P40763', 'STAT3', '15735682') in i
    
    
    def test_get_dip(self):
        
        example = [
            'P04637',
            'Q8N163',
            '18235501',
            'physical interaction',
            'anti tag coimmunoprecipitation',
            'DIP-189566E'
        ]
        
        d = dataio.get_dip()
        
        assert example in d
    
    
    def test_depod_interactions(self):
        
        example = (
            'P17706',
            'P51692',
            '11773439',
            '2',
            'dephosphorylation reaction'
        )
        
        d = dataio.depod_interactions()
        
        assert example in d
    
    
    def test_intact_interactions(self):
        
        example = [
            'Q99523',
            'Q9UJY4',
            '11331584;11390366',
            'two hybrid;pull down',
            ''
        ]
        
        i = dataio.intact_interactions()
        
        assert example in i
    
    
    def test_homologene_dict(self):
        
        h = dataio.homologene_dict(9606, 10090, 'GeneSymbol')
        
        assert 'Stard10' in h['STARD10']
    
    
    def test_mir2disease_interactions(self):
        
        example = [
            'hsa-miR-34',
            'SIRT1',
            '2009',
            'MiR-34, SIRT1 and p53: the feedback loop'
        ]
        
        m = dataio.mir2disease_interactions()
        
        assert example in m
    
    
    def test_mirdeathdb_interactions(self):
        
        m = dataio.mirdeathdb_interactions()
        
        assert ('MIMAT0004553', '1956', 9606, '22492316', 'autophagy_up') in m
    
    
    def test_lncdisease_interactions(self):
        
        example = (
            'SRA1',
            'TBL1X',
            'RNA',
            'Protein',
            'binding',
            'human',
            '20219889'
        )
        
        l = dataio.lncdisease_interactions()
        
        assert example in l
    
    
    def test_transmir_interactions(self):
        
        example = (
            'RUNX2',
            '860',
            'mir-27a',
            'Heart Failure; Autistic Disorder; Neoplasms; Stomach Neoplasms',
            'Repression',
            '20980664',
            'Human'
        )
        
        t = list(dataio.transmir_interactions())
        
        assert example in t
    
    
    def test_get_proteinatlas(self):
        
        a = dataio.get_proteinatlas()
        
        assert (
            a['normal']['adrenal gland:glandular cells']['O43657'] ==
            ('Not detected', 'Approved')
        )
    
    
    def test_get_tfregulons_old(self):
        
        example = [
            'TP53',
            'FOXN3',
            '0',
            'B',
            True,
            True,
            False,
            False,
            'PAZAR',
            'ReMap',
            '',
            '',
            'PAZAR,ReMap'
        ]
        
        t = list(dataio.get_tfregulons_old())
    
    
    def test_stitch_interactions(self):
        
        example = ('20474474', 'ENSP00000371927', 'binding', '', 168)
        
        s = dataio.stitch_interactions()
        
        s_first1000 = list(map(next(s), xrange(1000)))
        
        assert example in s_first1000
    
    
    def test_get_membranome(self):
        
        example = (
            'P53171',
            'Mitochondrial inner membrane',
            'intermembrane space'
        )
        
        m = list(dataio.get_membranome())
        
        assert example in m
    
    
    def test_get_exocarta(self):
        
        example = (
            '6201',
            'RPS7',
            9606,
            ('25890246', (9606,), 'Colorectal cancer cells'),
        )
        
        e = list(dataio.get_exocarta())
        
        assert example in e
    
    
    def test_get_vesiclepedia(self):
        
        example = (
            '2',
            'A2M',
            9606,
            (
                '29127410',
                (9606,),
                'Bronchial epithelial cells',
                ('Extracellular vesicles',)
            )
        )
        
        v = list(dataio.get_vesiclepedia())
        
        assert example in v
    
    
    def test_cellphonedb_protein_annotations(self):
        
        example = (
            None, True, False, False, False,
            None, False, True, False, False,
        )
        
        c = dataio.cellphonedb_protein_annotations()
        
        assert c['P17302'] == example
    
    
    def test_corum_complexes(self):
        
        c = dataio.corum_complexes()
        
        assert 'O95251-Q8WYH8-Q9HAF1-Q9NQC1' in c
        assert '16387653' in c['O95251-Q8WYH8-Q9HAF1-Q9NQC1'].references
    
    
    def test_havugimana_complexes(self):
        
        h = dataio.havugimana_complexes()
        
        assert 'Q8TF72-Q96ER3-Q96PV0-Q9Y6N9' in h
    
    
    def test_compleat_complexes(self):
        
        c = dataio.compleat_complexes()
        
        assert 'Q9HCE7-Q9NPB6' in c
        assert '15761148' in c['Q9HCE7-Q9NPB6'].references
    
    
    def test_complexportal_complexes(self):
        
        c = dataio.complexportal_complexes()
        
        assert 'Q9NS61-Q9NZV8' in c
    
    
    def test_cellphonedb_interactions(self):
        
        example = (
            'Q16671',
            'O00238',
            'CellPhoneDB',
            '',
            'receptor-receptor',
            'receptor',
            'receptor',
        )
        
        c = dataio.cellphonedb_interactions()
        
        assert example in c
    
    
    def test_adhesome_annotations(self):
        
        a = dataio.adhesome_annotations()
        
        assert ('Adaptor', True) in a['Q9NQ75']
    
    
    def test_adhesome_interactions(self):
        
        a = dataio.adhesome_interactions()
