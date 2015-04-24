#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import sys
import time
import threading

# from bioigraph:
import mysql
import mapping
import progress
import data_formats

class Chembl(object):
    
    def __init__(self, chembl_mysql, ncbi_tax_id = 9606, 
                 mapping_mysql = None, mapper = None):
        self.mysql = mysql.MysqlRunner(chembl_mysql)
        self.ncbi_tax_id = ncbi_tax_id
        if mapper.__class__.__name__ != 'Mapper':
            self.mapper = mapping.Mapper(ncbi_tax_id,mapping_mysql)
            # self.mapper.load_mappings(maps=data_formats.mapListUniprot)
        else:
            self.mapper = mapper
        self.chembl_uniprot_table()
        # constant elements:
        self.extra_fields = ['compound_names','action_type','target_domains',
                     'predicted_binding_domains','activities','pchembl']
        self.pbd_join = '''
        /* predicted binding domains */
        LEFT JOIN predicted_binding_domains AS pbd 
            ON ac.activity_id = pbd.activity_id 
        LEFT JOIN site_components AS sc 
            ON pbd.site_id = sc.site_id 
        LEFT JOIN domains AS pdm 
            ON sc.domain_id = pdm.domain_id '''
        self.pbd_select = ''',
        GROUP_CONCAT(
            DISTINCT(pdm.source_domain_id) 
            SEPARATOR ";") AS predicted_binding_domains'''
        self.dom_join = '''
        /* domains of the target */
        LEFT JOIN component_domains AS cd 
            ON tc.component_id = cd.component_id 
        LEFT JOIN domains AS tdm 
            ON cd.domain_id = tdm.domain_id '''
        self.dom_select = ''',
        GROUP_CONCAT(
            DISTINCT(tdm.source_domain_id) 
            SEPARATOR ";") AS target_domains'''
        self.atype_join = '''
        LEFT JOIN drug_mechanism AS dm 
            ON (dm.molregno = md.molregno AND dm.tid = td.tid)'''
        self.atype_select = ''',
            GROUP_CONCAT(DISTINCT(dm.action_type) SEPARATOR ';') AS action_type'''
        self.act_join = '''
        /* this is for activity values, if available */
        LEFT JOIN assays AS ay 
            ON dm.tid = ay.tid 
        LEFT JOIN activities AS ac 
            ON (ac.assay_id = ay.assay_id AND ac.molregno = dm.molregno)'''
        self.act_join = '''
        /* this is for activity values, if available */
        LEFT JOIN (
            SELECT 
                ay.tid,
                ac.standard_type,
                ac.standard_value,
                ac.standard_units,
                ac.pchembl_value,
                ac.molregno,
                ac.activity_id 
            FROM assays AS ay 
            INNER JOIN activities AS ac 
            ON ac.assay_id = ay.assay_id 
        ) AS ac ON (ac.molregno = dm.molregno AND ac.tid = dm.tid)'''
        self.act_select = ''',
            GROUP_CONCAT(DISTINCT(CONCAT(
                    ac.standard_type,'=',ac.standard_value,'=',ac.standard_units
                )) SEPARATOR ';') AS activities'''
        self.pchembl_select = ''',
            GROUP_CONCAT(DISTINCT(ac.pchembl_value) SEPARATOR ';') AS pchembl'''
        self.group_concat_len = '''SET group_concat_max_len=18446744073709551615;'''
        self.comp_syn = '''SELECT 
                ms.synonyms,
                md.chembl_id 
            FROM 
                molecule_synonyms AS ms 
            LEFT JOIN molecule_dictionary AS md 
                ON ms.molregno = md.molregno 
            WHERE ms.synonyms IN (%s) 
            GROUP BY 
                ms.synonyms,md.chembl_id 
            ORDER BY NULL;'''
    
    def compounds_targets(self, id_list, id_type='uniprot', assay_types=['B','F'],
                          relationship_types=['D','H'],domains=False,pred_bind_d=False,
                          action_type=False,activities=False,pchembl=False,
                          one_query=False):
        '''
        Same as compounds_targets(), but queries each id by separate mysql query.
        Better performance expected in case the batch query requires disk_tmp_table.
        '''
        if id_type == 'uniprot':
            compound_lookup = True
            id_list = self.get_chembl_uniprots(id_list)
        self.result = []
        if one_query:
            query_thread = threading.Thread(
                target=self.compound_target,
                args=[id_list],
                kwargs={
                    'id_type': id_type,
                    'assay_types':assay_types,
                    'relationship_types': relationship_types,
                    'domains': domains,
                    'pred_bind_d': pred_bind_d,
                    'action_type': action_type,
                    'activities': activities,
                    'pchembl': pchembl
                })
            query_thread.start()
            sys.stdout.write('\n')
            sys.stdout.flush()
            while query_thread.isAlive():
                self.mysql.print_status()
                time.sleep(1)
            sys.stdout.write('\r'+' '*90)
            sys.stdout.write('\r\t:: MySQL: ready.')
            sys.stdout.write('\n\n')
            sys.stdout.flush()
        else:
            prg = progress.Progress(total=len(id_list),
                                    name='Searching for compounds/targets', interval=5)
            for identifier in id_list:
                prg.step()
                self.result += self.compound_target(identifier, id_type = id_type, 
                                assay_types = assay_types, 
                                relationship_types = relationship_types,
                                domains = domains, pred_bind_d = pred_bind_d,
                                action_type = action_type, activities = activities,
                                pchembl = pchembl)
            prg.terminate()
    
    def compound_target(self, id_list, id_type = 'uniprot', assay_types = ['B','F'],
                          relationship_types = ['D','H'], 
                          domains = False, pred_bind_d = False,
                          action_type = False, activities = False, pchembl = False):
        '''Get compounds for list of targets or targets for a list of compounds
        
        Inputs:
        id_list -- list of uniprot ids or compound chembl ids
        id_type -- if 'uniprot', search compounds for targets;
            to search targets for compounds, 
        relationship_types -- list of assay relationship types,
        possible values:
        -- D - Direct protein target assigned
        -- H - Homologous protein target assigned
        -- M - Molecular target other than protein assigned
        -- N - Non-molecular target assigned
        -- S - Subcellular target assigned
        -- U - Default value - Target has yet to be curated
        assay_types -- list of assay types, possible values:
        -- A - ADME
        -- B - Binding
        -- F - Functional
        -- P - Physicochemical
        -- T - Toxicity
        -- U - Unassigned
        '''
        fields = {
            'uniprot': 'cs.accession',
            'chembl': 'md.chembl_id',
            'compound_synonym': 'ms.synonyms'
        }
        if type(id_list) is list:
            where = fields[id_type] + ' IN ('+ ','.join(['"%s"'%i for i in id_list]) + ')'
        else:
            where = fields[id_type] + ' = ' + '"%s"'%id_list
        assay_types = ','.join(['"%s"'%i for i in assay_types])
        relationship_types = ','.join(['"%s"'%i for i in relationship_types])
        select_extra = ''
        join_extra = ''
        if domains:
            select_extra += self.dom_select
            join_extra += self.dom_join
        if pred_bind_d:
            select_extra += self.pbd_select
            join_extra += self.pbd_join
        if action_type:
            select_extra += self.atype_select
            join_extra += self.atype_join
        if activities:
            select_extra += self.act_select
        if pchembl:
            select_extra += self.pchembl_select
        q = '''SELECT 
            GROUP_CONCAT(DISTINCT(ms.synonyms) SEPARATOR ';') AS compound_names,
            cs.accession AS target_uniprot,
            md.chembl_id AS compound_chembl,
            td.tax_id AS tax_id,
            tt.target_type AS target_type %s 
        FROM component_sequences AS cs 
        INNER JOIN target_components AS tc ON tc.component_id = cs.component_id 
        INNER JOIN target_dictionary AS td ON tc.tid = td.tid 
        INNER JOIN target_type AS tt ON td.target_type = tt.target_type 
        INNER JOIN assays AS ay ON td.tid = ay.tid 
        INNER JOIN activities AS ac ON ay.assay_id = ac.assay_id 
        INNER JOIN molecule_dictionary AS md ON ac.molregno = md.molregno 
        LEFT JOIN molecule_synonyms AS ms ON md.molregno = ms.molregno 
        %s 
        WHERE 
            tt.parent_type = "PROTEIN" AND 
            ac.potential_duplicate IS NULL AND 
            ac.pchembl_value IS NOT NULL AND 
                (
                data_validity_comment IS NULL OR 
                data_validity_comment = 'manually validated'
                ) AND
            ay.assay_type IN (%s) AND 
            ay.relationship_type IN (%s) AND 
            %s 
        GROUP BY cs.accession,md.chembl_id 
        ORDER BY NULL;''' % (select_extra, join_extra,
                             assay_types,relationship_types,
                             where)
        self.mysql.run_query(self.group_concat_len,silent=True)
        self.mysql.run_query(q,silent=True)
        if type(id_list) is list:
            self.result = self.mysql.result
        else:
            return self.mysql.result
    
    def compounds_targets_mechanism(self, id_list, id_type = 'uniprot',
                                    domains = False, pred_bind_d = False,
                                    activities = False, pchembl = False, 
                                    one_query = False):
        if id_type == 'uniprot':
            compound_lookup = True
            id_list = self.get_chembl_uniprots(id_list)
        self.result = []
        if one_query:
            query_thread = threading.Thread(
                target=self.compound_target_mechanism,
                args=[id_list],
                kwargs={
                    'id_type': id_type,
                    'domains': domains,
                    'pred_bind_d': pred_bind_d,
                    'activities': activities,
                    'pchembl': pchembl
                })
            query_thread.start()
            sys.stdout.write('\n')
            sys.stdout.flush()
            while query_thread.isAlive():
                self.mysql.print_status()
                time.sleep(1)
            sys.stdout.write('\r'+' '*90)
            sys.stdout.write('\r\t:: MySQL: ready.')
            sys.stdout.write('\n\n')
            sys.stdout.flush()
        else:
            prg = progress.Progress(total=len(id_list),
                                    name='Searhing for compounds/targets',interval=5)
            for identifier in id_list:
                prg.step()
                self.result += self.compound_target_mechanism(
                                identifier, id_type = id_type,
                                domains = domains, pred_bind_d = pred_bind_d,
                                activities = activities, pchembl = pchembl)
            prg.terminate()
    
    def compound_target_mechanism(self, id_list, id_type = 'uniprot',
                                    domains = False, pred_bind_d = False,
                                    activities = False, pchembl = False):
        fields = {
            'uniprot': 'cs.accession',
            'chembl': 'md.chembl_id',
            'compound_synonym': 'ms.synonyms'
        }
        if type(id_list) is list:
            where = fields[id_type] + ' IN ('+ ','.join(['"%s"'%i for i in id_list]) + ')'
        else:
            where = fields[id_type] + ' = ' + '"%s"'%id_list
        select_extra = ''
        join_extra = ''
        if domains:
            select_extra += self.dom_select
            join_extra += self.dom_join
        if activities:
            select_extra += self.act_select
        if pchembl:
            select_extra += self.pchembl_select
        if pchembl or activities or pred_bind_d:
            join_extra += self.act_join
        if pred_bind_d:
            select_extra += self.pbd_select
            join_extra += self.pbd_join
        q = '''SELECT 
            md.chembl_id AS compound_chembl, 
            GROUP_CONCAT(DISTINCT(dm.action_type) SEPARATOR ";") AS action_type,
            GROUP_CONCAT(DISTINCT(ms.synonyms) SEPARATOR ";") AS compound_names,
            cs.accession AS target_uniprot %s 
        FROM drug_mechanism AS dm 
        /* target side */
        INNER JOIN target_dictionary AS td 
            ON dm.tid = td.tid 
        INNER JOIN target_components AS tc 
            ON td.tid = tc.tid 
        INNER JOIN component_sequences AS cs 
            ON tc.component_id = cs.component_id 
        /* target optionals */
        LEFT JOIN target_type AS tt 
            ON td.target_type = tt.target_type 
        /* compound side */
        INNER JOIN molecule_dictionary AS md 
            ON dm.molregno = md.molregno 
        LEFT JOIN molecule_synonyms AS ms 
            ON md.molregno = ms.molregno 
        %s 
        WHERE 
            tt.parent_type = "PROTEIN" AND 
            %s 
        GROUP BY 
            compound_chembl,
            target_uniprot
        ORDER BY NULL;
        ''' % (select_extra,join_extra,where)
        self.mysql.run_query(self.group_concat_len,silent=True)
        self.mysql.run_query(q,silent=True)
        if type(id_list) is list:
            self.result = self.mysql.result
        else:
            return self.mysql.result
    
    def synonyms2chembl(self,synonyms):
        self.result = {}
        syn_lst = ','.join(['"'+syn+'"' for syn in synonyms])
        query_thread = threading.Thread(
            target=self.mysql.run_query,
            kwargs={
                'query': self.comp_syn % syn_lst,
                'silent': True
            })
        query_thread.start()
        sys.stdout.write('\n')
        sys.stdout.flush()
        while query_thread.isAlive():
            self.mysql.print_status()
            time.sleep(1)
        sys.stdout.write('\r'+' '*90)
        sys.stdout.write('\r\t:: MySQL: ready.')
        sys.stdout.write('\n\n')
        sys.stdout.flush()
        for r in self.mysql.result:
            if r['synonyms'] not in self.result:
                self.result[r['synonyms']] = []
            if r['chembl_id'] is not None:
                self.result[r['synonyms']].append(r['chembl_id'])
    
    def get_chembl_uniprots(self,originals):
        chembls = []
        for u in originals:
                if u in self.uniprot_chembl:
                    chembls += self.uniprot_chembl[u]
        return chembls
    
    def result_table(self):
        header = [self.result[0].keys()]
        s = '\t'.join(header) + '\n'
        for r in self.result:
            s += '\t'.join(r) + '\n'
        return s
    
    def compounds_by_target(self, update_uniprots = True):
        ct = {}
        for r in self.result:
            if r['target_uniprot'] in self.chembl_uniprot and update_uniprots:
                uniprots = self.chembl_uniprot[r['target_uniprot']]
            else:
                uniprots = [r['target_uniprot']]
            for u in uniprots:
                this_compound = {}
                this_compound['chembl'] = r['compound_chembl']
                for e in self.extra_fields:
                    this_compound[e] = [] if e not in r or r[e] is None \
                        else r[e].split(';')
                if u not in ct:
                    ct[u] = []
                ct[u].append(this_compound)
        self.compounds = ct
    
    def targets_by_compound(self):
        tc = {}
        for r in self.result:
            this_target = {}
            this_target['uniprot'] = r['target_uniprot']
            for e in self.extra_fields:
                this_target[e] = [] if e not in r or r[e] is None else r[e].split(';')
            if r['compound_chembl'] not in tc:
                tc[r['compound_chembl']] = []
            tc[r['compound_chembl']].append(this_target)
        self.targets = tc
    
    def chembl_uniprot_table(self):
        q = '''SELECT DISTINCT(accession) AS ac 
        FROM component_sequences 
        WHERE 
        component_type = 'PROTEIN' AND 
        tax_id = %u ORDER BY NULL;''' % (self.ncbi_tax_id)
        self.mysql.run_query(q)
        if self.mysql.result is not None:
            chembl_uniprots = []
            for r in self.mysql.result:
                chembl_uniprots.append(r['ac'])
            self.chembl_uniprot = {}
            self.uniprot_chembl = {}
            for u in chembl_uniprots:
                umapped = self.mapper.map_name(u,'uniprot','uniprot')
                self.chembl_uniprot[u] = umapped
                for w in umapped:
                    if w not in self.uniprot_chembl:
                        self.uniprot_chembl[w] = []
                    self.uniprot_chembl[w].append(u)
    