#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  This file is part of the `pypath` python module
#  MySQL client for a big Chembl database (~30GB size).
#  It was never really used. Current suggestion: to drop for pypath
#  release.
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
from past.builtins import xrange, range, reduce

import sys
from itertools import chain
import time
import threading

# from pypath:
try:
    import pypath.mysql as mysql
except:
    sys.stdout.write('\t:: No module `mysql` available.\n')
    sys.stdout.flush()

import pypath.mapping as mapping
import pypath.share.progress as progress
import pypath.share.common as common

class Chembl(object):
    def __init__(self,
                 chembl_mysql=(None, 'chembl_ebi'),
                 ncbi_tax_id=9606,
                 mapping_mysql=None,
                 mapper=None):
        self.mysql = mysql.MysqlRunner(chembl_mysql)
        self.ncbi_tax_id = ncbi_tax_id
        if mapper.__class__.__name__ != 'Mapper':
            self.mapper = mapping.Mapper(ncbi_tax_id, mapping_mysql)
            # self.mapper.load_mappings(maps=data_formats.mapListUniprot)
        else:
            self.mapper = mapper
        self.chembl_uniprot_table()
        self.result = None
        # constant elements:
        self.mandatory_fields = set(['compound_chembl',
                                     'target_uniprot',
                                     'tax_id',
                                     'target_type',
                                     'potential_duplicate'])
        self.extra_fields = [
            'compound_names', 'action_type', 'target_domains',
            'predicted_binding_domains', 'activities', 'pchembl'
        ]
        self.set_group_concat_len = '''SET group_concat_max_len=18446744073709551615;'''
        self.group_concat_len_increased = False
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
        self.cprop_select = ',\n            CAST(cp.%s AS CHAR) AS %s'
        self.cprop_join = '''
        /* various properties of the compounds */
        LEFT JOIN compound_properties AS cp
            ON md.molregno = cp.molregno'''
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
                ms.synonyms AS syn,
                md.chembl_id 
            FROM 
                molecule_synonyms AS ms 
            LEFT JOIN molecule_dictionary AS md 
                ON ms.molregno = md.molregno 
            WHERE ms.synonyms IN (%s) 
            GROUP BY 
                ms.synonyms, md.chembl_id 
            ORDER BY NULL;'''
        self.comp_rec = '''SELECT 
                cr.compound_name AS syn,
                md.chembl_id 
            FROM 
                compound_records AS cr 
            LEFT JOIN molecule_dictionary AS md 
                ON cr.molregno = md.molregno 
            WHERE cr.compound_%s %s 
            GROUP BY 
                cr.compound_name, md.chembl_id 
            ORDER BY NULL;'''

    def huge_group_concat(self):
        if not self.group_concat_len_increased:
            qid = self.mysql.get_qid(self.set_group_concat_len)
            self.mysql.send_query(self.set_group_concat_len)
            self.mysql.wait_results([qid])
            self.group_concat_len_increased = True

    def compounds_targets(self,
                          id_list,
                          id_type='uniprot',
                          assay_types=['B', 'F'],
                          relationship_types=['D', 'H'],
                          compound_props = [],
                          domains=False,
                          pred_bind_d=False,
                          action_type=False,
                          activities=False,
                          pchembl=False,
                          one_query=False,
                          client_side=False):
        '''
        Same as compounds_targets(), but queries each id by separate mysql query.
        Better performance expected in case the batch query requires disk_tmp_table.
        '''
        if id_type == 'uniprot':
            compound_lookup = True
            id_list = self.get_chembl_uniprots(id_list)
        self.result = []
        id_list = id_list if type(id_list) is list else [id_list]
        if one_query:
            query_thread = threading.Thread(
                target=self.compound_target,
                args=[id_list],
                kwargs={
                    'id_type': id_type,
                    'assay_types': assay_types,
                    'relationship_types': relationship_types,
                    'compound_props': compound_props,
                    'domains': domains,
                    'pred_bind_d': pred_bind_d,
                    'action_type': action_type,
                    'activities': activities,
                    'pchembl': pchembl
                })
            query_thread.daemon = True
            query_thread.start()
            sys.stdout.write('\n')
            sys.stdout.flush()
            while query_thread.isAlive():
                self.mysql.print_status()
                time.sleep(1)
            self.mysql_ready()
            if client_side:
                self.result = list(self.result)
        else:
            prg = progress.Progress(
                total=len(id_list), name='Starting queries', interval=5)
            qids = []
            for identifier in id_list:
                prg.step()
                qids.append(
                    self.compound_target(
                        identifier,
                        id_type=id_type,
                        assay_types=assay_types,
                        relationship_types=relationship_types,
                        compound_props=compound_props,
                        domains=domains,
                        pred_bind_d=pred_bind_d,
                        action_type=action_type,
                        activities=activities,
                        pchembl=pchembl,
                        wait=False))
            prg.terminate()
            self.mysql_ready(qids)
            for qid in qids:
                self.result.extend(list(self.mysql.get_result(qid)))

    def compound_target(self,
                        id_list,
                        id_type='uniprot',
                        assay_types=['B', 'F'],
                        relationship_types=['D', 'H'],
                        compound_props = [],
                        domains=False,
                        pred_bind_d=False,
                        action_type=False,
                        activities=False,
                        pchembl=False,
                        wait=True,
                        select_extra = '',
                        join_extra = '',
                        where_extra = ''):
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
            where = fields[id_type] + \
                ' IN (' + ','.join(['"%s"' % i for i in id_list]) + ')'
        else:
            where = fields[id_type] + ' = ' + '"%s"' % id_list
        assay_types = ','.join(['"%s"' % i for i in assay_types])
        relationship_types = ','.join(['"%s"' % i for i in relationship_types])
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
        if len(compound_props):
            for prop in compound_props:
                select_extra += self.cprop_select % (prop, prop.lower())
            join_extra += self.cprop_join
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
            ac.potential_duplicate != 1 AND 
            ac.pchembl_value IS NOT NULL AND 
                (
                data_validity_comment IS NULL OR 
                data_validity_comment = 'manually validated'
                ) AND
            ay.assay_type IN (%s) AND 
            ay.relationship_type IN (%s) AND 
            %s 
        GROUP BY cs.accession,md.chembl_id 
        ORDER BY NULL;''' % (select_extra, join_extra, assay_types,
                             relationship_types, where)
        self.huge_group_concat()
        qid = self.mysql.get_qid(q)
        cursor = self.mysql.ss_cursor if type(
            id_list) is list else self.mysql.cs_cursor
        self.mysql.send_query(q, cursor=cursor, silent=True)
        if wait:
            self.mysql_ready([qid])
        else:
            return qid
        result = self.mysql.get_result(qid)
        self.result = list(result)

    def compounds_targets_mechanism(self,
                                    id_list,
                                    id_type='uniprot',
                                    domains=False,
                                    pred_bind_d=False,
                                    activities=False,
                                    pchembl=False,
                                    one_query=False,
                                    client_side=False):
        if id_type == 'uniprot':
            compound_lookup = True
            id_list = self.get_chembl_uniprots(id_list)
        self.result = []
        id_list = id_list if type(id_list) is list else [id_list]
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
            query_thread.daemon = True
            query_thread.start()
            sys.stdout.write('\n')
            sys.stdout.flush()
            while query_thread.isAlive():
                self.mysql.print_status()
                time.sleep(1)
            self.mysql_ready()
            if client_side:
                self.result = list(self.result)
        else:
            prg = progress.Progress(
                total=len(id_list), name='Sending queries', interval=5)
            qids = []
            for identifier in id_list:
                prg.step()
                qids.append(
                    self.compound_target_mechanism(
                        identifier,
                        id_type=id_type,
                        domains=domains,
                        pred_bind_d=pred_bind_d,
                        activities=activities,
                        pchembl=pchembl,
                        wait=False))
            prg.terminate()
            self.mysql_ready(qids)
            for qid in qids:
                self.result += list(self.mysql.get_result(qid))

    def compound_target_mechanism(self,
                                  id_list,
                                  id_type='uniprot',
                                  domains=False,
                                  pred_bind_d=False,
                                  activities=False,
                                  pchembl=False,
                                  wait=True):
        fields = {
            'uniprot': 'cs.accession',
            'chembl': 'md.chembl_id',
            'compound_synonym': 'ms.synonyms'
        }
        if type(id_list) is list:
            where = fields[id_type] + \
                ' IN (' + ','.join(['"%s"' % i for i in id_list]) + ')'
        else:
            where = fields[id_type] + ' = ' + '"%s"' % id_list
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
        ''' % (select_extra, join_extra, where)
        self.huge_group_concat()
        qid = self.mysql.get_qid(q)
        self.mysql.send_query(q, silent=True)
        if wait:
            self.mysql_ready([qid])
        else:
            return qid
        result = self.mysql.get_result(qid)
        self.result = result

    def synonyms2chembl(self, synonyms, like=True):
        self.result = {}
        syn_lower = dict(zip([s.lower() for s in synonyms], synonyms))
        syn_lst = ','.join(['"%s"' % syn for syn in synonyms])
        synq = self.comp_syn % syn_lst
        recq = self.comp_rec % ('name', ' IN (%s)' % syn_lst)
        synqid = self.mysql.get_qid(synq)
        recqid = self.mysql.get_qid(recq)
        self.mysql.send_query(synq, silent=True)
        self.mysql.send_query(recq, silent=True)
        self.mysql.wait_results([synqid, recqid])
        self.mysql_ready()
        for r in chain(
                self.mysql.get_result(synqid), self.mysql.get_result(recqid)):
            syn = syn_lower[r['syn'].lower()]
            if syn not in self.result:
                self.result[syn] = []
            if r['chembl_id'] is not None:
                self.result[syn].append(r['chembl_id'])
        if like:
            like_results = {}
            notfound = [
                n for n in list(set(synonyms) - set(self.result.keys()))
                if not n.isdigit()
            ]
            qids = {}
            trds = []
            for field in ['name', 'key']:
                for syn in notfound:
                    q = self.comp_rec % (field, ' LIKE "%%%s%%"' % syn)
                    qid = self.mysql.get_qid(q)
                    qids[qid] = syn
                    self.mysql.send_query(q, silent=True)
            self.mysql.wait_results(qids.keys())
            self.mysql_ready()
            for qid, syn in iteritems(qids):
                res = self.mysql.get_result(qid)
                this_result = []
                for r in res:
                    if r['chembl_id'] is not None:
                        this_result.append(r['chembl_id'])
                if syn not in like_results:
                    like_results[syn] = []
                like_results[syn].append(this_result)
            for syn, results in iteritems(like_results):
                # choosing the shortest returned list of ChEMBL IDs
                if len(results) > 0 and syn not in self.result:
                    results = [common.uniqList(r) for r in results]
                    self.result[syn] = reduce(
                        lambda x, y: x if len(y) == 0 or len(x) < len(y) and len(x) > 0 else y,
                        results)
        self.result = dict([(k, common.uniqList(v))
                            for k, v in iteritems(self.result)])

    def get_chembl_uniprots(self, originals):
        chembls = []
        for u in originals:
            if u in self.uniprot_chembl:
                chembls += self.uniprot_chembl[u]
        return chembls

    def result_table(self):
        """
        Returns result as a tab separated string which
        can be written to a file. First row is header.
        """
        header = [self.result[0].keys()]
        s = '\t'.join(header) + '\n'
        for r in self.result:
            s += '\t'.join(r) + '\n'
        return s

    def compounds_by_target(self, update_uniprots=True):
        """
        This processes the result and saves
        the compounds to the `compounds` attribute as a dictionary.
        """
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
                for ee in sorted(r.keys()):
                    if ee not in self.mandatory_fields and \
                       ee not in self.extra_fields:
                        this_compound[ee] = r[ee]
                if u not in ct:
                    ct[u] = []
                ct[u].append(this_compound)
        self.compounds = ct

    def targets_by_compound(self):
        """
        This processes the result and saves
        the targets to the `targets` attribute as a dictionary.
        """
        tc = {}
        for r in self.result:
            this_target = {}
            this_target['uniprot'] = r['target_uniprot']
            for e in self.extra_fields:
                this_target[e] = [] if e not in r or r[e] is None else r[
                    e].split(';')
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
        qid = self.mysql.get_qid(q)
        self.mysql.send_query(q)
        self.mysql.wait_results([qid])
        chembl_uniprots = []
        for r in self.mysql.get_result(qid):
            chembl_uniprots.append(r['ac'])
        self.chembl_uniprot = {}
        self.uniprot_chembl = {}
        for u in chembl_uniprots:
            umapped = self.mapper.map_name(u, 'uniprot', 'uniprot')
            self.chembl_uniprot[u] = umapped
            for w in umapped:
                if w not in self.uniprot_chembl:
                    self.uniprot_chembl[w] = []
                self.uniprot_chembl[w].append(u)

    def mysql_ready(self, qids=None):
        if qids is not None:
            sys.stdout.write('\t:: Waiting for MySQL...')
            sys.stdout.flush()
            self.mysql.wait_results(qids)
        sys.stdout.write('\r' + ' ' * 90)
        sys.stdout.write('\r\t:: MySQL: ready.')
        sys.stdout.write('\n')
        sys.stdout.flush()


'''
some tests:

from pypath import chembl
c = chembl.Chembl((None, 'chembl_ebi'))
#c.compound_target('CHEMBL105', id_type = 'chembl')
#c.compounds_targets('CHEMBL105', id_type = 'chembl')
c.compound_target_mechanism('CHEMBL19', id_type = 'chembl')
c.compounds_targets_mechanism(['CHEMBL19'], id_type = 'chembl')

'''
