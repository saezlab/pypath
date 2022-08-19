#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Tennur Kılıç
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from typing import List, Dict

import csv
import collections

import pypath.utils.mapping as map
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.formats.obo as obo

def hpo_gene_annotations() -> Dict[str, list]:
    """
    Retrieves Gene-HPO relationships from HPO.

    Returns:
        namedtuple.
    """

    url = urls.urls['hpo']['gene']
    c = curl.Curl(url, large = True, silent = False)

    gene = list(csv.DictReader(c.result, delimiter = ','))

    fields = ('entrez_gene_id','entrez_gene_symbol','HPO_Term_ID')

    HPOGeneAnnotations = collections.namedtuple('HPOGeneAnnotations', fields,defaults = ("",) * len(fields))

    annotations = collections.defaultdict(list)

    for rec in gene:

        values = rec.values()
        values = list(values)[0].replace('\t',',').split(',')
        id = map.map_name(values[1], 'genesymbol', 'uniprot')
        id = list(id)

        if id:

            annotations[id[0]].append(
                HPOGeneAnnotations(
                    entrez_gene_id = values[0],
                    entrez_gene_symbol = values[1],
                    HPO_Term_ID = values[2],
                    )
            )

    return annotations

def hpo_disease_annotations() -> List[tuple] :
    """
    Retrieves Disease-HPO relationships from HPO.

    Returns:
        namedtuple.
    """

    url = urls.urls['hpo']['disease']
    c = curl.Curl(url, large = True, silent = False)

    disease = list(csv.DictReader(c.result, delimiter = '\t'))

    fields = ('DatabaseID', 'DiseaseName', 'Qualifier', 'HPO_ID', 'Reference', 'Evidence', 'Aspect')

    HPODiseaseAnnotations = collections.namedtuple('HPODiseaseAnnotations', fields,defaults = ("",) * len(fields))

    result = []

    for i in range(4,len(disease)):

        values = disease[i].values()
        values = list(values)

        result.append(
            HPODiseaseAnnotations(
                DatabaseID = values[0],
                DiseaseName = values[1][0],
                Qualifier = values[1][1],
                HPO_ID = values[1][2],
                Reference = values[1][3],
                Evidence = values[1][4],
                Aspect = values[1][9],
                )
            )


    return result

def hpo_ontology() -> List[tuple] :
    """
    Retrieves ontology from HPO.

    Returns:
        namedtuple.
    """

    url = urls.urls['hpo']['ontology']
    reader = obo.Obo(url)
    hpo_ontology = [i for i in reader]


    fields = ('hpo_id','term_name','synonyms','xrefs','is_a')

    Ontology = collections.namedtuple('Ontology', fields,defaults = ("",) * len(fields))


    result = []

    for rec in hpo_ontology:

        syn_lst = []
        xref_lst = []
        isa_lst = []

        if rec[2][1]:

            name = rec[2][0] + " " + rec[2][1]

        else:

            name = rec[2][0]

        result.append(
            Ontology(
                hpo_id = rec[1][0],
                term_name = name,
            )
        )

        if rec[5].get('synonym'):

            synonym = list(rec[5].get('synonym'))

            for i in synonym:

                syn = i[0] + " " + i[1]
                syn_lst.append(syn)

            result[-1] = result[-1]._replace(
                synonyms = syn_lst
            )

        if rec[5].get('xref'):

            xref = list(rec[5].get('xref'))

            for i in xref:

                xref_lst.append(i[0])

            result[-1] = result[-1]._replace(
                xrefs = xref_lst
            )

        if rec[5].get('is_a'):

            is_a = list(rec[5].get('is_a'))

            for i in is_a:

                isa_lst.append(i[0] + " : " + i[2])

            result[-1] = result[-1]._replace(
                is_a = isa_lst
            )

    return result