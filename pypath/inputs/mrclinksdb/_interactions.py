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

from collections import namedtuple
from collections.abc import Generator
import collections
import urllib
import re

from pypath.share import curl
from pypath.inputs import uniprot
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
import pypath.internals.intera as intera
from . import _records

__all__ = [
    'mrclinksdb_raw',
    'mrclinksdb_interaction',
    'mrclinksdb_metabolite_cell',
]


def mrclinksdb_raw(
        organism: int | str = 'human',
    ) -> Generator[list, None, None]:

    latin_name = taxonomy.ensure_latin_name(organism)
    url = urls.urls["mrclinksdb"]["url"] % urllib.parse.quote(latin_name)

    c = curl.Curl(url, large = True, silent = False)

    for line in c.result:

        yield _records.MrclinksdbRaw(*line.split("\t"))


def mrclinksdb_interaction(
        organism: int | str = 'human',
    ) -> Generator[list, None, None]:
    """
    Ligand-receptor interactions from MRClinksDB.

    Args:
        organism:
            Name or NCBI Taxonomy ID of an organism.
    """

    lines = list(mrclinksdb_raw(organism))
    uniprot_id = [row.receptor_uniprot_id for row in lines[1:]]
    all_uniprots = sorted({
        u
        for uu in uniprot_id
        for u in uu.split('_')
    })
    uniprot_locations = {}
    chunk_size = 98

    for i in range(0, len(all_uniprots), chunk_size):

        uniprot_locations.update(
            uniprot.uniprot_locations(
                accession=all_uniprots[i: i + chunk_size],
                organism=None,
                reviewed=None,
            )
        )

    for rec in lines:

        if '_' in (uniprot_id := rec.receptor_uniprot_id):

            uniprots = uniprot_id.split('_')
            location = sorted(set.union(*(
                loc
                for u in uniprots
                if (loc := uniprot_locations.get(u))
            )))


            uniprot_id = intera.Complex(
                components = sorted(uniprots),
                ncbi_tax_id = 9606,
                sources = 'MRClinksDB',
            )

        else:

            location = uniprot_locations.get(uniprot_id)

        pubchem, pubchem_sid = rec.pubchem_cid_sid.split(';')
        pubmeds = (
            sorted(
                p.strip()
                for p in rec.pmid.split(';')
            )
                if rec.pmid else
            []
        )

        yield _records.MrclinksdbInteraction(
            mrid = rec.mrid,
            hmdb = rec.hmdb_id,
            name = rec.metabolite_name,
            pubchem = pubchem,
            pubchem_sid = pubchem_sid,
            formula = rec.molecular_formula,
            compound_kingdom = rec.kingdom,
            compound_superclass = rec.super_class_,
            compound_class = rec.class_,
            smiles = rec.canonical_smiles,
            receptor_entrez = rec.receptor_gene_id,
            receptor_uniprot = uniprot_id,
            receptor_genesymbol = rec.receptor_symbol,
            pmids = pubmeds,
            resource = rec.other_db_,
            receptor_location = location,
        )


def metabolite_cell() -> Generator[Metabolite_cell, None, None]:
    url = "https://www.cellknowledge.com.cn/mrclinkdb/download/Metabolite-cell%20interaction.txt"
    c = curl.Curl(url)
    metabolite_cell = c.result
    lines = metabolite_cell.strip('\n').split('\n')
    for line in lines:
        hmdb_id, metabolite, interaction, cell_type, experimental_subject, disease, effect, _, pubmed_id = line.split(
            '\t')
        yield Metabolite_cell(hmdb_id, metabolite, interaction, cell_type, experimental_subject, disease, effect,
                              pubmed_id)






