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
    url = urls.urls['mrclinksdb']['url'] % urllib.parse.quote(latin_name)

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

            all_locs = [
                loc
                for u in uniprots
                if (loc := uniprot_locations.get(u)) is not None and loc
            ]

            if all_locs:
                combined = set.union(*all_locs)
                filtered = [item for item in combined if item is not None]
                try:
                    location = sorted(filtered)
                except TypeError:
                    location = filtered
            else:
                location = []

            uniprot_id = intera.Complex(
                components = sorted(uniprots),
                ncbi_tax_id = 9606,
                sources = 'MRClinksDB',
            )

        else:
            loc = uniprot_locations.get(uniprot_id)
            if loc is None or not loc:
                location = []
            else:
                filtered = [item for item in loc if item is not None]
                try:
                    location = sorted(filtered)
                except TypeError:
                    location = filtered

        #pubchem, pubchem_sid = rec.pubchem_cid_sid.split(';')
        parts = rec.pubchem_cid_sid.split(';')
        pubchem = parts[0] if parts else ''
        pubchem_sid = parts[1] if len(parts) > 1 else ''

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


def mrclinksdb_metabolite_cell() -> Generator[_records.MrclinksdbMetaboliteCell, None, None]:

    url = urls.urls['mrclinksdb']['metabolite_cell']
    c = curl.Curl(url, large = True, silent = False)

    for line in c.result:

        line = line.split('\t')

        yield _records.MrclinksdbMetaboliteCell(*line)
