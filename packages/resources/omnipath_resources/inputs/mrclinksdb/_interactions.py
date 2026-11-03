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

from collections.abc import Generator
import urllib

from pypath.share import curl
from pypath.inputs import uniprot
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
import pypath.internals.intera as intera
from . import _records

__all__ = [
    'mrclinksdb_receptor_raw',
    'mrclinksdb_receptor_interaction',
    'mrclinksdb_transporter_raw',
    'mrclinksdb_transporter_interaction',
    'mrclinksdb_metabolite_cell',
    # Backward-compatible aliases
    'mrclinksdb_raw',
    'mrclinksdb_interaction',
]


def mrclinksdb_receptor_raw(
        organism: int | str = 'human',
    ) -> Generator[_records.MrclinksdbReceptorRaw, None, None]:
    """
    Raw rows from the MRCLinksDB metabolite ligand-receptor interaction file.

    Args:
        organism:
            Name or NCBI Taxonomy ID of an organism.

    Yields:
        :class:`~_records.MrclinksdbReceptorRaw` records, including the
        header row as the first element.
    """

    latin_name = taxonomy.ensure_latin_name(organism)
    url = urls.urls['mrclinksdb']['url'] % urllib.parse.quote(latin_name)

    c = curl.Curl(url, large = True, silent = False)

    rows = iter(c.result)

    # Read header to detect column order.
    # In some organisms (e.g. mouse) the 'Receptor gene ID' and
    # 'Receptor uniprot id' columns are swapped relative to the
    # MrclinksdbReceptorRaw namedtuple definition (gene ID at pos 9,
    # UniProt at pos 10).
    header_line = next(rows, None)

    if header_line is None:
        return

    header = header_line.split('\t')
    swap_gene_uniprot = (
        len(header) > 10 and
        header[9].strip().lower() == 'receptor uniprot id'
    )

    yield _records.MrclinksdbReceptorRaw(*header)

    for line in rows:

        fields = line.split('\t')

        if swap_gene_uniprot and len(fields) > 10:
            fields[9], fields[10] = fields[10], fields[9]

        yield _records.MrclinksdbReceptorRaw(*fields)


def mrclinksdb_receptor_interaction(
        organism: int | str = 'human',
    ) -> Generator[_records.MrclinksdbReceptorInteraction, None, None]:
    """
    Ligand-receptor interactions from MRCLinksDB.

    Args:
        organism:
            Name or NCBI Taxonomy ID of an organism.

    Yields:
        :class:`~_records.MrclinksdbReceptorInteraction` records.
    """

    lines = list(mrclinksdb_receptor_raw(organism))
    uniprot_id = [row.receptor_uniprot_id for row in lines[1:]]
    all_uniprots = sorted({
        u
        for uu in uniprot_id
        for u in uu.split('_')
    })
    uniprot_locs = {}
    chunk_size = 98

    for i in range(0, len(all_uniprots), chunk_size):

        uniprot_locs.update(
            uniprot.uniprot_locations(
                accession = all_uniprots[i: i + chunk_size],
                organism = None,
                reviewed = None,
            )
        )

    for rec in lines:

        if '_' in (uid := rec.receptor_uniprot_id):

            uniprots = uid.split('_')

            all_locs = [
                loc
                for u in uniprots
                if (loc := uniprot_locs.get(u))
            ]

            if all_locs:
                try:
                    location = sorted([item for item in set.union(*all_locs) if item is not None])
                except TypeError:
                    location = [item for item in set.union(*all_locs) if item is not None]
            else:
                location = []

            uid = intera.Complex(
                components = sorted(uniprots),
                ncbi_tax_id = 9606,
                sources = 'MRClinksDB',
            )

        else:
            loc = uniprot_locs.get(uid)

            if loc is None or not loc:
                location = []
            else:
                filtered = [item for item in loc if item is not None]
                try:
                    location = sorted(filtered)
                except TypeError:
                    location = filtered

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

        yield _records.MrclinksdbReceptorInteraction(
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
            receptor_uniprot = uid,
            receptor_genesymbol = rec.receptor_symbol,
            pmids = pubmeds,
            resource = rec.other_db_,
            receptor_location = location,
        )


def mrclinksdb_transporter_raw(
        organism: int | str = 'human',
    ) -> Generator[_records.MrclinksdbTransporterRaw, None, None]:
    """
    Raw rows from the MRCLinksDB transporter protein file.

    The human file has 8 columns including a HMDB protein ID (``HAMDBP_ID``);
    the mouse file has 7 columns with that field absent.  Both are normalised
    to the 8-field :class:`~_records.MrclinksdbTransporterRaw` schema with
    ``hmdbp_id`` set to ``''`` for mouse rows.

    Args:
        organism:
            Name or NCBI Taxonomy ID of an organism.

    Yields:
        :class:`~_records.MrclinksdbTransporterRaw` records, including the
        header row as the first element.
    """

    latin_name = taxonomy.ensure_latin_name(organism)
    url = urls.urls['mrclinksdb']['transporter'] % urllib.parse.quote(latin_name)

    c = curl.Curl(url, large = True, silent = False)

    rows = iter(c.result)

    header_line = next(rows, None)

    if header_line is None:
        return

    header = [h.strip() for h in header_line.split('\t')]

    # Human has HAMDBP_ID at position 2; mouse omits it entirely.
    # The column is spelled 'HAMDBP_ID' (with 'A'), so detect by 'hamdbp'.
    has_hmdbp = any('hamdbp' in h.lower() for h in header)

    # Yield a normalised header row (8 fields).
    if has_hmdbp:
        yield _records.MrclinksdbTransporterRaw(*header)
    else:
        yield _records.MrclinksdbTransporterRaw(
            header[0],
            header[1],
            'hmdbp_id',     # synthetic column name for the header row
            *header[2:],
        )

    for line in rows:

        fields = [f.strip() for f in line.split('\t')]

        if not has_hmdbp:
            # Insert empty hmdbp_id at position 2 to match the 8-field schema.
            fields = fields[:2] + [''] + fields[2:]

        yield _records.MrclinksdbTransporterRaw(*fields)


def mrclinksdb_transporter_interaction(
        organism: int | str = 'human',
    ) -> Generator[_records.MrclinksdbTransporterInteraction, None, None]:
    """
    Transporter-metabolite interactions from MRCLinksDB.

    Args:
        organism:
            Name or NCBI Taxonomy ID of an organism.

    Yields:
        :class:`~_records.MrclinksdbTransporterInteraction` records with
        UniProt subcellular location data.
    """

    lines = list(mrclinksdb_transporter_raw(organism))

    # Skip the header row (first element).
    data_rows = lines[1:]

    all_uniprots = sorted({
        row.uniprot_id
        for row in data_rows
        if row.uniprot_id
    })

    uniprot_locs = {}
    chunk_size = 98

    for i in range(0, len(all_uniprots), chunk_size):

        uniprot_locs.update(
            uniprot.uniprot_locations(
                accession = all_uniprots[i: i + chunk_size],
                organism = None,
                reviewed = None,
            )
        )

    for rec in data_rows:

        if not rec.uniprot_id or not rec.hmdb_id:
            continue

        loc = uniprot_locs.get(rec.uniprot_id)

        if loc is None or not loc:
            location = []
        else:
            filtered = [item for item in loc if item is not None]
            try:
                location = sorted(filtered)
            except TypeError:
                location = filtered

        yield _records.MrclinksdbTransporterInteraction(
            hmdb = rec.hmdb_id,
            metabolite_name = rec.metabolite_name,
            hmdbp_id = rec.hmdbp_id,
            enzyme_name = rec.enzyme_name,
            transporter_entrez = rec.gene_id,
            transporter_genesymbol = rec.gene_name,
            transporter_uniprot = rec.uniprot_id,
            transporter_location = location,
        )


def mrclinksdb_metabolite_cell() -> Generator[_records.MrclinksdbMetaboliteCell, None, None]:

    url = urls.urls['mrclinksdb']['metabolite_cell']
    c = curl.Curl(url, large = True, silent = False)

    for line in c.result:

        line = line.split('\t')

        yield _records.MrclinksdbMetaboliteCell(*line)


# Backward-compatible aliases — callers using the old names continue to work.
mrclinksdb_raw = mrclinksdb_receptor_raw
mrclinksdb_interaction = mrclinksdb_receptor_interaction
