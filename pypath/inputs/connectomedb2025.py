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

"""
ConnectomeDB2025: human ligand-receptor interactions.

https://connectomedb.org/

Separate from the legacy `connectomedb` (2020/NATMI) module. The human release
is a single CSV; PubMed references are embedded in the per-interaction
``AI summary`` URL. Both `Direct` and `Inferred` evidence rows are loaded, with
the evidence type kept as an edge attribute.
"""

from __future__ import annotations

import collections
import csv
import re
import itertools

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session

_log = session.Logger(name = 'connectomedb2025_input')._log

# PubMed IDs are embedded in the `AI summary` URL, e.g. `...Pubmed-ID-123, 456`
_repmid = re.compile(r'Pubmed-ID-([\d,\s]+)')

# ConnectomeDB2025 location wording -> intercell annotation vocabulary
_LOCATIONS = {
    'secreted': 'secreted',
    'cell membrane': 'plasma membrane',
}


Connectomedb2025Interaction = collections.namedtuple(
    'Connectomedb2025Interaction',
    (
        'ligand',             # 0: ligand gene symbol
        'receptor',           # 1: receptor gene symbol
        'references',         # 2: ';'-joined PubMed IDs
        'evidence',           # 3: `Direct` or `Inferred`
        'interaction_id',     # 4: ConnectomeDB interaction id (e.g. CDB15:...)
        'ligand_location',    # 5: raw ligand location
        'receptor_location',  # 6: raw receptor location
    ),
)


def _symbol(value: str) -> str:
    """Primary gene symbol (drops parenthesised aliases)."""

    return (value or '').split('(')[0].strip()


def _pmids(ai_summary: str) -> list[str]:

    m = _repmid.search(ai_summary or '')

    if not m:

        return []

    return [p.strip() for p in m.group(1).split(',') if p.strip().isdigit()]


def connectomedb2025_raw() -> list[Connectomedb2025Interaction]:
    """
    Parses the ConnectomeDB2025 human ligand-receptor CSV.
    """

    c = curl.Curl(
        urls.urls['connectomedb2025']['url'],
        silent = False,
        large = True,
    )

    result = []

    for row in csv.DictReader(c.result):

        ligand = _symbol(row['Ligand Symbols'])
        receptor = _symbol(row['Receptor Symbols'])

        if not ligand or not receptor:

            continue

        result.append(
            Connectomedb2025Interaction(
                ligand = ligand,
                receptor = receptor,
                references = ';'.join(_pmids(row['AI summary'])),
                evidence = row['Evidence'] or None,
                interaction_id = row['Interaction ID'] or None,
                ligand_location = row['Ligand Location'] or None,
                receptor_location = row['Receptor Location'] or None,
            )
        )

    _log(f'ConnectomeDB2025: {len(result)} interactions loaded.')

    return result


def connectomedb2025_interactions() -> list[Connectomedb2025Interaction]:
    """
    Human ligand-receptor interactions from ConnectomeDB2025, with gene symbol
    identifiers (mapped to UniProt by the network reader).
    """

    return connectomedb2025_raw()


def _norm_locations(raw: str) -> set:

    locations = set()

    for part in (raw or '').split(','):

        norm = _LOCATIONS.get(part.strip().lower())

        if norm:

            locations.add(norm)

    return locations


def connectomedb2025_annotations() -> dict:
    """
    Intercell role/location annotations from ConnectomeDB2025.

    Returns a dict mapping UniProt IDs to sets of `Connectomedb2025Annotation`
    named tuples carrying the `role` (`ligand`/`receptor`) and the normalised
    cellular `location`.
    """

    import pypath.utils.mapping as mapping

    Connectomedb2025Annotation = collections.namedtuple(
        'Connectomedb2025Annotation',
        ('role', 'location'),
    )

    result = collections.defaultdict(set)

    for rec in connectomedb2025_raw():

        for symbol, role, loc_raw in (
            (rec.ligand, 'ligand', rec.ligand_location),
            (rec.receptor, 'receptor', rec.receptor_location),
        ):

            locations = _norm_locations(loc_raw) or {None}
            uniprots = mapping.map_name(symbol, 'genesymbol', 'uniprot')

            for uniprot, location in itertools.product(uniprots, locations):

                result[uniprot].add(
                    Connectomedb2025Annotation(role = role, location = location)
                )

    return dict(result)
