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

from __future__ import annotations

import collections

__all__ = [
    'G2PInteraction',
    'G2PLigand',
    'G2PTarget',
]


def _record_factory(
        name: str,
        fields: list[str],
        defaults: tuple = (),
    ) -> collections.namedtuple:

    RecordBase = collections.namedtuple(name, fields)
    RecordBase.__new__.__defaults__ = defaults

    class Record(RecordBase):

        def __new__(cls, *args, **kwargs):
            """
            Makes sure empty strings are converted to `None`.
            """

            args = [arg if arg != '' else None for arg in args]
            kwargs = {k: (v if v != '' else None) for k, v in kwargs.items()}

            return super().__new__(cls, *args, **kwargs)

    return Record


G2PInteraction = _record_factory(
    "G2PInteraction",
    [
        "ligand",
        "target",
        "action",
        "ligand_type",
        "is_stimulation",
        "is_inhibition",
        "endogenous",
        "affinity_high",
        "affinity_low",
        "affinity_median",
        "affinity_units",
        "primary_target",
        "pubmed",
    ],
)

G2PLigand = _record_factory(
    "G2PLigand",
    [
        "name",
        "uniprot",
        "pubchem",
        "iupac",
        "chembl",
        "smiles",
        "inchi",
        "organism",
        "entity_type",
        "subtype",
        "role",
    ],
    ('ligand',),
)


G2PTarget = _record_factory(
    "G2PTarget",
    [
        "uniprot",
        "symbol",
        "entrez",
        "ensembl",
        "refseq",
        "refseqp",
        "family",
        "target_type",
        "organism",
        "entity_type",
        "role",
    ],
    ('target',),
)
