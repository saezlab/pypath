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
CollecTRI2: TF-target transcriptional regulatory interactions.

Two files are combined:

- the *combined* file provides, per interaction, the direction (TF -> target),
  the consensus sign (+1 activation / -1 repression) and the TF type;
- the *main* (detail) file provides, per interaction, which primary sources
  support it (with their PMIDs and source-specific variables) and the TF
  classification flags.

The primary sources (ExTRI2, HTRI, TRRUST, TFactS, GOA, IntAct, SIGNOR, CytReg,
GEREDB, NTNUcuration, Pavlidis2021, DoRothEA-A) are recorded as secondary
(`via`) resources of CollecTRI2; their PMIDs become the interaction references
and their source-specific variables are kept as an `evidence` edge attribute.
"""

from __future__ import annotations

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session

_log = session.Logger(name = 'collectri2_input')._log


# primary source -> column indices (0-based) in the main (detail) file:
#   `present` flag column, `pmid` column, and source-specific `vars`
_SOURCES = {
    'ExTRI2': {'present': 3, 'pmid': 4, 'vars': {'sign': 5, 'tf_type': 6}},
    'HTRI': {'present': 7, 'pmid': 9, 'vars': {'technique': 8, 'confidence': 10}},
    'TRRUST': {'present': 11, 'pmid': 13, 'vars': {'regulation': 12}},
    'TFactS': {
        'present': 14, 'pmid': 18,
        'vars': {'sign': 15, 'species': 16, 'source': 17, 'confidence': 19},
    },
    'GOA': {'present': 20, 'pmid': 22, 'vars': {'sign': 21}},
    'IntAct': {'present': 23, 'pmid': 24, 'vars': {'method_id': 25}},
    'SIGNOR': {'present': 26, 'pmid': 29, 'vars': {'effect': 27, 'sign': 28}},
    'CytReg': {
        'present': 30, 'pmid': 34,
        'vars': {
            'assay_type': 31,
            'species': 32,
            'activation_repression': 33,
            'year': 35,
        },
    },
    'GEREDB': {'present': 36, 'pmid': 38, 'vars': {'effect': 37}},
    'NTNUcuration': {'present': 39, 'pmid': 41, 'vars': {'sign': 40}},
    'Pavlidis2021': {'present': 42, 'pmid': 43, 'vars': {'mor': 44}},
    'DoRothEA-A': {'present': 45, 'pmid': 46, 'vars': {'directed': 47, 'effect': 48}},
}

# TF-classification columns (0-based) in the main file
_TFCLASS_COLS = {
    'lambert': 49,
    'lovering': 50,
    'go_0003700': 51,
    'go_0140223': 52,
    'go_0003712': 53,
    'tfclass': 54,
}

_AUTOREG_COL = 55

_N_MAIN_COLS = 56


CollecTRI2Interaction = collections.namedtuple(
    'CollecTRI2Interaction',
    (
        'tf',               # 0: TF gene symbol
        'target',           # 1: target gene symbol
        'sign',             # 2: +1 / -1 (consensus sign from the combined file)
        'tf_type',          # 3: dbTF / coTF / coTF candidate
        'resources',        # 4: ';'-joined secondary (via) source names
        'references',       # 5: ';'-joined PubMed IDs
        'auto_regulation',  # 6: auto-regulation flag
        'evidence',         # 7: dict of per-source variables + TF classification
    ),
)


def collectri2_raw() -> list[CollecTRI2Interaction]:
    """
    Parses the CollecTRI2 combined and main (detail) files and joins them on
    the (TF, target) pair, returning records with gene symbols (not yet mapped
    to UniProt).
    """

    # combined file: ID, source (TF), target, weight, TF Type, Auto-regulation
    c_comb = curl.Curl(
        urls.urls['collectri2']['combined'],
        silent = False,
        large = True,
    )

    combined = {}

    for line in c_comb.result:

        if not line or line.startswith('#'):

            continue

        f = line.rstrip('\n').split('\t')

        if len(f) < 5:

            continue

        combined[(f[1], f[2])] = (
            int(f[3]) if f[3] in ('1', '-1') else 0,
            f[4] or None,
        )

    # main (detail) file
    c_main = curl.Curl(
        urls.urls['collectri2']['main'],
        silent = False,
        large = True,
    )

    result = []
    n_unmatched = 0

    for line in c_main.result:

        if not line or line.startswith('#'):

            continue

        f = line.rstrip('\n').split('\t')

        if len(f) < _N_MAIN_COLS:

            continue

        tf, target = f[1], f[2]

        sources = []
        pmids = set()
        evidence = {}

        for name, spec in _SOURCES.items():

            if not f[spec['present']]:

                continue

            sources.append(name)

            pmid_val = f[spec['pmid']]

            if pmid_val:

                pmids.update(
                    p.strip() for p in pmid_val.split('|') if p.strip()
                )

            ev = {
                vname: f[vidx]
                for vname, vidx in spec['vars'].items()
                if f[vidx]
            }

            if ev:

                evidence[name] = ev

        tfclass = {
            k: f[idx]
            for k, idx in _TFCLASS_COLS.items()
            if f[idx]
        }

        if tfclass:

            evidence['tf_classification'] = tfclass

        sign, tf_type = combined.get((tf, target), (0, None))

        if (tf, target) not in combined:

            n_unmatched += 1

        result.append(
            CollecTRI2Interaction(
                tf = tf,
                target = target,
                sign = sign,
                tf_type = tf_type,
                resources = ';'.join(sources),
                references = ';'.join(sorted(pmids)),
                auto_regulation = f[_AUTOREG_COL] or None,
                evidence = evidence,
            )
        )

    if n_unmatched:

        _log(
            f'CollecTRI2: {n_unmatched} interactions from the main file had no '
            'match in the combined file; their sign is set to unknown.'
        )

    _log(f'CollecTRI2: {len(result)} interactions loaded.')

    return result


def collectri2_interactions() -> list[CollecTRI2Interaction]:
    """
    TF-target interactions from CollecTRI2, with gene symbol identifiers.

    Identifier translation to UniProt is performed by the network reader
    (the corresponding `NetworkInput` uses `id_type = 'genesymbol'`).
    """

    return collectri2_raw()


def collectri2_annotations() -> dict:
    """
    Transcription factor type annotations from CollecTRI2.

    Returns a dict mapping UniProt IDs to sets of `CollecTRI2Annotation` named
    tuples carrying the TF type (`dbTF`, `coTF`, `coTF candidate`). The TF type
    is consistent per TF in the source.
    """

    import pypath.utils.mapping as mapping

    CollecTRI2Annotation = collections.namedtuple(
        'CollecTRI2Annotation',
        ('tf_type',),
    )

    result = collections.defaultdict(set)
    seen = set()

    for rec in collectri2_raw():

        if not rec.tf_type or rec.tf in seen:

            continue

        seen.add(rec.tf)

        for uniprot in mapping.map_name(rec.tf, 'genesymbol', 'uniprot'):

            result[uniprot].add(CollecTRI2Annotation(tf_type = rec.tf_type))

    return dict(result)
