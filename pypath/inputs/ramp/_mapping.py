#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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

import pandas as pd

from ._sqlite import ramp_sqlite


def ramp_mapping(
        id_type_a: str,
        id_type_b: str,
        return_df: bool = False,
        curies: bool = False,
    ) -> dict[str, set[str]] | pd.DataFrame:
    """
    Retrieve the mapping between two identifiers.

    Args:
        id_type_a:
            The identifier type of the first identifier.
        id_type_b:
            The identifier type of the second identifier.
        return_df:
            Return a pandas DataFrame instead of a dictionary.
        curies:
            Do not remove CURIEs from the identifiers.

    Returns:
        A dictionary with the mapping between the two identifiers.
    """

    query = (
        'SELECT DISTINCT a.sourceId as id_type_a, b.sourceId as id_type_b '
        'FROM '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_a}") a '
        'JOIN '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_b}") b '
        'ON a.rampId = b.rampId;'
    )

    con = ramp_sqlite()
    df = pd.read_sql_query(query, con)
    con.close()

    if not curies:

        df[df.columns] = df[df.columns].apply(
            lambda y: [x.split(':', maxsplit = 1)[-1] for x in y],
        )

    return (
        df
            if return_df else
        df.groupby('id_type_a')['id_type_b'].apply(set).to_dict()
    )


def ramp_synonym_mapping(
        id_type: str,
        return_df: bool = False,
        curies: bool = False,
    ) -> dict[str, set[str]] | pd.DataFrame:
    """
    Retrieve the mapping between an identifier type and synonyms.

    Joins the ``source`` and ``analytesynonym`` tables on RaMP IDs to
    produce an ``<ID type> -> synonym`` translation dictionary.

    Args:
        id_type:
            The identifier type, e.g. ``'hmdb'``, ``'swisslipids'``,
            ``'chebi'``, etc.
        return_df:
            Return a pandas DataFrame instead of a dictionary.
        curies:
            Do not remove CURIEs from the identifiers.

    Returns:
        A dictionary mapping identifiers to sets of synonyms, or a
        pandas DataFrame with ``source_id`` and ``synonym`` columns.
    """

    query = (
        'SELECT DISTINCT s.sourceId AS source_id, a.Synonym AS synonym '
        'FROM source s '
        'JOIN analytesynonym a ON s.rampId = a.rampId '
        f'WHERE s.geneOrCompound = "compound" AND s.IDtype = "{id_type}"'
    )

    con = ramp_sqlite()
    df = pd.read_sql_query(query, con)
    con.close()

    if not curies:

        df['source_id'] = df['source_id'].apply(
            lambda x: x.split(':', maxsplit = 1)[-1],
        )

    return (
        df
            if return_df else
        df.groupby('source_id')['synonym'].apply(set).to_dict()
    )


def ramp_synonyms_chebi() -> dict[str, str]:
    """
    Compound synonym → ChEBI ID mapping derived from RaMP.

    Inverts :func:`ramp_synonym_mapping` for ``'chebi'`` to produce a
    lookup dict keyed by lowercase synonym strings.  RaMP aggregates
    synonyms from HMDB, ChEBI, KEGG, WikiPathways, and Reactome, so
    this provides broader coverage than any single source.

    Returns:
        Dict mapping lowercase synonym strings to ChEBI ID strings
        (e.g. ``'adenosine triphosphate'`` → ``'CHEBI:30616'``).
        When multiple ChEBI IDs share a synonym, the first encountered
        mapping is kept.
    """

    chebi_to_syns = ramp_synonym_mapping('chebi', curies = True)

    result: dict[str, str] = {}

    for chebi_curie, syns in chebi_to_syns.items():

        # curies=True gives 'chebi:15422'; normalise to 'CHEBI:15422'
        chebi_id = str(chebi_curie).upper()

        for syn in syns:

            if syn:
                result.setdefault(str(syn).lower(), chebi_id)

    return result
