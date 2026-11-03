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

"""
MetaNetX/MNXref metabolite cross-reference mappings.

MetaNetX (https://www.metanetx.org) is a reconciled namespace for
genome-scale metabolic networks.  The MNXref cross-reference file maps
external database identifiers (ChEBI, HMDB, KEGG, etc.) to MetaNetX
compound identifiers (``MNXM*``).

Primary use-case in pypath: bridging database-specific compound IDs
(e.g. BiGG base metabolite IDs via their ``metanetx.chemical``
annotations) to ChEBI identifiers when direct cross-references are
unavailable.
"""

from __future__ import annotations

__all__ = ['metanetx_metabolite_chebi']

import io

import pypath.resources.urls as urls
import pypath.share.session as session_mod
from pypath.share.downloads import dm, _resolve_data_dir


_logger = session_mod.Logger(name='inputs.metanetx')


def metanetx_metabolite_chebi() -> dict[str, str]:
    """
    Build a MetaNetX compound ID (MNXM) → ChEBI ID mapping.

    Downloads the MNXref chemical cross-reference file
    (``chem_xref.tsv``) once and caches it to disk under
    ``<cachedir>/metanetx/chem_xref.tsv``.  Subsequent calls in any
    process read the cached file without re-downloading.

    The cross-reference file has tab-separated columns
    ``source_id``, ``mnx_id``, ``description``.  Only rows where
    ``source_id`` starts with ``chebi:`` are used; all other
    external sources are ignored.  When a MetaNetX ID maps to
    multiple ChEBI entries, the first encountered is kept.

    Returns:
        Dict mapping MetaNetX compound IDs (e.g. ``'MNXM15'``) to
        ChEBI IDs (e.g. ``'CHEBI:30616'``).
    """

    url = urls.urls['metanetx']['chem_xref']
    local_path = _resolve_data_dir() / 'metanetx' / 'chem_xref.tsv'
    local_path.parent.mkdir(parents=True, exist_ok=True)

    if not local_path.exists():
        try:
            dm.download(url, dest=str(local_path))
        except Exception as exc:
            _logger._log(
                f'MetaNetX: download failed ({exc}); returning empty mapping.'
            )
            return {}

    try:
        raw = local_path.read_text(encoding='utf-8')
    except Exception as exc:
        _logger._log(f'MetaNetX: could not read cache ({exc}); returning empty mapping.')
        return {}

    result: dict[str, str] = {}

    for line in io.StringIO(raw):

        line = line.strip()

        if not line or line.startswith('#'):
            continue

        fields = line.split('\t')

        if len(fields) < 2:
            continue

        source_id = fields[0].strip()
        mnx_id = fields[1].strip()

        if (
            not source_id.startswith('chebi:') or
            not mnx_id.startswith('MNXM') or
            mnx_id in result
        ):
            continue

        result[mnx_id] = source_id.replace('chebi:', 'CHEBI:')

    _logger._log(
        f'MetaNetX: loaded {len(result):,} MNXM → ChEBI entries '
        f'from MNXref chem_xref.tsv.'
    )

    return result


# ------------------------------------------------------------------
# General cross-reference functions
# ------------------------------------------------------------------

# Map chem_xref.tsv source_id prefixes to canonical id_type names.
_PREFIX_TO_ID_TYPE: dict[str, str] = {
    "chebi:": "chebi",
    "hmdb:": "hmdb",
    "kegg.compound:": "kegg",
    "kegg.drug:": "kegg",
    "bigg.metabolite:": "bigg",
    "seed.compound:": "seed",
    "lipidmaps:": "lipidmaps",
    "metacyc.compound:": "metacyc",
    "reactome:": "reactome",
    "swisslipids:": "swisslipids",
}


def _strip_prefix(source_id: str) -> tuple[str | None, str]:
    """Return (id_type, bare_id) for a chem_xref source_id string."""

    for prefix, id_type in _PREFIX_TO_ID_TYPE.items():
        if source_id.startswith(prefix):
            bare = source_id[len(prefix):]

            if id_type == "chebi":
                bare = f"CHEBI:{bare}"

            return id_type, bare

    return None, source_id


def _load_chem_xref() -> str:
    """Download or read cached chem_xref.tsv, return raw text."""

    url = urls.urls["metanetx"]["chem_xref"]
    local_path = _resolve_data_dir() / "metanetx" / "chem_xref.tsv"
    local_path.parent.mkdir(parents=True, exist_ok=True)

    if not local_path.exists():
        dm.download(url, dest=str(local_path))

    return local_path.read_text(encoding="utf-8")


def metanetx_mapping(
    source_type: str,
    target_type: str,
) -> dict[str, set[str]]:
    """Pairwise metabolite ID mapping via MetaNetX bridge.

    Parses ``chem_xref.tsv`` and extracts all pairs where an MNXM
    entry has both *source_type* and *target_type* external references.
    For example, ``metanetx_mapping("bigg", "chebi")`` returns
    ``{"atp": {"CHEBI:30616"}, ...}``.

    Supported id_type values: ``chebi``, ``hmdb``, ``kegg``, ``bigg``,
    ``seed``, ``lipidmaps``, ``metacyc``, ``reactome``, ``swisslipids``.

    Args:
        source_type: Canonical id_type name for the source.
        target_type: Canonical id_type name for the target.

    Returns:
        Dict mapping source IDs to sets of target IDs.
    """

    from collections import defaultdict

    text = _load_chem_xref()

    # First pass: collect all external IDs per MNXM entry
    mnx: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )

    for line in io.StringIO(text):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")

        if len(fields) < 2:
            continue

        source_id = fields[0].strip()
        mnx_id = fields[1].strip()

        if not mnx_id.startswith("MNXM"):
            continue

        id_type, bare = _strip_prefix(source_id)

        if id_type is not None and id_type in (source_type, target_type):
            mnx[mnx_id][id_type].add(bare)

    # Second pass: build pairwise mapping
    result: dict[str, set[str]] = defaultdict(set)

    for _mnx_id, xrefs in mnx.items():
        src_ids = xrefs.get(source_type)
        tgt_ids = xrefs.get(target_type)

        if not src_ids or not tgt_ids:
            continue

        for sid in src_ids:
            result[sid].update(tgt_ids)

    _logger._log(
        f"MetaNetX: {source_type} -> {target_type}: "
        f"{len(result):,} entries from MNXref chem_xref.tsv."
    )

    return dict(result)
