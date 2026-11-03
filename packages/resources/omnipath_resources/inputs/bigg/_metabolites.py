"""BiGG Models metabolite cross-reference parser."""

from __future__ import annotations

__all__ = ["bigg_metabolite_mapping"]

import io
import re
from collections import defaultdict

import logging

import pypath.resources.urls as urls
from pypath.share.downloads import dm, _resolve_data_dir


_log = logging.getLogger(__name__)

# Mapping from BiGG database_links label to our canonical id_type names
_DB_LABEL_TO_ID_TYPE: dict[str, str] = {
    "CHEBI": "chebi",
    "Human Metabolome Database": "hmdb",
    "KEGG Compound": "kegg",
    "KEGG Drug": "kegg",
    "KEGG Glycan": "kegg_glycan",
    "MetaNetX (MNX) Chemical": "metanetx",
    "SEED Compound": "seed",
    "LipidMaps": "lipidmaps",
    "BioCyc": "biocyc",
    "Reactome Compound": "reactome",
    "InChI Key": "inchikey",
}

# Regex to extract the database ID from a URL
_RE_URL_ID = re.compile(r".*/([^/]+)$")


def _parse_db_links(db_links_str: str) -> dict[str, set[str]]:
    """Parse the semicolon-separated database_links field.

    Format: ``"DB_NAME: url; DB_NAME: url; ..."``
    Returns ``{id_type: {id_value, ...}}``.
    """

    result: dict[str, set[str]] = defaultdict(set)

    if not db_links_str:
        return result

    for entry in db_links_str.split(";"):
        entry = entry.strip()

        if ": " not in entry:
            continue

        db_name, url = entry.split(": ", 1)
        db_name = db_name.strip()
        url = url.strip()

        id_type = _DB_LABEL_TO_ID_TYPE.get(db_name)

        if id_type is None:
            continue

        m = _RE_URL_ID.match(url)

        if not m:
            continue

        result[id_type].add(m.group(1))

    return result


def bigg_metabolite_mapping(
    target_type: str = "chebi",
) -> dict[str, set[str]]:
    """Build a BiGG base metabolite ID -> target ID mapping.

    Downloads the BiGG universal metabolite TSV once and caches it to
    disk.  Parses cross-references from the ``database_links`` column.

    Args:
        target_type: Target ID type (e.g. ``"chebi"``, ``"hmdb"``,
            ``"kegg"``, ``"metanetx"``).

    Returns:
        Dict mapping BiGG base IDs (e.g. ``"atp"``) to sets of target
        IDs (e.g. ``{"CHEBI:30616", "CHEBI:15422"}``).
    """

    url = urls.urls["bigg"]["metabolites_tsv"]
    local_path = _resolve_data_dir() / "bigg" / "bigg_models_metabolites.txt"
    local_path.parent.mkdir(parents=True, exist_ok=True)

    if not local_path.exists():
        try:
            dm.download(url, dest=str(local_path))
        except Exception as exc:
            _log.warning("BiGG: download failed (%s); returning empty mapping.", exc)
            return {}

    try:
        text = local_path.read_text(encoding="utf-8")
    except Exception as exc:
        _log.warning("BiGG: could not read cache (%s); returning empty mapping.", exc)
        return {}

    result: dict[str, set[str]] = defaultdict(set)

    lines = text.strip().split("\n")

    # Skip header
    for line in lines[1:]:
        fields = line.split("\t")

        if len(fields) < 5:
            continue

        universal_id = fields[1].strip()     # base ID (e.g. "atp")
        db_links_str = fields[4].strip()     # semicolon-separated cross-refs

        if not universal_id or not db_links_str:
            continue

        xrefs = _parse_db_links(db_links_str)
        target_ids = xrefs.get(target_type)

        if target_ids:
            result[universal_id].update(target_ids)

    _log.info(
        "BiGG: %d universal metabolites -> %s from BiGG Models TSV.",
        len(result), target_type,
    )

    return dict(result)
