"""MetaNetX chemical identities, structures, and external cross-references."""

from __future__ import annotations

from collections.abc import Iterable

from ..schema import CHEMICAL_TAXON
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv"
PROPERTIES_URL = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv"

_PREFIX = {
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
    for prefix, id_type in _PREFIX.items():
        if source_id.startswith(prefix):
            bare = source_id[len(prefix) :]
            if id_type == "chebi":
                bare = f"CHEBI:{bare}"
            return id_type, bare
    return None, source_id


def emit(
    writer: HubParquetWriter,
    lines: Iterable[str] | None = None,
    property_lines: Iterable[str] | None = None,
) -> None:
    seen: set[str] = set()
    # Explicit fixture lines stay offline unless property lines are supplied.
    properties = (
        iter_url_lines(PROPERTIES_URL)
        if lines is None and property_lines is None
        else property_lines
        if property_lines is not None
        else ()
    )
    for line in properties:
        if writer.full:
            return
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 9 or not (fields[0].startswith("MNXM") or fields[0] == "WATER"):
            continue
        mnx_id = fields[0].strip()
        if mnx_id not in seen:
            seen.add(mnx_id)
            if not writer.add_identity("metanetx", mnx_id, CHEMICAL_TAXON, "metanetx"):
                return
        # chem_prop structures are MetaNetX-standardized, not a claim that all
        # external cross-references have identical protonation/stereochemistry.
        for position, namespace in ((1, "name"), (6, "inchi"), (7, "inchikey"), (8, "smiles")):
            value = fields[position].strip()
            if value and value not in ("NA", "-"):
                if not writer.add(namespace, value, mnx_id, CHEMICAL_TAXON, "metanetx"):
                    return
    for line in lines if lines is not None else iter_url_lines(URL):
        if writer.full:
            return
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        source_id = fields[0].strip()
        mnx_id = fields[1].strip()
        if not (mnx_id.startswith("MNXM") or mnx_id == "WATER"):
            continue
        if mnx_id not in seen:
            seen.add(mnx_id)
            if not writer.add_identity("metanetx", mnx_id, CHEMICAL_TAXON, "metanetx"):
                return
        source_type, bare = _strip_prefix(source_id)
        if not source_type or not bare:
            continue
        if source_type == "metanetx" and bare == mnx_id:
            continue
        if not writer.add(source_type, bare, mnx_id, CHEMICAL_TAXON, "metanetx"):
            return
