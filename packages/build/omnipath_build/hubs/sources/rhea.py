"""Rhea hub: master reactions, their directional ids, Rhea's own cross-references and equations.

One record per master reaction (the reaction without direction, the identity anchor). Its
left-to-right, right-to-left and bidirectional ids are ``rhea`` claims of the same record, so a
directional id resolves to its master. Cross-references are the ones Rhea itself asserts
(rhea2xrefs); EC numbers are kept as an attribute, not as an identity. The ``DEFINITION`` of
the master (an equation with participant names) is the record's ``name``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..schema import CHEMICAL_TAXON
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

BASE = "https://ftp.expasy.org/databases/rhea"
DIRECTIONS_URL = f"{BASE}/tsv/rhea-directions.tsv"
XREFS_URL = f"{BASE}/tsv/rhea2xrefs.tsv"
REACTIONS_URL = f"{BASE}/txt/rhea-reactions.txt.gz"

# Rhea database name -> hub source type.
XREF_TYPES = {
    "KEGG_REACTION": "kegg_reaction",
    "REACTOME": "reactome",
    "METACYC": "metacyc_reaction",
    "ECOCYC": "ecocyc_reaction",
    "MACIE": "macie",
    "EC": "ec",
}


def _rows(lines: Iterable[str]):
    """Tab-separated rows after the header line."""
    header = True
    for line in lines:
        if header:
            header = False
            continue
        line = line.rstrip("\r\n")
        if line:
            yield line.split("\t")


def emit(
    writer: HubParquetWriter,
    directions: Iterable[str] | None = None,
    xrefs: Iterable[str] | None = None,
    reactions: Iterable[str] | None = None,
) -> None:
    masters: set[str] = set()
    for fields in _rows(directions if directions is not None else iter_url_lines(DIRECTIONS_URL)):
        if len(fields) < 4 or not fields[0].strip():
            continue
        master = fields[0].strip()
        masters.add(master)
        if not writer.add_identity("rhea", master, CHEMICAL_TAXON, "rhea"):
            return
        for directional in fields[1:4]:
            directional = directional.strip()
            if directional and directional != master:
                if not writer.add("rhea", directional, master, CHEMICAL_TAXON, "rhea"):
                    return
    # rhea2xrefs: RHEA_ID, DIRECTION, MASTER_ID, ID, DB
    for fields in _rows(xrefs if xrefs is not None else iter_url_lines(XREFS_URL)):
        if writer.full:
            return
        if len(fields) < 5:
            continue
        master, value, source_type = (
            fields[2].strip(),
            fields[3].strip(),
            XREF_TYPES.get(fields[4].strip()),
        )
        if master in masters and value and source_type:
            if not writer.add(source_type, value, master, CHEMICAL_TAXON, "rhea"):
                return
    # rhea-reactions.txt: ENTRY RHEA:n / DEFINITION <equation with names> / ... / '///'
    entry = None
    for line in reactions if reactions is not None else iter_url_lines(REACTIONS_URL):
        if writer.full:
            return
        if line.startswith("ENTRY"):
            entry = line.split(None, 1)[1].strip().removeprefix("RHEA:") if " " in line else None
        elif line.startswith("DEFINITION") and entry in masters:
            definition = line.split(None, 1)[1].strip() if " " in line.strip() else ""
            if definition and not writer.add("name", definition, entry, CHEMICAL_TAXON, "rhea"):
                return
        elif line.startswith("///"):
            entry = None
