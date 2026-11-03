"""MetaNetX reaction hub: one record per MNXR reaction with the reaction ids MetaNetX reconciles.

MetaNetX lists Rhea master and directional ids; every Rhea id is written as its master so that
the one-step attach rule reaches the Rhea record (whose record id is the master). Only the
long source prefixes are read (``bigg.reaction:`` rather than the duplicate ``biggR:``).
"""

from __future__ import annotations

from collections.abc import Iterable

from ..schema import CHEMICAL_TAXON
from ..stream import iter_url_lines
from ..writer import HubParquetWriter
from .rhea import DIRECTIONS_URL, _rows

URL = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv"

PREFIX = {
    "rhea:": "rhea",
    "kegg.reaction:": "kegg_reaction",
    "metacyc.reaction:": "metacyc_reaction",
    "bigg.reaction:": "bigg_reaction",
    "seed.reaction:": "seed_reaction",
    "vmhreaction:": "vmh_reaction",
    "sabiork.reaction:": "sabiork_reaction",
}


def rhea_masters(lines: Iterable[str]) -> dict[str, str]:
    """Every Rhea id (master and directional) -> its master."""
    master_of = {}
    for fields in _rows(lines):
        if len(fields) >= 4 and fields[0].strip():
            for value in fields[:4]:
                if value.strip():
                    master_of[value.strip()] = fields[0].strip()
    return master_of


def emit(
    writer: HubParquetWriter,
    lines: Iterable[str] | None = None,
    directions: Iterable[str] | None = None,
) -> None:
    master_of = rhea_masters(
        directions if directions is not None else iter_url_lines(DIRECTIONS_URL)
    )
    seen: set[str] = set()
    for line in lines if lines is not None else iter_url_lines(URL):
        if writer.full:
            return
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 2:
            continue
        source, mnx_id = fields[0].strip(), fields[1].strip()
        if not mnx_id.startswith("MNXR"):
            continue  # EMPTY and other pseudo reactions
        if mnx_id not in seen:
            seen.add(mnx_id)
            if not writer.add_identity("metanetx_reaction", mnx_id, CHEMICAL_TAXON, "metanetx"):
                return
        for prefix, source_type in PREFIX.items():
            if source.startswith(prefix):
                value = source[len(prefix) :].strip()
                if source_type == "rhea":
                    value = master_of.get(value, value)
                if value and not writer.add(source_type, value, mnx_id, CHEMICAL_TAXON, "metanetx"):
                    return
                break
