"""Lipid name nodes and their generalization edges, from the identities a release publishes.

A lipid with no structure is identified by its Goslin name: ``goslin:<level>:<name>``, the
level at which the name was anchored (species, molecular species, sn position, ...). This
reads those identities, parses each name with the nomenclature grammar for the class, the
totals and the chains, and derives the graph: every name below species is a kind of its own
species. Nothing is asserted by a resource, and no service is read.

**Chain counts come from the parse object, not from counting ``C:D`` tokens in a string**: a
one-chain class renders empty slots (``MG 18:1/0:0/0:0``) and a sterol ester renders its
backbone as a chain (``CE 18:1`` is ``SE 27:1/18:1``). A chain is listed when its slot has
carbon. A name with some but not all chains listed (``TG 58:12_20:0``) is its own level,
``partially_specified``, between species and a fully resolved molecular species.
"""

from __future__ import annotations

from collections.abc import Iterable
from functools import cache
import re

PARTIALLY_SPECIFIED = "partially_specified"
_CHAIN_TOKEN = re.compile(r"\d+:\d+")
IDENTITY_PREFIX = "goslin:"


@cache
def _parser():
    from pygoslin.parser.Parser import LipidParser

    return LipidParser()


def _version() -> str:
    try:
        import importlib.metadata

        return importlib.metadata.version("pygoslin")
    except Exception:  # noqa: BLE001
        return ""


def parse_identity(reference_key: str | None) -> dict | None:
    """The node a published ``goslin:<level>:<name>`` identity stands for, with the fields of
    ``lipid_name_node``; ``None`` for any other key or a name the grammar cannot read."""
    if not reference_key or not reference_key.startswith(IDENTITY_PREFIX):
        return None
    parts = reference_key.split(":", 2)
    if len(parts) != 3 or not parts[2] or not _CHAIN_TOKEN.search(parts[2]):
        return None
    name = parts[2]
    try:
        from pygoslin.domain.LipidLevel import LipidLevel

        parsed = _parser().parse(name)
        lipid = parsed.lipid
        info = lipid.info
        level = info.level.name.lower()
        possible = info.poss_fa
        listed = sum(1 for fa in lipid.fa_list if fa.num_carbon > 0)
        species = parsed.get_lipid_string(LipidLevel.SPECIES)
    except Exception:  # noqa: BLE001 -- the grammar raises a family of parse errors
        return None
    if isinstance(possible, int) and possible > 0 and 0 < listed < possible:
        level = PARTIALLY_SPECIFIED

    def attribute(read):
        try:
            return read()
        except Exception:  # noqa: BLE001
            return None

    total_carbon = total_db = None
    chain = _CHAIN_TOKEN.search(species or "")
    if chain:
        total_carbon, total_db = (int(x) for x in chain.group().split(":"))
    return {
        "lipid_name": name,
        "lipid_level": level,
        "chains_possible": possible or 0,
        "chains_listed": listed,
        "lipid_category": attribute(lambda: lipid.headgroup.lipid_category.name),
        "lipid_class": attribute(lambda: lipid.headgroup.get_lipid_string()),
        "total_carbon": total_carbon,
        "total_db": total_db,
        "sum_formula": attribute(parsed.get_sum_formula),
        "parser_version": _version(),
    }


def node_key(node: dict) -> tuple[str, str, int, int]:
    """The primary key of ``lipid_name_node``."""
    return (node["lipid_name"], node["lipid_level"], node["chains_listed"], node["chains_possible"])


def lipid_graph(parsed: Iterable[dict]) -> tuple[dict[tuple, dict], dict[tuple, dict]]:
    """The nodes and ``is_a`` edges of parsed names.

    Every node below species is a kind of its own species: same class, same totals, no chain
    listed. The species node is made when no resource named it, so the edge always has a
    parent to attach to.
    """
    nodes = {node_key(node): dict(node) for node in parsed}
    edges: dict[tuple, dict] = {}
    for key, node in list(nodes.items()):
        if node["lipid_level"] == "species" or node["chains_listed"] == 0:
            continue
        lipid_class, carbon, db = node["lipid_class"], node["total_carbon"], node["total_db"]
        if not lipid_class or carbon is None or db is None:
            continue
        parent_name = f"{lipid_class} {carbon}:{db}"
        parent_key = (parent_name, "species", 0, node["chains_possible"])
        if parent_key not in nodes:
            nodes[parent_key] = {
                **node,
                "lipid_name": parent_name,
                "lipid_level": "species",
                "chains_listed": 0,
                "sum_formula": None,
            }
        if parent_key == key:
            continue  # a species-total rendering that round-tripped to itself
        edges[(key, parent_key)] = {
            "child": key,
            "parent": parent_key,
            "relation": "is_a",
            "derivation": "species_generalization",
        }
    return nodes, edges
