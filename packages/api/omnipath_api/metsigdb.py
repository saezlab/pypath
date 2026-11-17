"""The contract of ``/sets/metsigdb``: what a request may ask and how it is read.

The route validates a request with these rules and gives the step files
(``metsigdb_rows``, ``metsigdb_groups``, ``metsigdb_sets``) their parameters. The rules
belong here, with the field registry; the status codes belong to the route.
"""

from __future__ import annotations

from collections.abc import Iterable
import re

# The v1 enums. A filter value is matched to its published spelling at the edge, which is
# what makes the match case-insensitive without lowering a column.
RESOURCES = ("KEGG", "Reactome", "WikiPathways", "MACdb", "ClassyFire")
SET_TYPES = ("disease", "pathway", "chemical_class")
SET_SUB_TYPES = (
    "cancer",
    "phenotype",
    "medical intervention",
    "gene abnormality",
    "genotype",
    "overview_map",
    "metabolic_map",
)
# The classes of the exclusion product. A class never changes meaning, so the list only grows.
EXCLUSION_CLASSES = (
    "ions",
    "atoms",
    "small_molecules",
    "xenobiotics",
    "cofactors",
    "cofactors_extended",
)

DEFAULT_LIMIT = 1000
MAX_LIMIT = 100_000

# The row contract of the subset this route reads (``subset_layout.ROW_CONTRACTS``).
ROW_CONTRACT = 2

# Everything a response row may carry, in the order of the row contract, under the names
# the filters use: `entity` for the metabolite, `set` for the set's source id.
ALL_FIELDS = (
    "entity",
    "metabolite_label",
    "metabolite_entity_type",
    "metabolite_structure_key",
    "inchikey",
    "smiles",
    "hmdb",
    "pubchem",
    "chebi",
    "kegg",
    "has_structure",
    "exclusion_classes",
    "resource",
    "set",
    "set_entity",
    "set_label",
    "set_type",
    "set_sub_type",
    "organism",
    "set_size",
    "set_context",
    "provenance_source",
    "provenance_record",
    "build_id",
)
DEFAULT_FIELDS = (
    "entity",
    "metabolite_label",
    "inchikey",
    "hmdb",
    "pubchem",
    "chebi",
    "kegg",
    "resource",
    "set",
    "set_label",
    "set_type",
    "set_sub_type",
    "set_size",
)
ALL_KEYWORD = "all"

# A grouped response writes the set side once per set and the member side once per member.
# `set_context` describes a member's way into the set (a ClassyFire class is reached
# directly or through a descendant), so it is a member field.
SET_FIELDS = (
    "resource",
    "set",
    "set_entity",
    "set_label",
    "set_type",
    "set_sub_type",
    "organism",
    "set_size",
)
MEMBER_FIELDS = tuple(field for field in ALL_FIELDS if field not in SET_FIELDS)


def canonical(allowed: tuple[str, ...], value: str) -> str | None:
    """The published spelling of a vocabulary value, whatever case it arrived in."""
    for published in allowed:
        if published.lower() == value.lower():
            return published
    return None


def canonical_values(values: Iterable[str], allowed: tuple[str, ...], name: str) -> list[str]:
    """Each value in its published spelling; ValueError names what is not published."""
    values = list(values)
    found = [canonical(allowed, value) for value in values]
    unknown = sorted({value for value, spelled in zip(values, found) if spelled is None})
    if unknown:
        raise ValueError(
            f"Unsupported {name}: {', '.join(unknown)}. Supported: {', '.join(allowed)}."
        )
    return list(dict.fromkeys(found))


def exclusion_classes(values: Iterable[str]) -> list[str]:
    """The classes an ``exclude_metabolites`` request names, in published spelling.

    ``none`` and an absent parameter exclude nothing; ``all`` is every class. Neither may be
    combined with another value: a request that says both has not said what it means.
    """
    valid = ("none", "all", *EXCLUSION_CLASSES)
    named = list(dict.fromkeys(value.strip().lower() for value in values if value.strip()))
    unknown = [value for value in named if value not in valid]
    if unknown:
        raise ValueError(
            f"Unsupported exclude_metabolites: {', '.join(sorted(set(unknown)))}. "
            f"Supported: {', '.join(valid)}."
        )
    if not named or named == ["none"]:
        return []
    if "none" in named or "all" in named:
        if named == ["all"]:
            return list(EXCLUSION_CLASSES)
        raise ValueError(
            "exclude_metabolites: none and all stand alone and cannot be combined with another value."
        )
    return named


def resolve_fields(requested: Iterable[str] | None) -> tuple[str, ...]:
    """The fields one response carries, in contract order.

    ``fields`` adds to the default rather than replacing it, so a request naming one extra
    field does not lose the twelve it did not mention. ``all`` is every field. An unknown
    name raises: ignoring it would answer a different question with a 200.
    """
    requested = tuple(requested or ())
    if not requested:
        return DEFAULT_FIELDS
    if ALL_KEYWORD in requested:
        return ALL_FIELDS
    unknown = sorted({name for name in requested if name not in ALL_FIELDS})
    if unknown:
        raise ValueError(
            f"Unsupported field: {', '.join(unknown)}. "
            f'Supported: {", ".join(ALL_FIELDS)}, or "{ALL_KEYWORD}".'
        )
    wanted = set(DEFAULT_FIELDS) | set(requested)
    return tuple(field for field in ALL_FIELDS if field in wanted)


# The shape that names each identifier a metabolite may be given by. A caller pastes the
# identifier they hold and the route works out the namespace: a type parameter would make
# the easy case harder and the hard case no easier. SMILES is absent on purpose: it has no
# distinguishing shape, so a rule for it would swallow every value the others declined.
_ENTITY_KEY = re.compile(r"^[0-9a-f]{64}$", re.I)
_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", re.I)
_HMDB = re.compile(r"^HMDB(\d{5,11})$", re.I)
_CHEBI = re.compile(r"^CHEBI:(\d+)$", re.I)
_KEGG = re.compile(r"^C\d{5}$", re.I)
_INTEGER = re.compile(r"^\d+$")

IDENTIFIER_PARAMETERS = ("entity_keys", "inchikeys", "hmdb", "chebi", "kegg", "pubchem")


def entity_parameters(values: Iterable[str]) -> dict[str, list[str]]:
    """The step-file parameters that select the metabolites named by ``values``.

    A bare integer is two identifiers: ChEBI and PubChem numbers are both bare integers, and
    the value alone cannot decide between them. The union is the honest answer: guessing one
    namespace would drop the other's memberships from a result the caller believes complete.

    Raises ValueError for a value no rule claims, rather than answering with an empty page
    that reads as "no such metabolite".
    """
    found: dict[str, list[str]] = {name: [] for name in IDENTIFIER_PARAMETERS}
    for value in (v.strip() for v in values):
        if _ENTITY_KEY.match(value):
            found["entity_keys"].append(value.lower())
        elif _INCHIKEY.match(value):
            found["inchikeys"].append(value.upper())
        elif match := _HMDB.match(value):
            # Stored in the seven-digit form.
            found["hmdb"].append("HMDB" + match.group(1).lstrip("0").rjust(7, "0"))
        elif match := _CHEBI.match(value):
            found["chebi"].append(match.group(1))
        elif _KEGG.match(value):
            found["kegg"].append(value.upper())
        elif _INTEGER.match(value):
            found["chebi"].append(value)
            found["pubchem"].append(value)
        else:
            raise ValueError(
                f"Unrecognised entity identifier: {value}. Supported: an entity key, an "
                "InChIKey, or an HMDB, ChEBI, KEGG or PubChem identifier. SMILES is "
                "published but not searchable."
            )
    return {name: list(dict.fromkeys(items)) for name, items in found.items()}
