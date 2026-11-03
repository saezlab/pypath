"""Resolve seed identifiers to published chemical entities.

The release already holds every identifier an entity carries, so a seed is
found by reading them, not by translating anything. The match is on the
identifier type and the value together. A bare number such as ``15378`` is a
ChEBI id, a PubChem id and sometimes a MACdb trait id, and it must find only the
entity of the namespace the seed names.

Values are not spelled alike across sources. ChEBI keeps its ``CHEBI:`` prefix
in some places and drops it in others, and HMDB ships a five- and a seven-digit
accession. Each seed is looked up under every spelling the release may hold.
"""

from __future__ import annotations

from dataclasses import dataclass
import re

from omnipath_subsets.exclusion_metabolites.seeds import Seed

# The identifier type name the release uses for each chemical namespace.
IDENTIFIER_TYPES = {
    "chebi": "Chebi:MI:0474",
    "kegg": "Kegg Compound:MI:2012",
    "hmdb": "Hmdb:OM:0004",
    "pubchem": "Pubchem Compound:OM:0002",
    "lipidmaps": "Lipidmaps:OM:0003",
    "swisslipids": "Swisslipids:OM:0009",
    "chembl": "Chembl Compound:MI:0967",
    "drugbank": "Drugbank:MI:2002",
}

# Entity types a seed may resolve to: the published chemical families.
CHEMICAL_ENTITY_TYPES = ("chemical_entity", "small_molecule")

# An entity whose canonical identity is in one of these namespaces was resolved
# against a chemical database. Any other canonical identity (a name, a source
# accession) is a native fallback: the build had nothing better to anchor on.
STANDARD_IDENTITY_TYPES = tuple(IDENTIFIER_TYPES.values()) + ("Standard Inchi Key:MI:1101",)


def spellings(namespace: str, identifier: str) -> list[str]:
    """Every form in which the release may hold ``identifier``."""
    text = identifier.strip()
    if namespace == "chebi":
        match = re.fullmatch(r"(?i)(?:CHEBI:)?(\d+)", text)
        if match:
            number = match.group(1)
            return [f"CHEBI:{number}", number]
    elif namespace == "hmdb":
        match = re.fullmatch(r"(?i)HMDB(\d+)", text)
        if match:
            number = int(match.group(1))
            forms = [f"HMDB{number:07d}"]
            if number <= 99999:
                forms.append(f"HMDB{number:05d}")
            return forms
    elif namespace == "kegg":
        return [text.upper()]
    return [text]


@dataclass(frozen=True)
class Resolution:
    """What a seed found in the release."""

    seed: Seed
    entity_ids: tuple[str, ...]
    fallback_only: bool

    @property
    def resolved(self) -> bool:
        return bool(self.entity_ids)


def _has_column(cur, table: str, column: str) -> bool:
    cur.execute(
        "SELECT EXISTS (SELECT 1 FROM pg_attribute WHERE attrelid = to_regclass(%s) "
        "AND attname = %s AND NOT attisdropped)",
        [table, column],
    )
    return cur.fetchone()[0]


def resolve_seeds(cur, seeds: list[Seed]) -> list[Resolution]:
    """Resolve every seed in one query; the result follows the order of ``seeds``."""
    wanted = [
        (index, IDENTIFIER_TYPES[seed.namespace], spelling)
        for index, seed in enumerate(seeds)
        for spelling in spellings(seed.namespace, seed.identifier)
    ]
    # A release that stores a normalized value matches on it as well, so a
    # spelling this module did not foresee still finds its entity.
    matches_normalized = _has_column(cur, "identifier_evidence", "value_normalized")
    value_test = (
        "(i.value = w.spelling OR i.value_normalized = w.spelling)"
        if matches_normalized
        else "i.value = w.spelling"
    )
    cur.execute(
        f"""
        WITH wanted AS (
          SELECT * FROM unnest(%s::int[], %s::text[], %s::text[])
            AS w(seed_index, type_name, spelling)
        )
        SELECT DISTINCT w.seed_index, e.entity_id::text,
               coalesce(ct.name = ANY(%s::text[]), false) AS standard_identity
        FROM wanted w
        JOIN vocab_identifier_type t ON t.name = w.type_name
        JOIN identifier_evidence i
          ON i.identifier_type_id = t.identifier_type_id AND {value_test}
        JOIN entity_identifier_lookup l ON l.identifier_id = i.identifier_id
        JOIN entity e ON e.entity_id = l.entity_id
        JOIN vocab_entity_type vt
          ON vt.entity_type_id = e.entity_type_id AND vt.name = ANY(%s::text[])
        LEFT JOIN vocab_identifier_type ct ON ct.identifier_type_id = e.canonical_identifier_type_id
        """,
        [
            [row[0] for row in wanted],
            [row[1] for row in wanted],
            [row[2] for row in wanted],
            list(STANDARD_IDENTITY_TYPES),
            list(CHEMICAL_ENTITY_TYPES),
        ],
    )
    found: dict[int, list[tuple[str, bool]]] = {}
    for index, entity_id, standard in cur.fetchall():
        found.setdefault(index, []).append((entity_id, standard))
    resolutions = []
    for index, seed in enumerate(seeds):
        rows = sorted(found.get(index, []))
        resolutions.append(
            Resolution(
                seed,
                tuple(entity_id for entity_id, _ in rows),
                bool(rows) and not any(standard for _, standard in rows),
            )
        )
    return resolutions
