"""What a seed resolves to in a published release, and the spellings it is looked up under.

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

# The namespaces a seed may name, as the release spells them (the `ns` of an identifier).
NAMESPACES = (
    "chebi",
    "kegg",
    "hmdb",
    "pubchem",
    "lipidmaps",
    "swisslipids",
    "chembl",
    "drugbank",
)

# Entity types a seed may resolve to: the published chemical families.
CHEMICAL_ENTITY_TYPES = ("chemical_entity", "small_molecule")

# An entity whose canonical identity is in one of these namespaces was resolved against a
# chemical database. Any other canonical identity (a name, a source accession) is a native
# fallback: the build had nothing better to anchor on.
STANDARD_NAMESPACES = (*NAMESPACES, "inchikey")


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
