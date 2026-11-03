"""Canonical CV labels and accepted aliases used by the build pipeline.

Source records can expose entity and identifier types as enum values,
accessions, labelled accessions, or legacy strings. These constants normalize
the small set of type families that the resolver and canonicalization phases
need to recognize consistently.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    cv_term_label_accession,
)


PROTEIN_ENTITY_TYPE = "protein"
GENE_ENTITY_TYPE = "gene"
SMALL_MOLECULE_ENTITY_TYPE = "small_molecule"
CHEMICAL_ENTITY_TYPE = "chemical_entity"
LIPID_ENTITY_TYPE = "chemical_entity"
CV_TERM_ENTITY_TYPE = "ontology_class"

PROTEIN_ENTITY_TYPE_ALIASES = (
    PROTEIN_ENTITY_TYPE,
    str(EntityTypeCv.PROTEIN),
    "protein",
    "MI:0326:Protein",
    GENE_ENTITY_TYPE,
    str(EntityTypeCv.GENE),
    "gene",
    "MI:0250:Gene",
    "protein",
    "gene",
)
CHEMICAL_ENTITY_TYPE_ALIASES = (
    CHEMICAL_ENTITY_TYPE,
    str(EntityTypeCv.CHEMICAL),
    "chemical_entity",
    "OM:0037:Chemical",
    SMALL_MOLECULE_ENTITY_TYPE,
    str(EntityTypeCv.SMALL_MOLECULE),
    "small_molecule",
    "MI:0328:Small Molecule",
    LIPID_ENTITY_TYPE,
    str(EntityTypeCv.LIPID),
    "chemical_entity",
    "OM:0011:Lipid",
    "chemical",
    "small_molecule",
    "compound",
    "drug",
)
SUPPORTED_ENTITY_TYPE_ALIASES = PROTEIN_ENTITY_TYPE_ALIASES + CHEMICAL_ENTITY_TYPE_ALIASES

CV_TERM_ID_TYPE = cv_term_label_accession(IdentifierNamespaceCv.CV_TERM_ACCESSION)


def normalize_entity_type(value: object) -> str | None:
    """Normalize enum-like entity type values to nullable text."""

    if value is None:
        return None
    if value is EntityTypeCv.SMALL_MOLECULE:
        return CHEMICAL_ENTITY_TYPE
    text = cv_term_label_accession(value)
    if text:
        if text in CHEMICAL_ENTITY_TYPE_ALIASES:
            return CHEMICAL_ENTITY_TYPE
        return text
    if hasattr(value, "value"):
        value = value.value
    text = str(value).strip()
    if text in CHEMICAL_ENTITY_TYPE_ALIASES:
        return CHEMICAL_ENTITY_TYPE
    return text or None
