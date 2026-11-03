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


PROTEIN_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.PROTEIN)
GENE_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.GENE)
SMALL_MOLECULE_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.SMALL_MOLECULE)
CHEMICAL_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.CHEMICAL)
LIPID_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.LIPID)
CV_TERM_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.CV_TERM)

PROTEIN_ENTITY_TYPE_ALIASES = (
    PROTEIN_ENTITY_TYPE,
    str(EntityTypeCv.PROTEIN),
    'Protein:MI:0326',
    'MI:0326:Protein',
    GENE_ENTITY_TYPE,
    str(EntityTypeCv.GENE),
    'Gene:MI:0250',
    'MI:0250:Gene',
    'protein',
    'gene',
)
CHEMICAL_ENTITY_TYPE_ALIASES = (
    CHEMICAL_ENTITY_TYPE,
    str(EntityTypeCv.CHEMICAL),
    'Chemical:OM:0037',
    'OM:0037:Chemical',
    SMALL_MOLECULE_ENTITY_TYPE,
    str(EntityTypeCv.SMALL_MOLECULE),
    'Small Molecule:MI:0328',
    'MI:0328:Small Molecule',
    LIPID_ENTITY_TYPE,
    str(EntityTypeCv.LIPID),
    'Lipid:OM:0011',
    'OM:0011:Lipid',
    'chemical',
    'small_molecule',
    'compound',
    'drug',
)
SUPPORTED_ENTITY_TYPE_ALIASES = (
    PROTEIN_ENTITY_TYPE_ALIASES + CHEMICAL_ENTITY_TYPE_ALIASES
)

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
    if hasattr(value, 'value'):
        value = value.value
    text = str(value).strip()
    if text in CHEMICAL_ENTITY_TYPE_ALIASES:
        return CHEMICAL_ENTITY_TYPE
    return text or None
