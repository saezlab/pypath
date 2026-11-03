"""Stable integer IDs for identifier namespaces used in resolver data."""

from __future__ import annotations

import polars as pl

from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)

IDENTIFIER_TYPE_SCHEMA: dict[str, pl.DataType] = {
    'identifier_type_id': pl.UInt32,
    'name': pl.Utf8,
}

FALLBACK_IDENTIFIER_TYPE = 'Fallback'
UNRESOLVED_ID_TYPE = 'omnipath:unresolved_entity_key'
COMPLEX_MEMBER_HASH_ID_TYPE = 'omnipath:complex_member_hash'
REACTION_MEMBER_HASH_ID_TYPE = 'omnipath:reaction_member_hash'

IDENTIFIER_TYPE_NAMES: tuple[str, ...] = (
    cv_term_label_accession(IdentifierNamespaceCv.UNIPROT),
    cv_term_label_accession(IdentifierNamespaceCv.ENSEMBL),
    cv_term_label_accession(IdentifierNamespaceCv.ENTREZ),
    cv_term_label_accession(IdentifierNamespaceCv.HGNC),
    cv_term_label_accession(IdentifierNamespaceCv.GENE_NAME_PRIMARY),
    cv_term_label_accession(IdentifierNamespaceCv.GENE_NAME_SYNONYM),
    cv_term_label_accession(IdentifierNamespaceCv.UNIPROT_ENTRY_NAME),
    cv_term_label_accession(IdentifierNamespaceCv.CHEBI),
    cv_term_label_accession(IdentifierNamespaceCv.HMDB),
    cv_term_label_accession(IdentifierNamespaceCv.LIPIDMAPS),
    cv_term_label_accession(IdentifierNamespaceCv.SWISSLIPIDS),
    cv_term_label_accession(IdentifierNamespaceCv.PUBCHEM_COMPOUND),
    cv_term_label_accession(IdentifierNamespaceCv.STANDARD_INCHI_KEY),
    cv_term_label_accession(IdentifierNamespaceCv.CV_TERM_ACCESSION),
    FALLBACK_IDENTIFIER_TYPE,
    cv_term_label_accession(IdentifierNamespaceCv.NAME),
    cv_term_label_accession(IdentifierNamespaceCv.CHEMBL_COMPOUND),
    UNRESOLVED_ID_TYPE,
    COMPLEX_MEMBER_HASH_ID_TYPE,
    REACTION_MEMBER_HASH_ID_TYPE,
    cv_term_label_accession(IdentifierNamespaceCv.RAMP_ID),
    cv_term_label_accession(IdentifierNamespaceCv.REFMET),
    cv_term_label_accession(IdentifierNamespaceCv.KEGG_COMPOUND),
    cv_term_label_accession(IdentifierNamespaceCv.CAS),
    # miRNA accessions — appended to keep existing IDs stable
    cv_term_label_accession(IdentifierNamespaceCv.MIRBASE_PRECURSOR),
    cv_term_label_accession(IdentifierNamespaceCv.MIRBASE_MATURE),
    # Standardized lipid nomenclature name (spec 011 T116) — appended to keep
    # existing IDs stable, same convention as the miRNA accessions above.
    cv_term_label_accession(IdentifierNamespaceCv.LIPID_NAME),
    # DrugBank and Reactome, widened into the resolver_chemical structure
    # route (spec 011 T092) — appended here too, the same append-only
    # convention, since every resolver setup call validates its namespace
    # slugs against this registry (found live: T092 registered both in
    # RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE, duckdb_load.py, but not
    # here, so every dataset's resolver setup failed identically with
    # "Unknown resolver identifier type: 'Drugbank:MI:2002'" during the
    # cycle's final reload).
    cv_term_label_accession(IdentifierNamespaceCv.DRUGBANK),
    cv_term_label_accession(IdentifierNamespaceCv.REACTOME_STABLE_ID),
)

IDENTIFIER_TYPE_IDS: dict[str, int] = {
    name: index for index, name in enumerate(IDENTIFIER_TYPE_NAMES, start=1)
}


def identifier_type_id(name: str) -> int:
    """Return the stable integer ID for an identifier namespace."""

    try:
        return IDENTIFIER_TYPE_IDS[name]
    except KeyError as error:
        raise ValueError(
            f'Unknown resolver identifier type: {name!r}'
        ) from error


def identifier_type_rows(
    names: set[str] | None = None,
) -> list[dict[str, object]]:
    """Return resolver identifier type rows for all or selected namespaces."""

    selected = set(IDENTIFIER_TYPE_NAMES if names is None else names)
    return [
        {'identifier_type_id': identifier_type_id(name), 'name': name}
        for name in IDENTIFIER_TYPE_NAMES
        if name in selected
    ]
