"""Stable integer IDs for identifier namespaces used in resolver data."""

from __future__ import annotations


# Exact append-only main namespace registry, rendered using the vocabulary at
# omnipath-build 9f9bb709c764's pypath gitlink. No runtime parser/resolver imports.
VOCABULARY_SOURCE_COMMIT = "33f37fbaab59993d24f5c32bb4b3e7587085bcf4"
VOCABULARY_SOURCE_PATH = "pypath/internals/cv_terms/identifiers.py"
VOCABULARY_SOURCE_BLOB = "247c83e9223dcfb76ade356c8ce0f132ebc7aab3"


FALLBACK_IDENTIFIER_TYPE = "Fallback"
UNRESOLVED_ID_TYPE = "omnipath:unresolved_entity_key"
COMPLEX_MEMBER_HASH_ID_TYPE = "omnipath:complex_member_hash"
REACTION_MEMBER_HASH_ID_TYPE = "omnipath:reaction_member_hash"

IDENTIFIER_TYPE_NAMES: tuple[str, ...] = (
    "Uniprot:MI:1097",
    "Ensembl:MI:0476",
    "Entrez:MI:0477",
    "Hgnc:MI:1095",
    "Gene Name Primary:OM:0200",
    "Gene Name Synonym:OM:0201",
    "Uniprot Entry Name:OM:0221",
    "Chebi:MI:0474",
    "Hmdb:OM:0004",
    "Lipidmaps:OM:0003",
    "Swisslipids:OM:0009",
    "Pubchem Compound:OM:0002",
    "Standard Inchi Key:MI:1101",
    "Cv Term Accession:OM:0204",
    "Fallback",
    "Name:OM:0202",
    "Chembl Compound:MI:0967",
    "omnipath:unresolved_entity_key",
    "omnipath:complex_member_hash",
    "omnipath:reaction_member_hash",
    "Ramp Id:OM:0132",
    "Refmet:OM:0137",
    "Kegg Compound:MI:2012",
    "Cas:MI:2011",
    "Mirbase Precursor:OM:0127",
    "Mirbase Mature:OM:0128",
    "Lipid Name:OM:0209",
    "Drugbank:MI:2002",
    "Reactome Stable Id:OM:0130",
)

IDENTIFIER_TYPE_IDS: dict[str, int] = {
    name: index for index, name in enumerate(IDENTIFIER_TYPE_NAMES, start=1)
}


def identifier_type_id(name: str) -> int:
    """Return the stable integer ID for an identifier namespace."""

    try:
        return IDENTIFIER_TYPE_IDS[name]
    except KeyError as error:
        raise ValueError(f"Unknown resolver identifier type: {name!r}") from error


def identifier_type_rows(
    names: set[str] | None = None,
) -> list[dict[str, object]]:
    """Return resolver identifier type rows for all or selected namespaces."""

    selected = set(IDENTIFIER_TYPE_NAMES if names is None else names)
    return [
        {"identifier_type_id": identifier_type_id(name), "name": name}
        for name in IDENTIFIER_TYPE_NAMES
        if name in selected
    ]


# Vocabulary-only metadata previously imported from the old ingest module.
# No mapping/resolution code is imported or executed by the Parquet backend.
STANDARD_INCHI_KEY_TYPE = "Standard Inchi Key:MI:1101"
RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE = {
    "pubchem": "Pubchem Compound:OM:0002",
    "chembl": "Chembl Compound:MI:0967",
    "chebi": "Chebi:MI:0474",
    "hmdb": "Hmdb:OM:0004",
    "lipidmaps": "Lipidmaps:OM:0003",
    "swisslipids": "Swisslipids:OM:0009",
    "kegg": "Kegg Compound:MI:2012",
    "drugbank": "Drugbank:MI:2002",
    "reactome": "Reactome Stable Id:OM:0130",
}
