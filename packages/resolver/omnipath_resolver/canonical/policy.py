"""Declarative canonicalization policies, one per entity class.

A policy says which reference **library** an entity class is matched against
(``gene_protein`` or ``chemical``; ``None`` for classes that are never matched),
which identifier namespaces may **vote** during matching, which are
low-confidence **symbol** votes that need a taxon, and how to pick a **label**.

The matching algorithm itself lives in :mod:`.match` and is identical for all
classes; only these parameters differ.  There is no fallback cascade: an
observation either matches a library node or keeps its own primary identifier.
See ``canonical/README.md``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable
from biolink_model.datamodel import model
from omnipath_core.biolink import entity_type as biolink_entity_type, schema as biolink_schema

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

GENE_PROTEIN = "gene_protein"
CHEMICAL = "chemical"
REACTION = "reaction"
REACTION_NAMESPACES = frozenset(
    {
        "rhea",
        "kegg_reaction",
        "reactome",
        "metacyc_reaction",
        "ecocyc_reaction",
        "macie",
        "bigg_reaction",
        "vmh_reaction",
        "seed_reaction",
        "sabiork_reaction",
        "metanetx_reaction",
    }
)


@dataclass(frozen=True)
class EntityPolicy:
    entity_class: str
    entity_types: frozenset[str]
    library: str | None
    match_namespaces: frozenset[str]
    symbol_namespaces: frozenset[str]
    rewrites: tuple[tuple[str, str], ...]
    label_priority: tuple[str, ...]
    label_filters: tuple[Callable[[str], bool], ...]
    alias_namespaces: frozenset[str]

    @property
    def matched(self) -> bool:
        return self.library is not None


_NONE: frozenset[str] = frozenset()


RNA_ENTITY_TYPES = frozenset(
    biolink_entity_type(name)
    for name in biolink_schema().class_descendants("transcript", reflexive=True)
)
PROTEIN_ENTITY_TYPES = frozenset(
    biolink_entity_type(name)
    for name in biolink_schema().class_descendants("protein", reflexive=True)
)

GENE_PROTEIN_POLICY = EntityPolicy(
    entity_class=GENE_PROTEIN,
    entity_types=RNA_ENTITY_TYPES | PROTEIN_ENTITY_TYPES | {biolink_entity_type(model.Gene)},
    library=GENE_PROTEIN,
    match_namespaces=frozenset(
        {
            "uniprot",
            "uniprot-sec",
            "uniprot_entry",
            "entrez",
            "ensg",
            "enst",
            "ensp",
            "hgnc",
            "refseq",
            "refseq_protein",
            "genbank",
            "ramp_gene",
            "kegg_gene",
        }
    ),
    symbol_namespaces=frozenset({"genesymbol", "genesymbol-syn"}),
    rewrites=(),
    label_priority=("genesymbol", "name"),
    label_filters=(lambda s: len(s) <= 120,),
    alias_namespaces=frozenset(
        {
            "uniprot",
            "uniprot-sec",
            "uniprot_entry",
            "entrez",
            "ensg",
            "enst",
            "ensp",
            "hgnc",
            "refseq",
            "refseq_protein",
            "genbank",
            "ramp_gene",
            "kegg_gene",
            "genesymbol",
            "genesymbol-syn",
        }
    ),
)

CHEMICAL_POLICY = EntityPolicy(
    entity_class=CHEMICAL,
    entity_types=frozenset(
        biolink_entity_type(cls) for cls in (model.ChemicalEntity, model.SmallMolecule)
    ),
    library=CHEMICAL,
    match_namespaces=frozenset(
        {
            "inchikey",
            "inchi",
            "chebi",
            "chembl",
            "pubchem",
            "hmdb",
            "kegg",
            "cas",
            "lipidmaps",
            "swisslipids",
            "drugbank",
            "bigg",
            "metanetx",
            "goslin",
            "refmet",
            "ramp",
        }
    ),
    symbol_namespaces=_NONE,
    rewrites=(),
    label_priority=("name", "chebi", "pubchem"),
    label_filters=(lambda s: not INCHIKEY_RE.match(s) and len(s) <= 120,),
    alias_namespaces=frozenset(
        {
            "inchikey",
            "chebi",
            "chembl",
            "pubchem",
            "hmdb",
            "kegg",
            "cas",
            "lipidmaps",
            "swisslipids",
            "drugbank",
            "bigg",
            "metanetx",
            "goslin",
            "refmet",
            "ramp",
            "name",
        }
    ),
)

REACTION_POLICY = EntityPolicy(
    entity_class=REACTION,
    entity_types=frozenset({biolink_entity_type(model.MolecularActivity)}),
    library=REACTION,
    match_namespaces=REACTION_NAMESPACES,
    symbol_namespaces=_NONE,
    rewrites=(),
    label_priority=("name",),
    label_filters=(lambda s: len(s) <= 300,),
    alias_namespaces=REACTION_NAMESPACES | {"name"},
)

CV_TERM_POLICY = EntityPolicy(
    entity_class="cv_term",
    entity_types=frozenset({biolink_entity_type(model.OntologyClass)}),
    library=None,
    match_namespaces=_NONE,
    symbol_namespaces=_NONE,
    rewrites=((r"^KW-\d+$", "uniprot_keyword"),),
    label_priority=("uniprot_keyword", "name"),
    label_filters=(),
    alias_namespaces=_NONE,
)

COMPLEX_POLICY = EntityPolicy(
    entity_class="complex",
    entity_types=frozenset(
        {biolink_entity_type(cls) for cls in (model.MacromolecularComplex, model.ProteinFamily)}
    ),
    library=None,
    match_namespaces=_NONE,
    symbol_namespaces=_NONE,
    rewrites=(),
    label_priority=("name", "genesymbol"),
    label_filters=(),
    alias_namespaces=_NONE,
)

GENERIC_POLICY = EntityPolicy(
    entity_class="generic",
    entity_types=frozenset(),
    library=None,
    match_namespaces=_NONE,
    symbol_namespaces=_NONE,
    rewrites=(),
    label_priority=("name", "genesymbol"),
    label_filters=(),
    alias_namespaces=_NONE,
)

POLICIES: dict[str, EntityPolicy] = {}
for _policy in (GENE_PROTEIN_POLICY, CHEMICAL_POLICY, REACTION_POLICY, CV_TERM_POLICY, COMPLEX_POLICY):
    for _etype in _policy.entity_types:
        POLICIES[_etype] = _policy


def get_policy(entity_type: str) -> EntityPolicy:
    """Policy for a raw entity type."""
    norm = biolink_entity_type(entity_type)
    return POLICIES.get(norm, GENERIC_POLICY)
