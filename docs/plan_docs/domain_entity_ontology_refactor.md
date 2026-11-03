# Domain Entities and Ontology Relations Refactor

## Summary

OmniPath integrates biological knowledge from many sources into a canonical
entity and relation graph. Some source namespaces are also ontologies. ChEBI,
Reactome pathways, WikiPathways, KEGG pathways, GO, HPO, MONDO, and similar
resources provide identifiers, labels, hierarchy, and sometimes domain records.

The current model partly treats ontology-backed records as separate `CV_TERM`
entities and partly treats them as domain entities such as small molecules,
pathways, or reactions. This creates ambiguity for mixed resources such as
ChEBI, where concrete structures and abstract chemical classes share the same
namespace.

The target model is:

- Domain entities keep their domain type, for example `Chemical` or `Pathway`.
- Source identifiers keep their source namespace, for example ChEBI or Reactome
  stable ID.
- Ontology hierarchy is stored separately in `entity_ontology_relation`.
- Evidence-level ontology/pathway annotations stay simple and use
  `associated_with`.
- The dedicated `ontology_terms` table can be removed once ontology metadata is
  represented as entity identifiers, entity annotations, and ontology
  relations.

## Current Situation

### Entity Canonicalization

The DuckDB load pipeline canonicalizes source evidence into entities.

- Proteins and genes resolve through protein resolver tables, usually to UniProt.
- Chemicals resolve through chemical resolver tables, currently toward standard
  InChIKey.
- `CV_TERM` evidence resolves directly by `CV_TERM_ACCESSION`.
- Unresolved evidence falls back to a synthetic unresolved identifier hash.

The canonical entity key includes:

```text
entity_type
taxonomy_id
canonical_identifier_type
canonical_identifier
```

This means entity type is semantically important. A `CV_TERM` identified by
`CHEBI:46217` is not the same entity as a `Chemical` identified by
`CHEBI:46217`.

### Ontology Terms

Ontology datasets are currently projected into a dedicated `ontology_terms`
table. The load path also creates synthetic `CV_TERM` entities for ontology
terms. This creates a separate ontology island.

This works for pure vocabularies, but it is awkward for resources where ontology
terms are also domain objects.

### ChEBI

ChEBI mixes:

- concrete chemical structures with standard InChI/InChIKey
- residues and fragments without standard InChIKey
- chemical classes and hierarchy terms

Current chemical resolver behavior works for structural terms, for example a
ChEBI ID resolves to a standard InChIKey. It does not resolve non-structural
ChEBI participants such as `CHEBI:46217` (`L-alanine residue`), because no
standard InChIKey exists.

Current ChEBI ontology export also filters structural leaf terms out of the
ontology export. That is useful for avoiding a huge ontology browser full of
concrete structures, but terms that cannot be structure-canonicalized still need
to be represented.

### Reactome Pathways

Reactome pathway data has two uses:

- domain pathway entities identified by Reactome stable IDs
- pathway hierarchy suitable for browsing and ancestor traversal

Current pathway hierarchy is exported as ontology data, while reactions and
participants can carry pathway accessions as annotations. This makes pathway
hierarchy separate from the domain pathway entities used elsewhere in the graph.

## Target Model

### Domain Entity Types

Use domain entity types for ontology-backed domain objects:

```text
Chemical
Pathway
Protein
Reaction
Complex
CV_TERM
```

`CV_TERM` remains available for vocabulary concepts that do not have a better
domain entity type.

For ChEBI, `Chemical` means chemical entity, residue, fragment, or chemical
class. It does not imply a concrete structure.

### Identifier Types

Keep source-specific identifier namespaces:

```text
ChEBI
Reactome Stable ID
WikiPathways
KEGG Pathway
GO
HPO
MONDO
Standard InChIKey
```

Do not replace source namespaces with `CV_TERM_ACCESSION` simply because a
namespace has ontology metadata.

Examples:

```text
CHEBI:16449 alanine
entity_type = Chemical
canonical_identifier_type = Standard InChIKey
canonical_identifier = QNAYBMKLOCPYGJ-UHFFFAOYSA-N

CHEBI:46217 L-alanine residue
entity_type = Chemical
canonical_identifier_type = ChEBI
canonical_identifier = 46217

R-HSA-123456
entity_type = Pathway
canonical_identifier_type = Reactome Stable ID
canonical_identifier = R-HSA-123456
```

### Evidence-Level Annotations

Evidence-level ontology/pathway annotations should stay simple:

```text
source entity --associated_with--> annotation target entity
```

Do not infer `part_of` for evidence-level annotations. `part_of`, `is_a`, and
similar hierarchy predicates belong to ontology relations, not source evidence
annotations.

### Ontology Relations

Add a separate ontology hierarchy table:

```text
entity_ontology_relation
  source_id
  subject_entity_id
  predicate_id
  object_entity_id
  ontology_id
```

`relation_type` is not needed if `predicate_id` captures the relation semantics.
UI grouping can be derived from predicate metadata if required later.

This table is not part of the default relation graph and should not appear in
default predicate filters.

## Refactoring Plan

### Phase 1: Stabilize ChEBI Parsing Boundaries

Goal: make ChEBI parsing consistent before changing database schema.

Tasks:

- Preserve ChEBI `alt_id` values in parsed records.
- Emit molecule rows only for terms with a standard InChIKey.
- Keep non-structural ChEBI terms, such as residues without InChIKey, available
  for ontology/domain processing.
- Keep concrete structure leaves out of the ontology browser by default.

Acceptance checks:

- `CHEBI:16449` appears in ChEBI molecule rows.
- `CHEBI:46217` does not appear in molecule rows.
- `CHEBI:46217` remains available in the ChEBI ontology/domain term stream.
- `CHEBI:2539` can be recognized as an alias of `CHEBI:16449`.

### Phase 2: Introduce Domain-Backed Ontology Metadata

Goal: stop requiring ontology-backed domain records to become `CV_TERM`
entities.

Tasks:

- Treat database-load ontology metadata as normal entity output; the old
  separate ontology dataset abstraction has been removed.
- Add optional `ontology_relations` to silver `Entity` records and entity
  builders.
- Remove OBO export infrastructure from the active input/load path.
- Define which ontology-like entity streams are domain-backed:
  - ChEBI -> `Chemical`
  - Reactome pathway ontology -> `Pathway`
  - WikiPathways pathway ontology -> `Pathway`
  - KEGG pathway ontology -> `Pathway`
- Keep pure vocabularies as `CV_TERM` unless a better domain type exists.
- Add mapping metadata to entity streams so ontology-backed entities can
  declare:
  - entity type
  - identifier type
  - canonical identifier normalization

Acceptance checks:

- ChEBI ontology records can materialize as `Chemical` entities.
- Reactome pathway ontology records can materialize as `Pathway` entities.
- GO/HPO/MONDO behavior remains explicit and does not change accidentally.

### Phase 3: Extend Resolver Semantics

Goal: allow resolvers to canonicalize to non-structural source identifiers, not
only structure identifiers.

Tasks:

- Allow chemical resolver rows to have mixed canonical identifier types.
- For ChEBI structural terms:

```text
key_identifier_type = ChEBI
key_value = 16449
canonical_identifier_type = Standard InChIKey
canonical_identifier = QNAYBMKLOCPYGJ-UHFFFAOYSA-N
```

- For ChEBI non-structural terms:

```text
key_identifier_type = ChEBI
key_value = 46217
canonical_identifier_type = ChEBI
canonical_identifier = 46217
```

- Keep ambiguity detection based on `(key identifier type, key value,
  canonical identifier type)`.
- Ensure resolver lookup views do not assume chemical canonical identifiers are
  always standard InChIKey.

Acceptance checks:

- `CHEBI:46217` resolves as a `Chemical` with canonical identifier type ChEBI.
- Structural ChEBI terms still collapse by standard InChIKey across HMDB,
  ChEMBL, PubChem, etc.
- Existing chemical cross-resource resolution does not regress.

### Phase 4: Add Entity Ontology Relation Projection

Goal: store ontology hierarchy separately from evidence relations.

Tasks:

- Add `entity_ontology_relation` table to the PostgreSQL schema.
- Add a DuckDB staging table for ontology relations.
- Convert ontology `is_a` and supported relationship edges to
  `entity_ontology_relation`.
- Resolve ontology relation endpoints through the same canonical entity key
  model as normal entities.
- Keep `entity_ontology_relation` out of `relation`, `relation_evidence`, and
  `relation_evidence_relation`.

Acceptance checks:

- ChEBI hierarchy edges load into `entity_ontology_relation`.
- Reactome pathway hierarchy edges load into `entity_ontology_relation`.
- Default graph relation counts and predicate filters do not include ontology
  hierarchy edges.

### Phase 5: Simplify Evidence-Level Ontology Annotation Projection

Goal: make source annotations predictable and avoid over-interpreting
annotation semantics.

Tasks:

- Project ontology/pathway annotations as `associated_with`.
- Do not emit evidence-level `part_of` for pathway annotations.
- Resolve annotation targets to domain entities when the namespace is
  domain-backed.
- Resolve annotation targets to `CV_TERM` only for pure vocabulary namespaces.

Acceptance checks:

- Reaction/pathway membership annotations use `associated_with`.
- Protein/chemical/pathway annotations do not introduce hierarchy predicates.
- UI can distinguish evidence-level associations from ontology hierarchy.

### Phase 6: Replace `ontology_terms`

Goal: remove the dedicated ontology island once domain-backed ontology metadata
is represented through normal entity tables and ontology relation tables.

Tasks:

- Move ontology labels, definitions, synonyms, obsolete flags, and ontology IDs
  into entity annotations and identifiers.
- Keep a compatibility view named `ontology_terms` during transition if API/UI
  code depends on it.
- Update UI queries to use:
  - `entity`
  - `entity_identifier`
  - `entity_annotation`
  - `entity_ontology_relation`
- Drop physical `ontology_terms` after compatibility consumers are migrated.

Acceptance checks:

- Ontology browser can render ChEBI and Reactome hierarchy without
  `ontology_terms`.
- Search can find ontology-backed entities by name, synonym, and source
  identifier.
- Existing API endpoints either continue through a compatibility view or are
  intentionally migrated.

## Open Decisions

- Whether to rename `Small Molecule` to `Chemical`, or introduce `Chemical` as a
  broader parent term while keeping `Small Molecule` as a subtype.
- Which vocabularies are pure `CV_TERM` and which are domain-backed.
- Whether ChEBI concrete structural leaves should be hidden from ontology
  browser queries by default or excluded from ontology relation materialization.
- How much predicate normalization is needed for ontology relationships beyond
  `is_a`, `part_of`, and `has_part`.
- Whether ontology relation endpoints should store only entity IDs or also
  denormalized identifier fields for easier debugging.

## Non-Goals

- Do not merge domain entities with `CV_TERM` entities solely because they share
  an accession.
- Do not put ontology hierarchy edges into the default relation graph.
- Do not infer precise hierarchy predicates from evidence-level annotations.
- Do not remove source-specific identifier namespaces in favor of generic
  `CV_TERM_ACCESSION`.
