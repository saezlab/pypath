# Published Parquet to main PostgreSQL mapping

Baseline: current main `9f9bb709c764e2721a0bcd881cf7ef0eb22bed9a`.
The loader preserves published Parquet files, entity keys, statement keys and
upstream resolution. This document describes the aligned PostgreSQL projection; it
does not claim that unavailable original resolver diagnostics were recovered.

`aligned_projection.prepare_aligned_release(duckdb_connection, pinned_release,
dimension_rows=..., retain_published_provenance=False)` stages typed input, validates canonical conflicts, prepares
main dictionaries and returns ordered `CopyQuery(table, columns, query)` records.
The caller owns PostgreSQL DDL, source partitions, transactions and streaming
COPY. Existing dimension IDs are reused. New names receive deterministic IDs
after the largest existing ID; sequence advancement belongs to the loader.
Large flattenings and final reductions are materialized once in caller-owned
on-disk DuckDB. Python receives only dimension rows and aggregate diagnostics.

## Main base columns

| PostgreSQL table/columns | Published input and treatment |
| --- | --- |
| `data_source(source_id,name)` | Resource names and declared relation-evidence source names. The publishing resource is separate from the evidence's declared source. Null declared source uses the publishing resource for ownership, with the original null retained. |
| `dataset(dataset_id,source_id,name)` | Declared evidence dataset. A missing dataset uses the explicit `omnipath:published_statement` aggregate namespace. Entity aggregate evidence uses `omnipath:published_entities`. Original nulls and strings remain recoverable. |
| `vocab_entity_type` | Published type strings, including Biolink values. No protein-to-gene recanonicalization. |
| `vocab_identifier_type` | Explicit namespace spelling adapters to main's existing namespace labels, including KEGG reaction `Kegg Reaction:MI:2013` and SMILES `Smiles:MI:0239`. Unrecognized namespaces remain their published spelling. No identifier resolution is repeated. |
| `vocab_relation_predicate` | Published predicate names; Biolink remains the vocabulary. |
| `vocab_relation_category` | Published category; a null is represented by the explicit `omnipath:unspecified` category and remains null in the statement companion. |
| `vocab_annotation_scope` | Main's `relation`, `subject`, `object` IDs remain 1/2/3. Published `evidence` scope maps to main relation scope because the owner is already a specific evidence record; its original scope remains in occurrence provenance. Other explicit scopes are registered, not discarded. |
| `entity.entity_id` | Deterministic content UUID of the complete published entity key. Keys are not truncated. One UUID per distinct published identity across the release. |
| `entity.entity_type_id` | Dimension reference for the published type. Conflicting type, namespace or identifier for the same key is rejected. |
| `entity.taxonomy_id` | A single distinct known published taxon is retained. More than one known taxon becomes canonical NULL, with every source occurrence's taxon retained in `entity_evidence` and `parquet_entity`. This represents lack of a release-wide taxon consensus; it does not alter identity. |
| `entity.canonical_identifier_type_id`, `canonical_identifier` | Published namespace/identifier, with spelling-only namespace adaptation. For colliding source-scoped name identities, use main's existing `omnipath:unresolved_entity_key` namespace and the full published key as the storage identity. The original name remains a searchable identifier and in the entity crosswalk. Distinct entities are never silently merged. Other main natural-key collisions fail before COPY. |
| `entity.resolution_status_id`, `resolution_mechanism` | New explicit status 5 `published` and mechanism `published_parquet`. These describe the import boundary; they do not assert old `resolved`, `matched`, `ambiguous`, or `unresolved` decisions. |
| `entity.created_at` and other optional identity/classification fields | Main defaults and downstream derivations. No molecular state, gene anchoring, resolver explanation or resolution detail is fabricated. |
| `identifier_evidence.identifier_id` | Content UUID of the main namespace spelling plus identifier value; globally deduplicated by that pair. Main's explicit Standard InChI exclusion applies to generic dictionary and link rows; full original values remain in pinned Parquets and optional inspection copies with a NULL generic dictionary reference. No general length cutoff or truncation applies. Exact accepted published values remain. Null namespace/value occurrences have no invented dictionary row and remain in occurrence provenance. |
| `identifier_evidence.identifier_type_id`, `value` | Typed namespace reference and unchanged identifier string, plus safe bare lookup aliases for known outer CURIE prefixes and well-formed suffixes (for example `chebi` / `CHEBI:15377` also has `15377`). Intrinsic prefixes such as `CHEMBL`, `HMDB` and `SLM:` stay in the ID. Unknown namespaces and arbitrary colon-containing names are never stripped. Canonical identifiers are included even when an input identifier list is empty. |
| `identifier_evidence.value_normalized` | NULL at import: the published serving contract does not retain the old normalized alias field. Existing downstream normalization can populate supported search representations without changing entity identity. |
| `entity_identifier(source_id,entity_id,identifier_id)` | Deduplicated resource-owned canonical entity/alias links. Includes the published canonical identifier, safe bare aliases and the explicit storage canonical identifier for scoped-name fallbacks. Original duplicate occurrences are retained separately; generated lookup aliases do not fabricate published occurrences or alter entity identity. Eligible source-owned ontology-class/pathway IDs also have explicit controlled-vocabulary lookup aliases, retaining their exact strings. Genuine published scalar ontology labels may supply Name lookups only when no name/synonym exists on that entity and the label is not an identifier fallback. These aliases carry no invented occurrence or authority claim. |
| `entity_evidence` | One explicit aggregate per publishing resource/version/entity. UUID identifies that aggregate; source is the publishing resource; dataset is the aggregate dataset; row ID is a deterministic bigint aggregate surrogate; role is parent; original parent/member trees are unavailable. Type and taxonomy reflect the specific resource occurrence. |
| `entity_evidence_identifier` | Deduplicated links between the aggregate entity evidence and each valid published identifier. Computed bare lookup aliases receive no evidence link, preventing fabricated source occurrences or authority/reference roles. The original occurrence ordinal, canonical flag and declared source remain separate. |
| `entity_evidence_resolution` | Links the aggregate evidence to its published canonical entity using status 5. `reason_id` and `molecular_type_id` remain NULL when absent. Timestamp uses main's default. No downstream resolver is called. |
| `annotation.annotation_key` | Content UUID of term, original value and every quantity component. Distinct prefixes, comparators, source fields, units and all-null-versus-absent quantity structs remain distinct. |
| `annotation.term`, `value`, `unit` | Published term. Numeric quantities use the numeric value as main's scalar text `value` and the explicit unit as `unit`; exact original value is retained in `annotation_quantity.published_value`. Scalar/symbolic values are unchanged. A null-term occurrence is retained in provenance instead of inventing a main non-null term. |
| `entity_evidence_annotation` | Deduplicated annotation links on aggregate entity evidence. Main's dictionary/link primary keys are preserved; duplicate published occurrences and their attribution remain in provenance. |
| `relation.relation_id` | Main's canonical graph triple: subject UUID, predicate, object UUID. Multiple qualified published statements can map to one triple; their identities and qualifiers remain distinct. Ontology axioms do not enter this graph table. |
| `relation.subject_entity_id`, `predicate_id`, `object_entity_id` | Published endpoints and predicate via typed references. Missing endpoints fail validation. Direction and sign are not inferred from endpoint order. |
| `relation.relation_category_id` | The one agreed category for a graph triple, or NULL when qualified claims disagree. Per-claim category is always retained. |
| `relation_evidence.relation_evidence_id` | Deterministic resource/version/statement/ordinal UUID. Duplicate evidence positions remain distinct. A statement with no evidence receives explicitly synthetic aggregate claim evidence, marked as such in provenance. |
| `relation_evidence.source_id`, `dataset_id`, `row_id` | Source-owned declared context. Representable numeric row strings become bigint while preserving original text, including leading zeroes. Other strings use a deterministic bigint surrogate and retain exact text. Source/dataset row representation collisions are checked before COPY. |
| `relation_evidence` endpoint columns | Canonical endpoint UUID columns are populated; occurrence endpoint columns are NULL. This uses main's existing exactly-one-endpoint representation and avoids presenting synthetic aggregate entities as original subject/object occurrences. |
| `relation_evidence.predicate_id`, `relation_category_id` | Typed published claim references. |
| `relation_evidence_relation` | Links each graph evidence assertion to the corresponding main graph triple. Ontology assertions retain evidence and ontology links without an invented graph relation. |
| `relation_evidence_annotation` | Evidence annotations link to their exact occurrence. The writer's statement annotation array unions occurrence annotations; values already present in occurrence evidence are not broadcast into other events. True unmatched statement attributes respect declared source/dataset ownership. Scope is mapped as above. Original attribution and occurrence order remain separate. |
| `entity_ontology_relation` | Source-owned published `ontology` statements (and legacy `ontology_axiom`), including non-hierarchy predicates. Known input providers have an explicit ontology family independent of post-resolution canonical namespace. Graph statements keep their original namespace gates; unsupported ontology providers retain a declared published-scope fallback and capability count. |
| `ontology_terms` | Exact main legacy base table remains empty on a fresh build. Main's current projector has no writer for `ontology_terms_raw`; its coordinator records zero. Do not add a separate published-term enrichment that guesses original term identity from resolved aliases. |
| `entity_ontology_term` | Unchanged main derivation from source-owned ontology edges, relational identifiers and annotations. It selects term IDs, labels, definitions, aliases and synonyms with main's existing priorities and builds the existing search fields/indexes. |
| `entity_annotation_relation` | No invented resolved annotation target. Published entity-valued associations already represented as statements remain graph relations. Additional old annotation-resolution bridges require a demonstrably unambiguous retained target; unavailable links are reported instead of guessed. |

## Required molecular context

Serving schema 4 adds three required companions alongside the normalized main
tables. They are copied with or without `retain_published_provenance`:

| PostgreSQL table | Published input and treatment |
| --- | --- |
| `entity_reference_context` | Resource/version/entity key and main entity UUID, exact `reference_entity_key` and ordered `gene_reference_keys`. A virtual gene reference need not have an entity row. |
| `statement_reference_context` | Resource/version/relation key and both endpoint reference keys. |
| `molecular_evidence_context` | One row per ordered relation or standalone entity evidence occurrence. Nested `occurrence_json` keeps the paired forms, specific identifiers, features, coordinates and annotations together; relation occurrences link to the existing relation-evidence UUID. Nulls, empty lists and repetition survive. |

The loader performs no additional resolution, product selection or state
combination expansion. Source types and normalized graph triples retain their
existing identity rules. COSMOS keeps canonical GeneID catalysts as GeneID
instead of selecting one of their protein aliases.

## Quantity metadata and optional inspection copies

Default loads create the three molecular companions and `annotation_quantity` alongside the fourteen main COPY tables. Quantity metadata is required by main's measurement readers. The five `parquet_*` tables below are optional inspection copies: pass `retain_published_provenance=True` to the projector/loader and `companion_ddl`, or use `--retain-published-provenance` for a fresh CLI load. Both modes produce identical normalized scientific and required context tables.

References above to original occurrence provenance mean the pinned immutable Parquets by default, plus these PostgreSQL copies when enabled. Main evidence ownership, taxonomy, qualifiers and measurements remain relational; original keys, raw spelling, repetition, order and list-presence flags are inspection details. No downstream scientific reader consumes the five optional tables.

Available companion definitions:

- `parquet_entity`: resource/version/published key → compact entity and aggregate
  evidence UUIDs, original namespace/identifier/type/taxon, label, hierarchy
  counts, and identifier/annotation list-presence flags.
- `parquet_statement`: resource/version/published statement key → main graph
  relation where applicable, endpoints, predicate, original endpoint labels and
  types, taxon, directedness, sign, category, interaction class, evidence count,
  ordered source array including null entries, and list-presence flags.
- `parquet_evidence`: original ordinal, source-owned evidence UUID, declared
  source/dataset, exact original row string, upstream ID, annotation list-presence
  and explicit synthetic flag.
- `parquet_identifier_occurrence`: duplicate/order-preserving published alias
  occurrences, exact namespace/value, canonical flag and declared source.
- `parquet_annotation_occurrence`: owner kind/key, evidence ordinal, annotation
  ordinal/key, original source/dataset/scope. Null-term values remain in its
  nullable `untyped_value` field.
- `annotation_quantity`: typed numeric value, unit, prefix, binary relation,
  original source field, comparator and exact published value, keyed to the
  dictionary annotation identity. A present all-null quantity has a row; an
  absent quantity does not.

The molecular companion stores nested published evidence occurrences. Complete
entity/relation `record_json` and raw source bodies remain outside this aligned
projection. Null and empty lists remain distinguishable in molecular context,
the unchanged original Parquets and, when enabled, their inspection copies.

Loaded counts describe only copied tables. Compatibility metadata separately
records published input counts and whether provenance lives in pinned Parquets
alone or also in PostgreSQL. PostgreSQL-only finish preserves older twenty-table
checkpoints and their stored counts; it does not remove existing inspection copies.

## Audited compatibility boundary

The full existing release scalar audit found 2,188 shared entity keys with
different published taxonomy fields: 2,181 have one known value plus empty
unknown values; seven have two distinct known taxonomy IDs. It also found
39 main natural-key collision groups, all in the
source-scoped `name` namespace. The import policies above preserve their
published identities and source-specific information while keeping main's
strict natural-key index. The projector reports the affected entity counts.

Original entity occurrence trees, matched flags, resolver reasons, gene/state
projections and external resolution diagnostics are absent from the serving
files. They cannot be reconstructed by labeling all published entities as old
matches. Capability metadata must state this explicitly. The previously agreed
MetSigDB policy uses published chemical entities.

The serving contract also drops the original explicit ontology ID. Known
source/namespace mappings restore IDs such as `gene_ontology`,
`reactome_pathways`, `kegg_pathways`, `uniprot_keywords`, `enzyme_classification`,
`chebi`, and `swisslipids`. Unknown scopes use an explicitly named
`omnipath:published-scope:<resource>:<namespace>` scope and increment a reported
compatibility count; this is not claimed to be the original ontology ID.

Downstream adapters must use retained source-specific taxonomy where canonical
taxonomy is unknown, use actual member evidence for reaction direction and
compartments, and preserve distinct qualified assertions when building main's
interaction facts. `associated_with` gene rules are not inferred catalysis.
Those are consumer alignment requirements, not changes to the published files.

## Verification

The focused DuckDB fixture suite exercises duplicate aliases/annotations and
evidence, shared identities, scoped-name conflicts, known-taxonomy disagreement,
qualified statements sharing graph triples, source ownership, nonnumeric and
zero-padded row strings, numeric quantities with distinct comparators/prefixes,
null and empty lists/struct fields, known ontology metadata and source-scoped
ontology statements. PostgreSQL DDL/index and downstream consumer parity are
checked independently against the reused main backend.
