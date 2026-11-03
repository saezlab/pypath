# Entity resolution and molecular forms

Decision record · 4 October 2026 · grouping clarified 5 October

This document captures the decisions from our schema and explorer discussion.
The implementation follows these decisions; production validation is tracked
separately. The latest clarification takes precedence over earlier alternatives: **the reference entity
key identifies the gene, while the main-level `entity_type` retains the type
reported by the source, such as protein, gene or RNA.**

The implemented record layout and compatibility rules are specified in section 9.

## 1. Agreed direction

- Own entity resolution within the omnipath-build system, including reference
  preparation. Resource builds should not depend on a second resolution step in
  omnipath-utils or in the serving database.
- Use NCBI Gene as the shared reference for genes and their products wherever a
  supported gene mapping is available.
- Preserve the source's entity type at the main level. Resolving a protein to a
  gene reference does not change its `entity_type` to `gene`.
- Resolve protein identity separately to primary UniProt, retaining it in the
  molecular form. Preserve specific transcript identity when available.
- Standardize a structured `molecular_form` in `inputs_v2` and carry it through
  the build into the relation evidence (`relation_evidence`).
- Keep reusable protein references, but do not create an entity for every
  observed combination of isoforms, variants and modifications.
- Present a gene-centered explorer from which users can inspect products,
  observed forms and the knowledge attached to them. Advanced molecular-form
  search can follow later.

## 2. Reference identity, entity type and molecular form

These fields answer different questions:

| Concept                | Meaning                                                                                                                        |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `reference_entity_key` | Which gene does this entity refer to? Use the NCBI Gene reference where resolved.                                              |
| `entity_type`          | What does the source describe: a gene, protein, RNA, transcript, etc.? This remains a main-level field.                        |
| `molecular_form`       | Which resolved product and, where specified, which isoform, sequence, variants or modifications does this occurrence describe? |

A record can therefore have a gene reference identifier and
`entity_type = protein`. These are intentional, compatible statements. The
identifier's namespace does not determine the source's molecular type.

The source type should use our standardized vocabulary while preserving the
source's meaning. A generic or unknown source type must not acquire unsupported
specificity merely because its identifiers resolve successfully.

We are **not** adopting the earlier suggestion to replace the main type with
`gene` and move the source type into a new `participant_type` field or into
`molecular_form`.

### Illustrative records

The keys below are readable identifiers, not a specification of the final
serialized keys. The internal namespace `entrez` denotes NCBI Gene. A product
key may ultimately be the existing hashed key of a reusable protein record.

```yaml
# Source describes a specific protein isoform.
reference_entity_key: entrez:7157
entity_type: protein
molecular_form:
  protein_entity_key: uniprot:P04637
  isoform_identifier:
    ns: uniprot
    id: P04637-2
```

```yaml
# Source describes the gene itself.
reference_entity_key: entrez:7157
entity_type: gene
molecular_form: null
```

```yaml
# Source describes a protein, but supplies only a gene identifier.
reference_entity_key: entrez:7157
entity_type: protein
molecular_form: null
```

The third record remains protein-level information even though no specific
protein reference has been established. An absent molecular form therefore does
not, by itself, mean gene-level evidence.

## 3. Resolution at gene and product levels

Resolution should produce a coordinated result containing the gene reference
and any supported product identity. These are two levels of one build process.

1. Preserve the source entity type and specific molecular information before
   identifier normalization removes suffixes or versions.
2. Resolve gene identity to NCBI Gene and protein identity to primary UniProt
   using evidence appropriate to each level. Preserve transcript identities
   through the relevant transcript mappings.
3. Connect products to genes through explicit gene–product mappings, checking
   consistency with any gene identifiers supplied by the source.
4. Preserve unresolved identities and ambiguity rather than inventing a more
   specific participant.

The reference data already provides useful links: UniProt-to-GeneID,
Ensembl gene/protein/transcript mappings, HGNC cross-references and RefSeq
accessions. Gene symbols require organism information and ambiguity handling.

The relationships must remain distinguishable. A secondary UniProt accession
can identify the same protein entry; a GeneID identifies its gene. Sharing a
gene does not make two protein entries identical.

### Projection rules

- A supported protein-to-gene mapping supplies the common gene reference while
  preserving the protein identity in the molecular form.
- Gene-level evidence must not be expanded onto every protein associated with
  that gene.
- A representative or reviewed protein selected for display must not become
  the asserted participant of a source record.
- A source that describes a protein but provides only a gene identifier can
  remain protein-typed without a resolved protein identity.
- A transcript can map directly to a gene; it need not resolve through a
  protein. This accommodates noncoding products as well.
- Gene–product mappings must support missing and multiple mappings. A single
  gene reference should only be assigned when justified by the available
  evidence. Conflicting identifiers must remain diagnosable.
- When a gene cannot be established, retain a supported primary UniProt or
  source identity as an explicit fallback. Do not manufacture a gene identity
  or lose the reported type.

Fallback and ambiguous resolution use the representation in section 9.

## 4. What belongs in `molecular_form`

The structure should support:

- A resolved protein or transcript reference, where established.
- Isoform or sequence identifiers that specifically identify the reported
  molecular participant.
- Structured post-translational modifications.
- Structured variant or mutation information.

The source entity type stays outside this structure. General labels and alias
collections continue to use the existing entity identifier and annotation
mechanisms.

Modifications and variants must be individually addressable so that future
queries can find a particular modification or variant across different observed
combinations. Those combinations remain attached to their original occurrence:
evidence for variants X and Y must not be combined with evidence for X and Z to
invent an X/Y/Z form.

Missing information means unspecified. An absent isoform identifier does not
assert the canonical isoform; an absent modification or variant does not assert
an unmodified or wild-type molecule. A partial description is not necessarily
the identity of a fully specified molecular species.

We will not enumerate possible molecular forms or mint an entity for each
PTM/variant combination. A dedicated reusable state model can be reconsidered
if future pathway or state-transition use cases require it.

### Identifier retention

The distinction is biological specificity, not whether an identifier happens to
be present in the reference catalogue.

- Preserve identifiers in the molecular form when they specifically establish
  that occurrence's protein, transcript, isoform or sequence identity.
- Preserve a specific isoform identifier even when the reference protein's
  identifier collection already contains it.
- Keep form-specific identity attached to the occurrence. A catalogue link to
  an isoform or transcript must not imply that it is interchangeable with the
  general gene, or that every occurrence of that gene describes that form.
- Do not use “source identifiers minus canonical identifiers” as the retention
  rule. That would discard specificity whenever the catalogue already knows
  the identifier.
- Do not interpret every source cross-reference as an assertion that the
  corresponding transcript or isoform participated.
- General entity identifiers can retain their existing handling. Keeping a
  complete duplicate of every original identifier list on every evidence item
  is not a requirement of this design.
- Raw payloads remain useful for source provenance, but information needed to
  identify the molecular participant must survive in the structured output.

Sequence-specific positions retain their supplied coordinate reference, including
its version, coordinate system and indexing convention. An unknown reference
stays unknown; it is not replaced with the selected primary protein.

## 5. Parquet outputs

The published tables (see `core_documentation/entities.md`):

- `entity` with `entity_identifier` and `entity_annotation`: reusable
  entity/reference records, identifiers and annotations.
- `relation` with `relation_annotation` and `relation_evidence`: relations and
  their source evidence, including paired subject/object molecular forms.
- `evidence_payloads`: original source payloads for deeper inspection.

We previously agreed to keep reusable protein reference records, independently
resolved to primary UniProt, with explicit gene–product links. Specific
transcript references can follow the same approach. These references provide a
place for product identity and annotations without creating state-combination
entities.

The latest decision also requires main-level source types on gene-referenced
records. **A shared gene reference key must not be assumed to be the complete
unique key of every source-typed or product-specific record.** How those records
and the reusable protein references coexist is specified in section 9.

For relations:

- Use the shared gene reference for biological relation grouping where the
  mapping is supported.
- Preserve the subject and object molecular forms together on each evidence
  occurrence. A union of all forms across a relation would lose which
  participants were observed together.
- Keep `entity_type` at the main level. Aggregation must not overwrite
  protein/RNA/gene distinctions with one arbitrarily selected type.
- Preserve existing predicates, direction, signs and meaningful qualifiers.
  Molecular form combinations should not automatically become new entity keys.
- Any symmetric endpoint reordering must move the corresponding molecular
  forms and endpoint context together.
- Explicit structural links such as gene–product mappings must retain their
  actual meaning; projecting both ends to one gene must not turn them into
  meaningless self-links.

Different source types produce separate entity and relation rows. Their shared
gene reference provides the common browsing group without erasing those types.

## 6. Explorer and query behavior

The intended interface has search, results, filters and entity/relation tabs.
The default experience is gene-centered: a user can open a gene and explore its
gene-level knowledge, known products and observed molecular forms.

- Preserve visible entity-type distinctions from the sources, including
  protein and RNA. Grouping by gene must not label all results as gene evidence.
- Use one Group toggle. When enabled, apply gene-reference grouping to
  biological records and connectivity grouping to chemicals on the same
  results page. There is no user-facing choice between these grouping rules.
- Distinguish catalogue-known products from forms actually supported by
  resource evidence.
- Opening a specific protein or isoform must select evidence connected to that
  protein or isoform. It must not silently include every relation attached to
  the gene.
- A future combination of form filters must match the intended endpoint in the
  same evidence occurrence, rather than unrelated annotations elsewhere in
  the aggregated relation.
- Exports must retain the molecular context and any referenced product records
  needed to interpret the selected evidence.

Advanced molecular-form search, filtering and dedicated indexes are deferred.
The HTML prototype is a UX experiment and includes earlier design ideas; it is
not the schema specification for these decisions.

## 7. Serving, PostgreSQL and performance

The current application serves Parquet through DuckDB. We discussed PostgreSQL
as a downstream representation, but did not decide to replace the current
serving architecture or agree on a concrete PostgreSQL migration.

The subsequent [migration-branch integration](reports/molecular-parquet-migration-integration-20261005.md)
adds downstream PostgreSQL context tables and a COSMOS label correction. The
Parquet web architecture and independent resource builds remain the same.

Resolve identities and compute gene–product mappings during the build. Store
the selected gene references in the outputs so ordinary browsing does not have
to repeat entity resolution per request. A downstream PostgreSQL representation
should consume these resolved identities and preserve the same semantics.

Structured molecular forms add evidence data without requiring an entity for
every state combination. Exact protein/form filtering may eventually benefit
from a flattened serving projection or database indexes. Such projections can
be derived from the authoritative evidence.

The earlier local performance experiment showed that a flattened annotation
projection could accelerate a specific filter. It did not benchmark this final
schema or establish a production latency or storage guarantee.

## 8. Implementation plan

Validate SIGNOR and IntAct end to end first, then extend coverage. Keep advanced
molecular-form search and state-combination entities outside this first release.

### 1. Shared schema

Define the gene reference, main-level source `entity_type`, product references
and `molecular_form` in `omnipath_core`. Settle record keys, relation grouping,
standalone entity evidence and missing/ambiguous values using the choices below.
Version the output schema.

### 2. `inputs_v2`

Add the optional molecular-form structure to input entities and builders. Keep
the existing source type. Populate specific identifiers, isoforms, PTMs and
variants before normalization loses their details. Include form information in
builder caching so different occurrences are not accidentally merged. Start
with SIGNOR and IntAct, then update other relevant mappers.

### 3. Resolver and reference construction

Build explicit gene–product mappings and lookup paths for each identity level.
Update the Rust decision policy and Python result handling to return a gene
reference and any supported protein/transcript identity together. Remove the
current gene-to-protein expansion for gene-level evidence. Preserve the source
type, unresolved identities and ambiguity; never select a participant merely
because it is the gene's preferred protein.

### 4. omnipath-build

Carry the new fields through extraction, resolution and Parquet writing. Write
reusable product references and their gene links to the agreed entity layout.
Keep subject/object forms paired within each relation evidence item, including
during endpoint reordering. Preserve standalone molecular evidence and specific
form identifiers. Keep the three existing output files.

### 5. Serving interfaces

Update API contracts, DuckDB queries, the web explorer, Python client and
exports to read the new schema. Group by gene reference while filtering by the main-level source
type. Show forms in evidence details and restrict product navigation to the
matching evidence. Preserve referenced product records in exports. Rebuild
affected serving projections. Any downstream PostgreSQL loader must retain the
same identities and context without resolving again.

### 6. Validate and release

Check gene-only records, protein/RNA records, isoforms, variants, missing and
ambiguous mappings, mixed source types, standalone evidence and endpoint swaps.
Verify that gene evidence is not copied onto proteins and distinct forms are
not merged. Compare output size and search, filter and evidence-query timings
against the current build. Rebuild pilot resources, then publish a compatible
release; older missing context must remain unknown.

## 9. Implemented schema choices

- **Unique rows:** retain the existing hash of `(entity_type, namespace,
  identifier)`, including existing source/taxon scoping for unresolved names.
  A protein source resolved to gene 7157 has type `protein`,
  namespace `entrez`, identifier `7157`; a gene source has a separate
  gene-typed row. Both have `reference_entity_key = entrez:7157`.
- **Reusable products:** primary UniProt proteins and specifically reported
  transcripts have their own typed entity rows. Their `gene_reference_keys`
  contain supported catalogue links. An observation's `molecular_form` refers
  to those rows through `protein_entity_key` or `transcript_entity_key`.
  If the catalogue cannot resolve an explicitly reported protein identifier,
  retain that exact source product, including its supplied version, isoform or
  chain suffix. Mark the occurrence `omnipath:protein_mapping_status = reported`;
  this does not assert a primary UniProt mapping or supply catalogue gene links.
  Ordinary cross-references alone do not establish a product participant.
  Without catalogue resolution, competing product assertions remain unselected,
  including a separately supplied entry and isoform; preserve their specific IDs.
- **Relations:** retain typed endpoint keys, predicate and qualifier identity.
  Store subject/object gene references as scalar fields and paired molecular
  forms inside each evidence item. Gene observations never fan out to products.
- **Forms:** optional isoform identifier, sequence identifiers, individual
  modifications and variants. Features retain exact positions only when supplied,
  plus coordinate identifier/version, system and position base. Partial or empty
  fields mean unspecified. Explicit gene-level genomic variants are preserved;
  they do not cause protein or transcript inference.
- **Fallbacks:** use a native product/source CURIE when no supported unique gene
  assignment exists. Catalogue links do not override conflicting occurrence
  evidence. If a fallback and reusable product share a row key, retain the
  fallback reference and keep known catalogue gene links separately. Exceptional
  outcomes remain on the affected evidence as `omnipath:gene_mapping_status`
  and `omnipath:gene_mapping_candidate` annotations.
- **Standalone observations:** `entity_evidence` retains each source
  occurrence, its annotations and resolved form. General entity identifiers
  exclude specific isoform/transcript/sequence aliases except when that identifier
  is itself the row's canonical identity.
- **Complexes:** preserve membership occurrence forms and distinguish explicitly
  different product/member compositions. No additional protein state entities
  are created from PTM or variant combinations.
- **Navigation:** reference view groups supported scalar gene references, or
  selects the exact typed row for a native fallback. Product view selects only
  occurrences whose molecular form explicitly names that product; an optional
  isoform filter applies to the same endpoint and occurrence. Choose the view
  explicitly, since a native reference and a product can share a row key.
  Older evidence without product pointers remains unspecified. Pagination
  covers both relations and standalone observations.
- **Compatibility:** build manifests retain `schema_version = 1` and use
  `serving_schema_version = 4`; disposable serving projections use
  `.serving/v2`. Older absent molecular context remains unknown. Build a fresh
  coordinated gene/product reference; old protein-projected identity components
  cannot be reused as the new gene identities.

The [resource audit](reports/molecular-form-input-coverage-20261004.md) records
which mappers retain molecular detail and which require further source work.

## 10. Code and discussion references

- [Current Parquet schemas](../packages/core/omnipath_core/schema.py)
- [Current entity and relation keys](../packages/core/omnipath_core/keys.py)
- [Shared input schemas](../packages/core/omnipath_core/silver_schema.py)
- [UniProt reference identifiers](../packages/build/omnipath_build/hubs/sources/uniprot.py)
- [Gene–product mappings and identity decisions](../packages/build/omnipath_build/identity/decisions.py)
- [Migration-branch integration](reports/molecular-parquet-migration-integration-20261005.md)

These code links describe the implementation. The coverage audit distinguishes
fixture validation from completed resource rebuilds; production performance
claims require measurements from those rebuilds.
