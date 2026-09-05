# Direct migration of legacy inputs to inputs_v2

Strategy agreed on 2026-09-05. Scope: complexes, enzyme–substrate relationships and interactions/networks. This document describes the implementation strategy; it does not claim that the listed resources have already been ported.

This strategy supersedes the execution architecture in the [legacy adapter proposal](inputs-v2-legacy-adapter-mapping.md). Keep the [resource inventory](inputs-v2-legacy-adapter-inventory.md) as the migration backlog, including its aliases, version distinctions and incomplete inputs. Field-level observations in the adapter proposal remain useful source-audit notes, but the implementation should follow this document and [Biolink modeling conventions](../pypath/inputs_v2/BIOLINK_MODELING.md).

## 1. Decision and target architecture

Port each resource into a native `inputs_v2` module that acquires source data, parses it into source-faithful rows, and explicitly maps those rows to the new schema:

```text
source file(s) / API responses
    → Download / shared download infrastructure
    → source parser → serializable parsed records
    → explicit Biolink mapper
    → Entity / Relation → silver extraction and existing build
```

The new module must not depend on calling the old public input function to obtain already mapped or aggregated objects. Legacy code is a reference for endpoints, file structure and documented interpretation. Reuse or extract sound parsing logic where appropriate, but remove its coupling to legacy ID mapping, `intera.Complex` / `DomainMotif` construction, process-global settings and core network aggregation.

This matters because a wrapper cannot recover source identifiers, isoforms, complex distinctions or individual evidence records that the old function has already removed. Direct ports preserve those distinctions before schema construction and use the native download/build lifecycle.

Each resource is modeled individually. Shared utilities should handle mechanics such as archive reading, publication normalization and aligned columns. Biological crosswalks, source repair rules and predicate decisions belong in the resource module or its parser, not in a universal resource-to-Biolink adapter or the shared schema engine.

## 2. Files and responsibilities

| Location | Responsibility |
| --- | --- |
| `pypath/inputs_v2/<resource>.py` | `ResourceConfig`, download definitions, explicit mapping definitions, datasets and module-level `resource`. |
| `pypath/inputs_v2/parsers/<resource>.py` | Optional source parser for custom formats, multiple files, nested data or source-specific preprocessing. |
| `pypath/inputs_v2/parsers/base.py` | Existing common readers such as CSV/TSV/JSON; use directly when sufficient. |
| `pypath/internals/tabular_builder.py` | Existing `FieldConfig`, `CV`, entity/relation/membership builders and generic execution machinery. Extend only for a demonstrated source-independent capability gap. |
| `omnipath_core` | Canonical schema, identifier registry and serialization/resolution contracts. Source-specific migration shortcuts do not belong here. |
| `test/` and fixture data | Parser, mapper and extractor tests using small representative source fixtures. |
| `docs/` | Per-resource coverage, semantic decisions, exclusions and snapshot validation results. |

A separate parser file is unnecessary for a plain table that an existing reader can expose faithfully. Conversely, do not hide archive traversal, XML parsing or multi-table joins inside a field-mapping callback. Use a parser subpackage only when resource complexity warrants it.

Keep the legacy module available for existing consumers during migration. Do not delete or change its public API as a side effect of adding the v2 resource. Removal or compatibility changes are separate work.

## 3. Audit a resource before writing its port

For each inventory entry, record:

1. **Resource identity:** actual upstream database, release/version, source/via provenance, aliases and whether a native v2 resource already exists. Add a missing dataset to the existing resource where appropriate instead of creating another wrapper for the same database.
2. **Acquisition:** files/API endpoints, authentication if required, archive members, source release and species options. Verify the current source or an available local snapshot before implementing acquisition; legacy URLs are leads, not proof of availability.
3. **Source contract:** row fields, missing-value conventions, identifier systems, entity types, evidence granularity and relationships among input tables.
4. **Legacy transformations:** ID translation, species filtering, isoform stripping, orthology inference, candidate expansion, complex merging, deduplication and discarded fields. Decide deliberately which operations are source parsing, which belong downstream and which should not be carried forward.
5. **Output contract:** named datasets, required participants, source-to-Biolink crosswalk, retained evidence and fields that cannot yet be represented in serving data.

Use `core/complex.py`, `resources/data/resources.json` and `resources/data_formats.py` to find intended legacy configurations. In particular, `NetworkInput` provides useful column, filter, reference and direction definitions. Treat those as migration specifications to inspect, not a runtime normalization engine to import into every new resource.

An unimplemented function or empty legacy output is not a source contract. For 3DComplex and HPMR complexes, audit source availability and biological evidence before planning an ordinary port.

## 4. Resource definitions and datasets

Create a real `ResourceConfig` with the existing metadata contract:

| Configuration | Requirement |
| --- | --- |
| `id`, `name` and resolved resource names | Reuse the canonical resource identity; add a resource registry entry only when needed. Follow the existing slug/short/full name conventions. |
| `url`, `description`, `pubmed` | Describe the upstream resource accurately; distinguish a resource citation from publications supporting individual records. |
| `license`, `update_category` | Use reviewed metadata and existing operational enums; do not guess permissions or update frequency. |
| `primary_category`, `annotation_ontologies` | Declare actual dataset scope and ontology dependencies when applicable. |

Operational resource/license enums remain acceptable in configuration. They must not become biological entity types or annotation vocabulary.

Expose a module-level `resource = Resource(config, ...)` with stable, meaningful dataset names such as `complexes`, `enzyme_substrate` and `interactions`. A dataset may emit `Entity` or `Relation` records as supported by the existing builders. Choose separate datasets when evidence meaning, source release, species selection or acquisition differs materially. Do not emit both a site-aware and a pair-only view as indistinguishable independent evidence.

For example, the following is the assembly pattern after the actual resource-specific configuration, download and mapper have been defined:

```python
from pypath.inputs_v2.base import Dataset, Resource
from pypath.inputs_v2.parsers.base import iter_tsv

# config, interactions_download and interactions_schema are defined above
# using this resource's reviewed metadata, source URL and biological mapping.
resource = Resource(
    config,
    interactions=Dataset(
        download=interactions_download,
        raw_parser=iter_tsv,
        mapper=interactions_schema,
    ),
)
```

This is a wiring example, not a complete resource implementation. Replace `iter_tsv` with the resource parser when necessary. Follow existing module-level discovery; do not build a second resource loader.

Versioned sources such as hu.MAP1/2/3 need explicit release identity and acquisition parameters. A current-release URL is not an immutable version. Preserve the acquired artifact/version in build provenance so outputs can be traced to the actual snapshot.

## 5. Downloads, parsers and parsed records

### Acquisition

Use `Download` and the shared download infrastructure for URLs, archive selection, encoding and refresh. Keep acquisition separate from biological mapping. Multiple-file/API resources may need a source-specific orchestration function using that same infrastructure; `download=None` is appropriate only when that orchestration genuinely owns acquisition, not as a shortcut to call the old input function.

Resource imports must not download data. Propagate `force_refresh` through every acquired artifact and any prepared cache. Do not copy a legacy pickle cache whose identity omits species, release, parser version or options. Verify that dataset argument changes select the correct artifacts and invalidate the applicable prepared data.

### Parsing

The parser follows the existing `(opener, **kwargs)` convention and yields dictionaries. It may decode formats, normalize source nulls, parse explicit accession prefixes, join source tables and expand genuinely separate source observations. Document each output field and its source.

Preserve original values alongside normalized values when conversion is significant: source action labels, numeric comparison operators, original identifiers, coordinates and evidence IDs. Keep records serializable; no legacy domain objects or arbitrary Python sets in persisted payloads. Sort unordered values deterministically while retaining meaningful source order and aligned member/coordinate arrays.

Reject malformed required fields with explicit diagnostics or counted exclusions. Do not silently turn a broken download, missing required archive member or changed column header into an empty successful dataset. Missing optional values remain missing; an unsupported biological value is distinct from a null.

### Identity and inference

Retain source identifier namespaces instead of forcing every identifier through the legacy UniProt mapper. Decode explicit source cross-references, but leave canonical resolution to the existing identity pipeline. A gene symbol, RefSeq nucleotide ID and UniProt protein accession are not interchangeable strings.

Preserve explicit species and isoforms. Do not reproduce orthology inference, candidate expansion or sequence-based reassignment as incidental parsing. If such a derived dataset is required, implement and label it separately with provenance. This includes the old enzyme–substrate processor's validation and PhosphoSite's non-strict cross-species behavior.

## 6. Explicit Biolink mapping

Use the installed pinned `biolink-model` and the current [modeling conventions](../pypath/inputs_v2/BIOLINK_MODELING.md), rather than assuming terms from a newer online model. Check class meaning, predicate domain/range, symmetry, enum values and deprecation status during implementation.

| Concern | Mapping rule |
| --- | --- |
| Biological classes | Pick the most specific native class supported by the source: `Protein`, `Gene`, `MicroRNA`, `MacromolecularComplex`, `ChemicalEntity`, etc. Do not infer a class solely from a display label. |
| Identifiers | Use `IdentifiersBuilder` with `Namespace`. Add namespaces only for actual identifier systems. Keep original accession identity and isoform distinctions. |
| Names and synonyms | Use `Namespace.NAME` / `Namespace.SYNONYM` identifiers, after real accessions, following the current application convention. |
| Relations | Use `RelationBuilder` with explicit subject, predicate and object. A broad association, physical interaction and causal effect are different assertions. |
| Qualifiers | Use native aspect, direction and mechanism slots/enums where source evidence supports them. Unknown qualifiers remain absent or require review; never invent a default aspect. |
| Membership | Use nested entities with `MembershipBuilder`, `Member` / `MembersFromList` and an appropriate native membership predicate. |
| Publications/taxonomy | Use `slots.publications` and `slots.in_taxon` with valid CURIEs. A publication count is not a PMID. Chemicals do not inherit the assay organism as their taxon. |
| Measurements | Use typed quantitative values and source-field/comparator context when appropriate, following the existing measurement contract. Preserve scale and units. |
| Unmapped structured fields | Retain in persisted parsed/evidence payloads and document the serving gap. Do not convert them to labeled description strings or mint convenience ontology terms. |
| Prose | Use `description` for actual source prose. |

Use `FieldConfig`, `CV`, `EntityBuilder`, `RelationBuilder` and membership builders for declarative extraction. Define reusable fields at module scope and keep transformations pure so the existing executor can safely cache them. Use source-local functions for more complex semantics; do not create biological fallback logic in the shared builder.

Each resource needs an auditable crosswalk of **source field/value → target class/namespace/slot/qualifier → rationale → unknown-value behavior**. Exact source mappings are acceptable; substring guessing is not. Review all observed categorical values in an available full snapshot, not only the few used in a fixture.

## 7. Category-specific ports

### Complexes

Parse source complexes directly instead of constructing and then unpacking `intera.Complex`.

| Source information | Target / decision |
| --- | --- |
| Source complex accession and name | `MacromolecularComplex` with real accession identifiers and a name identifier. Do not invent a source accession from a composition string. |
| Member identities | Nested typed members retaining source accessions, isoforms and species. Do not discard non-protein members simply because the old function did. |
| Membership | Native membership relation, usually `has_member`; keep the complex intact. No automatic all-pairs interaction expansion. |
| Explicit stoichiometric coefficient | Membership annotation with `slots.stoichiometry` when justified; verify its value survives extraction and serving. A member listed once does not establish experimentally known 1:1 stoichiometry. |
| References, methods, confidence and prediction status | Source evidence plus explicit supported attribute mappings. Preserve source-specific meanings and score scales. |

Retain distinct source complexes even if they have identical component sets; downstream consolidation must preserve source identities/evidence. For accessionless complexes, exercise the existing structural identity path before enabling the resource. Missing stoichiometry remains unknown.

Prioritize ComplexPortal, hu.MAP versions, PDB, Compleat, Havugimana, SPIKE and scConnect from the inventory. Inspect their underlying source tables, not just their legacy object outputs. In particular, direct parsing should avoid ComplexPortal's old isoform stripping and composition merges where the original data supports those distinctions. scConnect's annotation-derived groups still require a biological review before being asserted as physical complexes.

### Enzyme–substrate

Parse site-level source observations, keeping enzyme, substrate, identifier domains, both isoforms, modification, residue letter/position, sequence context, references and evidence source together. The internal field name should be `enzyme`, since phosphatases are included.

The initial supported relation is an explicitly qualified effect on the substrate. For example:

| Source assertion | Protein-level Biolink projection |
| --- | --- |
| Enzyme phosphorylates substrate | `affects`, phosphorylation aspect, increased direction. |
| Enzyme dephosphorylates substrate | `affects`, phosphorylation aspect, decreased direction. |
| Other modification | Exact reviewed modification/aspect/mechanism crosswalk; no global phosphorylation default. |

These directions describe phosphorylation, not protein activity. Do not infer activation or inhibition from the modification alone.

**Site representation is a completion requirement, not an incidental extra field.** The current silver relation lacks dedicated residue/site fields. A first implementation may provide qualified protein-level relations with preserved site evidence, but must retain distinct residue/position/isoform/motif occurrences through storage and consolidation. Do not call that full site-level serving support. Exposing sites as structured serving attributes or first-class sequence features requires a reviewed published schema mapping and end-to-end extractor/storage/API validation. If evidence persistence loses sites, the port is not ready to enable.

Source-specific checks include HPRD's RefSeq protein substrate field; MIMP's alternate nucleotide identifier, candidate expansion and publication count; iPTMnet's two isoforms and PTM type; DEPOD's dephosphorylation semantics; and ProtMapper's evidence/source attribution. Reassess the source files directly so the new implementation need not inherit legacy field omissions or inferred candidates. Source-wide constants such as phosphorylation must be justified for the selected dataset.

Do not equate a successful parse with sequence validation. Report whether coordinates are source-reported or independently validated, and keep isoform coordinates tied to their sequence.

### Interactions and networks

Translate each source's old `NetworkInput` specification into an explicit native parser/mapper. Preserve useful filters and source-defined dataset distinctions after reviewing their purpose. Do not instantiate the legacy core `Network` or blindly copy process-global filtering behavior.

| Legacy specification / source fact | Native port |
| --- | --- |
| Endpoint columns, ID domains and entity types | Explicit parsed fields and participant entity builders. Preserve nested complex endpoints. |
| Organism fields and filters | Source-aware species selection and participant taxonomy; reject unsupported species requests. |
| Sign and direction encodings | Exact source crosswalk to a reviewed predicate and qualifiers. Endpoint orientation, causal effect direction and predicate symmetry are separate concepts. |
| References and inclusion filters | Explicit dataset policy and counted exclusions; retain evidence occurrences. |
| Extra edge/node attributes | Reviewed slots, qualifiers or quantitative attributes; keep assertion-specific context off generic participant entities. |
| Source/via labels | Preserve aggregator provenance; do not double-count an aggregated observation as independent primary evidence. |
| Loop/deduplication/complex expansion options | Deliberate output policy; preserve qualified statements, sites and evidence. Do not inherit expansion or pair-only merging silently. |

For physical interaction evidence, select `interacts_with` or a more specific direct physical predicate only when warranted. Functional association datasets such as STRING links are not uniformly physical. Regulatory resources such as TRRUST and CollecTRI require source-supported effects/aspects; TF binding alone does not establish changed expression. Ligand–receptor labels alone do not establish agonism. RNA and drug–drug datasets need their actual participant domains and semantics.

Negatome must retain explicit negative evidence through the full storage/consumer path; it must never produce positive interactions or encode non-interaction as inhibition. Structural interface/domain APIs, drug-pair adverse-event data and process-derived reaction edges require dedicated mappings when a binary relation would discard essential context.

## 8. Evidence, reproducibility and intentional differences

Every port should trace a mapped observation back to source resource, dataset/release, artifact and source row/evidence identifier where available. Record parser/mapping versions and acquisition parameters in the build's provenance/cache identity. Retain source values for excluded or unmapped information in the appropriate audit/evidence artifacts, with counts and reasons.

Verify actual evidence persistence through the build. Extra keys returned by a parser do not establish that those fields reach persisted evidence. Likewise, two rows with the same endpoints but different sites, effects, methods or publications must not silently collapse into one indistinguishable occurrence.

Legacy output is a useful comparison baseline, not the correctness oracle. A direct port may retain more isoforms, preserve complex distinctions, remove unsupported inference or apply corrected parsing. Explain those differences using source evidence and representative examples. Do not force row-count parity by reproducing known information loss.

## 9. Implementation sequence and acceptance criteria

Implement small, reviewable resource ports using this sequence:

1. Select a resource from the inventory and confirm it is absent or identify the missing dataset in an existing v2 module.
2. Audit an available source snapshot and write its acquisition, parsed-row and semantic crosswalk contracts.
3. Add native downloads and the smallest necessary parser; test representative source-format fixtures.
4. Add the explicit mapper, metadata and `Resource`/`Dataset` definitions.
5. Run records through `SilverExtractor` and the relevant persistence/consolidation path, addressing schema gaps before enabling affected outputs.
6. Audit an available full snapshot and document counts, differences, exclusions and unresolved fields.
7. Integrate discovery/build configuration; create immutable resource versions and a local release when performing a completed migration under the existing modeling workflow. Publication is separate from writing or testing a port.

Suggested initial ports are one complex source (ComplexPortal or a hu.MAP release), one site-aware source (KEA, followed by DEPOD or iPTMnet), and one simple network source followed by a regulatory source (BioGRID subset, then TRRUST). Selection depends on available snapshots; do not block the whole migration on an inaccessible endpoint. Extend shared mechanics only after concrete ports demonstrate the repeated need.

A resource is ready only when the following are recorded and checked:

- **Acquisition:** imports are inert; source artifacts and refresh behavior are correct; species/release/options select the intended input.
- **Parser:** representative format fixtures, required-column failures, optional fields, nulls, archive members and multi-table joins are covered. Source identifiers/evidence survive parsing.
- **Mapper:** native types, namespaces, predicates and qualifiers match source meaning; unknown values are counted/reviewed; no invented direction, aspect, activity, taxonomy or identifier equivalence.
- **Category details:** complex membership/stoichiometry/accessionless identity; PTM sites/isoforms and opposite modification effects; network signs, filters, mixed participant types and negative evidence as applicable.
- **End to end:** extractor and persisted evidence preserve membership, qualified statements, distinct source observations and quantitative/site context. Tests must cover what must not be asserted as well as expected records.
- **Snapshot audit:** source rows, parsed rows, mapped observations, exclusions by reason, unresolved identifiers and serving-field gaps are reported. Identify partial snapshots and unavailable data explicitly.
- **Coverage:** the resource's promised datasets are listed; existing legacy APIs remain usable; differences from legacy behavior and any remaining schema work are documented.

The migration deliverable is a native, source-faithful resource module with evidence-backed mappings and verified build output. The inventory tracks coverage; a successful import or a wrapper around legacy output does not establish completion.
