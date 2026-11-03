# PostgreSQL parity with main

The existing resolved Parquet files are the input contract. Their schemas, resource versions and independent update schedule stay as they are. This work replaces the PostgreSQL projection and downstream implementation with the current main contract, allowing explicit Biolink vocabulary mappings and truthful provenance additions.

Reference: omnipath-build main `9f9bb709c764e2721a0bcd881cf7ef0eb22bed9a`, with pypath dependency `33f37fbaab59993d24f5c32bb4b3e7587085bcf4`. `main-reference/9f9bb709c764` is an untracked frozen source snapshot. The identical `legacy/postgres` copy can serve as the oracle when deploying the tests; the test verifies hashes of the main schema, indexes, interaction and COSMOS algorithms/SQL, identifier registry and classification map. The oracle uses its exact, privately imported controlled-vocabulary dependency under `tests/oracle_cv_terms` and never replaces migration pypath.

## What must stay equivalent

| Area | Required behavior | Evidence to compare |
| --- | --- | --- |
| Public database structure | Main columns, types, keys, foreign keys, partitions, indexes, operator classes, partial predicates, extended statistics and query functions | Definitions from the executed PostgreSQL catalogue; counts or object names alone are insufficient |
| Canonical entities | Keep each published resolved entity identity; make the main canonical key explicit and unique | Published-key crosswalk, canonical-key conflict report, endpoint coverage |
| Identifier dictionary | Main stable identifier type IDs; one global identifier per type/value and deduplicated source links | Exact type registry, identifier identity and source-owned link samples |
| Annotation dictionary | Shared term/value/unit identity, deduplicated links, preserved scope and source ownership | Duplicate-sensitive dictionary/link samples; differing quantities must remain distinct |
| Canonical relations | One graph triple per canonical subject/predicate/object | Collapse differing published qualifier keys onto the graph triple while preserving each evidence assertion |
| Resource facts | Main record identity: ordered endpoints, class, source **name**, and three-valued assertion signature | Full rows, deterministic UUIDs, opposite directions, two signs from one source, an unsigned neighboring source |
| Interaction headers | Main UUID over class and sorted participant **multiset**; source union | Header UUID, class, arity, complete participants and contributing sources |
| Reactions | Main chemistry grouping across resources, enzyme parties, per-role stoichiometry and compartments | Two resource event IDs for the same chemistry; same cargo consumed and produced with different numbers/locations |
| Transport | Join input and output spokes for the same event/member/source; avoid combining claims from different resources | Finished transporter fact orientation and `transport` attributes; same cargo in two compartments |
| Ontologies | Main ontology relation/term tables, non-hierarchical axioms, closure, labels and synonym lookup | Child/parent orientation, subclass and part-of examples, separate ontology identity, a non-hierarchical axiom |
| Search, counts and classification | Main source decomposition, labels, canonical type gates, classification rules and bitmap outputs | Representative result rows and bounded planner checks; flag missing original resolver inputs |
| MetSigDB | Main product structure and extraction semantics, using the accepted eligibility policy: published chemical entities | All five source fixtures, complete membership rows and source/version metadata |
| Presets and network views | Main registered defaults and stored resource/organism/license scopes | Complete registry rows/defaults; filtering and folding belong to the separate API consumer |
| COSMOS | Main chemistry grouping, reversible behavior, orphan events, translated connectors, labels, sources and namespaces | Complete edge/label multiset and statistics, rather than edge counts only |
| Publication and retries | Main provenance plus the migration's separate resource/release schedule; complete products checkpoint together | Base identity, release manifest and committed product metadata; retry only remaining products |

The deliberate resolution extension is status `5 = published`. It must not claim that the old resolver matched, rescued, or rejected an original occurrence. Original status IDs, reason fields and constraints remain; the entity-bearing status check explicitly admits `5`. The independent structural comparison normalizes only this addition and the declared gene/protein type vocabulary literals.

## Identity and provenance issues already observed

The current release audit found **2,181 shared published entity keys with one known taxonomy plus empty values, and 7 keys with genuinely conflicting known taxonomies**. Preserve those keys. Keep the sole known taxonomy when the other records are empty; canonical taxonomy is NULL for the 7 genuine conflicts. Retain source-owned organisms for facts and queries. Choosing an arbitrary taxonomy or treating an empty value as a conflict would misstate the published data.

The audit also found **39 canonical-key conflicts among resource-scoped fallback names**. Keep distinct published entities distinct. Encode the published scope in the fallback canonical identifier/key and retain the readable original name and identifiers. Dropping resource scope, silently merging entities or weakening the main unique index is unacceptable.

Published row IDs are text and may be prefixed; main's row ID is bigint. Use a deterministic surrogate with a crosswalk to the original published row/upstream identifiers. Casting to NULL or reusing zero discards provenance.

Published relation keys can differ only in qualifiers. Main's graph-triple key cannot represent that distinction, but its evidence rows can. Deduplicate the graph triple and keep each source's scoped qualifier, measurement, citation and assertion signature.

## Explicit Biolink mappings

| Published input | Main downstream reading |
| --- | --- |
| `chemical_entity`, `small_molecule` | Accepted published chemical eligibility; do not include every descendant of `molecular_entity` |
| `protein`, `gene`, `rna_product`, `transcript` | Use the concrete type and available identifiers; do not rerun gene anchoring |
| `molecular_activity` | Reaction-event type gate; a pathway with similarly named membership predicates stays a pathway |
| Event `has_input` / `has_output` member | Reactant / product, even when role annotations are absent |
| Event `enabled_by` enzyme | Catalyst; reverse only the **finished main fact** to enzyme → event, retaining the Biolink graph direction |
| Enzyme `catalyzes` event | Catalyst with existing enzyme → event orientation |
| `stoichiometry`, `biopax:cellularLocation` | Numbers and locations from the membership's own relation evidence and role |
| `biopax:conversionDirection` | Read the source-owned event or member claim; do not reverse already oriented KEGG participants a second time |
| Relation-scoped `object_direction_qualifier` | Increased → stimulation TRUE; decreased → inhibition TRUE; unasserted opposite flag stays NULL; include supported enum descendants |
| Subject/object-scoped `connectomedb:participant_role` | Ligand/receptor roles for that evidence assertion, respecting any symmetric endpoint swap |
| miRBase `derives_from` with mature/precursor namespaces | Source-specific precursor → mature maturation fact; other sources/namespaces retain their generic claim, including when the graph triple is shared |
| `publications = PMID:…` / `doi:…` | Main PubMed/DOI hot columns with their prefixes removed; retain the original relational annotation |
| BAO Ki/Kd/IC50/EC50 | Affinity endpoints only; preserve units and comparators; KON/KOFF are different quantities |
| `has_quantitative_value` with `source_field=pchembl_value` | Main pChEMBL hot column; pH, temperature and flux bounds must not enter it |
| `has_confidence_score` with `source_field=combined_score` | Main score hot column from STITCH; `action_score` and ChEMBL assay confidence remain relational annotations |
| ChEMBL evidence dataset `mechanisms` | Explicit mechanism-of-action origin; an activity merely mentioning a mechanism does not prove this flag |

The `annotation_quantity` companion retains numeric value, unit, prefix, relation, comparator, source field and original published value. Different measurement semantics must not share one dictionary identity merely because their visible numeric text is equal. Unit conversion for a hot column must be explicitly defined and tested; retain the original quantity.

## Genuine limits to report

These are input limitations to measure and disclose, rather than reasons to change Parquet outputs or fabricate evidence:

- Original occurrence-level resolver decisions, parent/member occurrence identities and match reasons are absent from aggregated published entities. Canonical relation endpoints and scoped published annotations remain usable. Any aggregate entity-evidence row must identify its published origin honestly.
- ChEMBL's folded increase/decrease qualifier cannot recover every original agonist/antagonist/orthosteric/allosteric class. Retain precise classes when a published controlled annotation supports them; do not infer a precise class from sign or prose alone.
- `OntologyRelation.ontology_id` is absent from the published ontology row. The twelve audited current providers each declare one explicit ontology family for their ontology records; recover that input contract independently of a resolved endpoint namespace. Unknown providers keep a visible fallback/capability gap. The legacy base `ontology_terms` table stays empty, matching fresh main; main derives serving `entity_ontology_term` from normalized edges, identifiers and annotations.
- Ontology descriptions include both source definitions and comments. A readable description can populate the supported main field, but the exact original definition/comment distinction is unavailable.
- Some original classification inputs are not published: HMDB's exact superclass/class rank priority and the old Recon3D/metatlas subsystem fields need a documented fallback or a capability gap. Generic `has_topic` EC references do not replace these concepts.
- Published gene references and molecular occurrences are retained in the required context companions. Their presence does not recover original main resolver/state diagnostics or external lipid utility results. Preserve supported published facts and mark unavailable capabilities with reasons.
- Main's existing pair/reaction identity collision behavior remains part of this parity target. Fixing that design is separate work.

## Bounded verification

The focused checks are optional development verification, not a restored mandatory full-release rebuild/rollback phase. They never download a resource, invoke a resolver, mutate existing release schemas or print source data or credentials.

Local static and oracle checks:

```sh
uv run --frozen --no-sync pytest packages/postgres/tests/test_main_parity_contract.py packages/postgres/tests/test_main_metsigdb_parity.py -q
```

For the bounded PostgreSQL checks, the parent task sets `OMNIPATH_MAIN_PARITY_DSN` privately to the **new isolated migration PostgreSQL**, then runs the same command. Do not use production or the prototype database. The fixture creates and removes only random `parity_main_*` / `parity_port_*` schemas. It compares the exact frozen main algorithm against the adapted backend over declared resolved SQL facts, with no parser or original resolver involved.

Current harness covers:

1. Exact stable identifier IDs, source partition/primary-key contracts and deterministic interaction identity SQL.
2. Main static DDL and the executed catalogue: columns/defaults/nullability, constraints and their validation state, indexes and their readiness/validity, partition definitions, extended statistics, views and query functions.
3. A six-entity fixture with shared citations, opposing signs from one source, another inhibitory source, a silent source, reverse orientation, an unsigned pair and scoped ligand/receptor roles.
4. Merged cross-resource reaction chemistry; both catalyst orientations; the same cargo on both sides with different stoichiometry/compartments; a pathway gate; a 72-party reaction exercising main's wide-key hash indexes.
5. Equivalent affinity concentration units, distinct generic quantities, source-field-specific pChEMBL, source-owned STITCH combined score/citations and an explicitly published ChEMBL mechanisms dataset with no generic annotations.
6. Complete COSMOS edge rows and statistics over shared chemistry, explicit member/event direction mappings, orphan events, reversibility, connectors, transport and existing chemical aliases, with external mappings disabled. Catalogue definitions are compared again after scientific derive to cover late interaction indexes and statistics.
7. Exact main ontology term/definition/synonym/alias output, two sources for one ontology, source-scoped shortest ClassyFire closure with preserved non-hierarchy axioms, and primary gene-symbol ranking by attestations, length and alphabetic order.
8. MiRBase maturation orientation/class with wrong-source and wrong-namespace guards, plus a shared Biolink graph triple retaining both source-specific MIR maturation and a generic claim. Source IDs above the smallint range test the main bigint identity contract.
9. All five MetSigDB extraction/publication algorithms over exactly 20 synthetic source rows: complete membership columns and provenance, duplicate record selection, both pathway orientations, KEGG input/output deduplication and overview maps, MACdb structured subtype, source-scoped ClassyFire hierarchy and direct-assignment precedence, set taxonomy/name selection, chemical identifier display and structure groups. The same fixture checks exact executed product DDL/indexes, idempotent publication and resource-local stale removal under an unchanged build stamp. Main matched occurrences are compared with published status 5; only the approved chemical eligibility/type vocabulary is normalized.
10. Exact main network preset descriptors and complete registry rows, with only the declared concrete chemical type gate adapted to Biolink. The separate `test_main_network_presets_parity.py` fixture covers old-registry migration/defaults, apply/refresh and complete upserts, stored evidence/license scopes, ordered composition, mandatory/default attributes, timestamp stamping and literal loaded-resource warnings. Five bounded facts retain human/mouse organisms, source-local signs and unknown-license metadata unchanged. Main expects `connectomedb2025` and `metatlas`; `connectomedb` is a module spelling and `humangem` a display label. The real pinned release already uses the two expected resource keys.

11. Actual synthetic Parquet quantities differing only in numeric value, unit, prefix, comparator, binary relation or source field, including absent versus all-null quantity and duplicate occurrences; copied ontology axioms, known scopes, non-hierarchy predicates, source-owned labels/synonyms and complete frozen-main serving terms/shortest closure over fifteen entity/statement records.

12. A bounded nine-record Parquet fixture compares the complete four-column relation annotation link multiset with the frozen pre-batching query, plus independent decoded expectations. It covers nullable versus empty source/dataset attribution, evidence ordinal ownership, normalized and custom scopes, publisher evidence-union suppression, synthetic claims, resource/version isolation and NULL-term observations. Batches of 1/2/7/128 complete-owner groups introduce no duplicate links; the real loader observer callback is exercised.

13. The real SilverExtractor/ParquetWriter ontology literal, source family independent of resolved namespace, and the exact fresh-main empty ontology_terms base. Complete serving rows match frozen main over eligible bare/IRI/CURIE identifiers, original definitions and multiple ontology roots. Negative source/type/namespace and existing chemical/protein CV cases retain main priorities; scalar Name fallback requires a meaningful globally absent name/synonym and excludes identifier-valued labels. Computed lookups add no occurrence or authority evidence.

14. Default eighteen-table versus retained twenty-three-table loads over eighteen synthetic entity/statement records: complete normalized DuckDB and PostgreSQL rows, nonempty serving derivation, typed quantities, source-scoped evidence, required molecular companions and exact main constraints/indexes. Optional copies retain duplicate/order/raw inspection facts; default mode does not prepare their final tables, DDL or COPY. Published input counts are separate from copied counts. An older retained checkpoint finishes without re-projection or COPY.

15. Exact SHA-verified main Standard InChI inclusion rule over five synthetic source records, default/retained dictionary and links, unchanged canonical identities/quantities and complete excluded audit occurrences with NULL dictionary references. Accepted long names/synonyms/SMILES/unknown strings remain intact; actual local COPY and complete MetSigDB publication exercise the exact `Smiles:MI:0239` namespace and declared eligibility adaptation.

16. Serving-schema-4 fixtures retain virtual gene references, ordered gene links,
paired subject/object forms, specific product/isoform identifiers, feature
coordinates, repeated evidence and standalone observations in DuckDB and actual
PostgreSQL JSONB COPY, including default loads without optional audit tables.
Content reset clears the companions. COSMOS keeps GeneID catalyst and connector
labels in both published-alias and optional utils paths; genuine UniProt products,
chemical compartments, reversible edges and chemical-only MetSigDB stay covered.

Compare duplicate-sensitive full rows and content identities. Preserve distinctions between NULL and FALSE, source-local and merged claims, and participant sets versus multisets. Equal counts or an empty oracle result are not acceptance.

For the real release, sample supported examples from multiple resources after projection and each downstream product. Use bounded aggregate queries and query plans, inspect FK/index validity and release metadata, and record observed storage and timings. If a fixture or real sample exposes a mismatch, fix the mapping before the full build. Do not rerun already committed products to diagnose a later phase.

The independent Parquet fixture now compares complete structured measurement UUID/quantity/occurrence rows and actual copied ontology axioms, followed by frozen-main serving ontology facts and shortest closure. Published definitions/comments remain subject to the documented information limit. The main preset/network-view build boundary is the exact registry: its current definitions materialize nothing, and main provides no production preset query or collapse executor. Resource/organism/license filtering, unknown-license exclusion and recomputation of summaries after scoping belong to the separate API consumer and remain outside this build oracle; test-side fold SQL is not a production serving implementation. COSMOS external identifier translation is tested separately by the existing main tests; its bounded oracle intentionally exercises the supported offline fallback. These precise remaining gaps must be reported; registry equality does not establish serving-query equivalence.

The separate MetSigDB oracle invokes unchanged `load_resource` extraction/publication and the exact frozen source hashes. It creates no resource artifacts and does not invoke the full-build coordinator, resolver or external mappings. Its SQL fixture is optional focused verification; the running pinned release build has no additional mandatory verification phase.
