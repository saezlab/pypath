# Current PostgreSQL relational backend

This namespace reuses the downstream code at omnipath-build commit
`9f9bb709c764`. The independent test oracle retains that source unchanged.
The stable identifier registry is rendered from its pypath gitlink
`33f37fbaab59993d24f5c32bb4b3e7587085bcf4`; its exact source blob is recorded in
`identifier_types.py`. Migration pypath need not expose newer enum members.

The published Parquets remain unchanged. The canonical loader projects their
resolved identities into main's relational substrate before these steps run.
No input parser, raw-record resolver, resource build or identifier download is
called by the shared derivation pipeline.

## Preserved main contracts

* Canonical entity/relation tables, source partitions, typed foreign keys,
  identifier links, evidence/annotation links, ontology links and main indexes.
* Main interaction content/record UUID expressions, assertion grain, tier
  precedence, source aggregation, nullable direction/signs, reaction chemical
  grouping, party ordering and the existing pair/reaction collision policy.
* Main query/count tables, ontology derivation, chemical classifications,
  label ranking/cascades, bitmap/facet tables, overlap summary, resource/license
  metadata and build-manifest shapes.
* Main MetSigDB membership table and indexes, all five source rules, bounded
  hierarchy walks and publication upsert shape.
* Main network preset registry and compositions. The old materialized-view SQL
  remains as reference code; main already treats those views as unmanaged and
  does not recreate them in a fresh build.
* Main COSMOS table, indexes, projection order, grouping and summary shape.

`run_main_derivations` follows main's shared-step order, seeding identifier
authority first as main's ingestion does. Resource descriptors come from pinned
release metadata. Missing `mints` declarations use the six exact declarations
from main's pinned pypath source; source commit and file hashes are recorded in
`resource_declarations.py`. No authority is inferred from data. Licenses use main's curated metadata precedence and unknown rules.

`omnipath_subsets.scientific.run_product` defers internal helper commits. The shared subsets runner commits each
complete product's tables, indexes and metadata together. Rollback remains real.
This intentionally replaces main's intermediate source commits with the user's
chosen complete-product checkpoints; an unfinished product cannot overwrite a
previously committed product.

## Vocabulary and published-identity adaptations

* Entity readers use `protein`, `gene`, `chemical_entity`, `small_molecule`,
  `ontology_class`, `pathway` and `molecular_activity`. Chemical gates include
  both published chemical families rather than assuming earlier normalization.
* An explicit activity subject with `has_input`/`has_output` supplies main
  reactant/product roles. Pathway membership does not become a reaction star.
  Stoichiometry and compartment come from the same source-owned membership
  evidence before main's aggregation. Separate input/output graph triples join
  by activity, member and source when deriving movement of the same cargo.
* `enabled_by` remains activity-to-enzyme in the graph, but its derived main
  interaction endpoint order is enzyme-to-activity. Explicit `catalyzes` already
  has that order. Generic regulation and gene association do not imply enzymes. Source-owned
  miRBase `mirbase_mature derives_from mirbase_precursor` becomes the same main
  precursor-to-mature maturation fact; other `derives_from` claims stay unchanged.
* Object direction qualifiers supply positive/negative flags. Absent opposite
  flags remain NULL. Published serving `sign=0` and symmetry flags are not
  treated as explicit biological assertions.
* Context-specific ConnectomeDB subject/object ligand/receptor annotations
  retain main's role gates. Direction alone does not invent pharmacological
  roles or orthosteric/allosteric classification.
* Publications populate main PMID/DOI fields after stripping their explicit
  prefixes. Quantity companions preserve units, prefixes, source fields and
  comparators. Only Ki/Kd/IC50/EC50 BAO endpoints populate affinity, with compatible
  concentrations normalized to nM; bounded or unknown-unit values remain in
  relational annotations and companions instead of becoming misleading scalars.
  pChEMBL requires the original `pchembl_value` quantity source field. Main's
  score retains STITCH's `combined_score` meaning; action score and ChEMBL assay
  confidence remain separate annotations/quantities. Mechanism
  curation requires ChEMBL's explicit `mechanisms` dataset.
* MetSigDB uses published chemicals, as approved by the user. Its ChEBI display
  column strips a recognized numeric `CHEBI:` prefix to retain main's bare-ID
  output; dictionary values and canonical/published identities remain exact. Status `5/published`
  states that upstream resolution is already complete and the original resolver
  outcome was not retained; it is not a fabricated `matched` status. All main
  resolution status/reason columns and foreign keys remain present.
* Primary gene symbols use main's attestation/shortest/alphabetical ranking for
  genes and published proteins. Proteins carry `published_gene_symbol` label
  provenance without gene recanonicalization.
* Other non-chemical entities without a symbol (ontology terms, pathways,
  complexes, reactions, miRNAs) take their most attested published name, ranked
  the same way (`name` rule), before the identifier fallback.
* Source-owned known taxonomy can populate interaction facts when a shared
  canonical identity has conflicting organism occurrences; ambiguous values are
  not guessed.
* COSMOS reads directions from available event and source-owned membership
  annotations. Oriented `right_to_left` input is not swapped a second time.
  Labels use published UniProt/ChEBI aliases; online translation is disabled.

## Information boundaries

Main's unsupported diagnostic tables/columns remain present, and unavailable
capabilities are recorded in both `build_manifest` and `build_capability`.
Serving Parquets aggregate entities, so original entity-occurrence identity,
original resolver flags/conflicts and old gene-anchored identities cannot be
recovered exactly. Original metabolic-domain source-priority fields and precise
pharmacological roles are absent. Some ligand/receptor sources retain only a
generic predicate, so their original occurrence-role flags cannot be guessed.
ChEMBL activity claims may lose a mechanism identifier even though explicit
mechanism datasets are identifiable. IntAct strips the original confidence
value's prefixed spelling, so the old numeric-filter outcome is unavailable. Biolink activities do not retain the old
Reaction/Transport type distinction; explicit same-cargo compartment movement
can establish transport, other events cannot be guessed. Ontology definitions
and comments share `description`, so their original distinction is unavailable.
External parsed-lipid mapping tables are not fetched again.

These are earlier-input information boundaries, not reasons to change the
Parquet schema or replace main's downstream physical model.

## Runtime requirements and checks

The code retains psycopg2 semantics (`psycopg2-binary`) and PyYAML. PostgreSQL must
provide main's roaringbitmap and pg_trgm extensions. Package SQL/YAML files must
be included in the wheel.

Focused tests cover independent frozen-main DDL/index/UUID/identifier contracts,
shared orchestration, deferred commits and cancellation. Opt-in bounded SQL
fixtures compare full main interaction results and protein/gene label ranking
inside a guarded private PostgreSQL container. They do not build resources.

Extensions are shared in `public` so new release schemas reuse the same main
bitmap type and trigram operators. New installs explicitly target `public`; an
existing extension elsewhere raises a setup error and is never silently moved
or replaced by a build. The miRBase namespace/source maturation adaptation is
assertion-local: a shared canonical `derives_from` triple produces separate
main maturation and generic facts with only their own contributing resources.
