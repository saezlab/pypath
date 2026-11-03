# Reactions as typed entities, relations as projections

Status: implemented on `parquet-migration` (2026-10-16), not yet released;
see "Implementation" at the end for what landed, where, and what differs.

## Why

The explorer's transport filter shows 1,467 relations (921 after collapsing),
all from Mebocost, MRClinksDB, Cellinker and Macrophage. TCDB, Recon3D,
Human-GEM, Rhea and Reactome contribute none. Causes found on release
2026.10.15.1 (nicesrv):

- **The transport type was lost in the Biolink move.** `main` published a
  `Transport` entity type. Biolink has no such class, so every transport became
  a plain `molecular_activity`. What is left is the dataset name
  (`transport_reactions`) and the compartment on each participant. Reactome
  maps BioPAX `Transport` and `TransportWithBiochemicalReaction` to the same
  `'reaction'` as `BiochemicalReaction`.
- **Only Postgres derives transport pairs.** `_stage_transport_pairs` builds
  transporter→metabolite pairs from reactions, but only inside Postgres, for
  MetaLinksDB. The explorer reads Parquet, so it never sees them.
- **TCDB:** the cached `tcdb_substrates.tsv` on nicesrv is a saved HTTP 500
  page, so the build published 0 relations. The modelling itself is right.
- **Rhea parser bug:** Rhea lists each ChEBI once per reaction, and the parser
  pairs that list position by position with the equation's terms. A
  metabolite's second occurrence is lost, and the input/output roles and
  in/out sides shift. None of the 1,485 transport reactions keeps its cargo
  on both sides.
- **Reaction identity ignores compartments.** Recon3D's 4,197 transport rows
  land on 3,445 reaction entities. Once merged, a reaction's compartments can
  no longer be paired back to the right participants.
- **Every reaction is emitted twice.** Recon3D, Human-GEM and Rhea publish
  `reactions` and also `metabolic_reactions` + `transport_reactions`, which
  together make up the same set (Recon3D: 10,600 = 6,403 + 4,197).
- **The source's reaction type is dropped.** Recon3D ships an SBO term on every
  reaction (biochemical, transport, exchange, demand, sink, biomass). Human-GEM
  has a subsystem per reaction and TCDB ids for its transports. None is kept.
  The 1,809 exchange, demand, sink and biomass reactions of Recon3D are
  modelling devices, yet they sit in the build as ordinary reactions.

## Decisions

1. **A reaction is an entity.** It is a Biolink `molecular_activity`: inputs
   and outputs through `has_input` / `has_output`, catalysts through
   `enabled_by`, and stoichiometry, direction and compartment on each
   participant link. It stays the full statement of what the source says. A
   reaction is not turned into binary relations in place.
2. **Reactions carry a type.** The type is an SBO term in the entity's `type`
   slot (Biolink `rdf:type`). The source's own label (BioPAX class, Human-GEM
   subsystem) is kept verbatim as an annotation.
3. **Input modules emit each reaction once.** One `reactions` dataset per
   resource; the type replaces the `metabolic` / `transport` dataset split.
   Input modules emit no projected relations.
4. **Binary relations are projected downstream** from reaction entities by one
   shared build step. The first projection is transport: transporter
   `affects` cargo, aspect `transport`, mechanism `relocalization`, plus the
   from/to compartments. That is the shape TCDB, Cellinker, Mebocost and
   MRClinksDB already use.
5. **Reaction identity includes compartments.** The same reaction in different
   compartments, or the same transport across different membranes, are
   different entities.

## 1. The reaction record (what an input module emits)

| Part | Shape |
|---|---|
| Entity | `molecular_activity`; identifiers: the source's reaction id plus cross-references (Rhea, MetaNetX, KEGG, Reactome, BiGG/VMH, EC as topic) |
| Type | `type` = SBO term (table below); annotation `reaction_type_source` = the source's own label |
| Participants | membership `has_input` / `has_output`, each with `stoichiometry`, compartment, conversion direction |
| Compartment | one per participant link: Biolink `anatomical_context_qualifier` with a GO cellular-component CURIE (see below) |
| Catalysts | membership `enabled_by` → protein, complex, gene or gene-rule group |
| Regulators | unchanged (Reactome `controls`): relations, not part of this plan |

Catalysts today:

| Resource | Today | Change |
|---|---|---|
| Rhea | `enabled_by` membership to protein | none |
| Recon3D, Human-GEM | `associated_with` to gene or gene-rule group (AND clause) | decide: `enabled_by`, or keep `associated_with` and let the projection read it (see open questions) |
| Reactome | separate `catalyzes` relation (controller → reaction) in `controls` | the projection reads `catalyzes` too, or Catalysis becomes `enabled_by` membership |

**Compartments, unified.** Biolink's `anatomical_context_qualifier` is the
term: "a sub-cellular location where the relationship took place", with
ontology CURIEs as values. Biolink's own example is GO:0005794 (Golgi
apparatus), and its `cellular component` class is GO CC. Today there are two
terms (`biopax:cellularLocation` with names, `rhea:transport_side` with
in/out); both are replaced:

- Reactome: BioPAX `CellularLocationVocabulary` carries the GO xref. The parser
  keeps the id instead of the name only.
- Recon3D, Human-GEM: a fixed table from the model's compartments to GO, e.g.
  c → GO:0005829 cytosol, e → GO:0005576 extracellular region, m → GO:0005739
  mitochondrion, i → GO:0005758 mitochondrial intermembrane space,
  n → GO:0005634 nucleus, r → GO:0005783 ER, g → GO:0005794 Golgi,
  l → GO:0005764 lysosome, x → GO:0005777 peroxisome. Check the full list per
  model.
- Rhea is the one exception. `(in)` / `(out)` are relative sides of a membrane
  Rhea does not name, so there is no GO term to map them to. They stay a
  membrane-side attribute, and the projection reads it only where no GO
  compartment exists.

The projected transport relation needs two compartments, and Biolink has no
from/to qualifier, so it carries `transport_from` / `transport_to` with the
same GO CURIEs (Rhea: in/out).

A shared helper in `pypath/inputs_v2/` (e.g. `_reactions.py`) holds the SBO map,
the GO compartment map and the catalyst predicate, so the modules cannot drift.

Modules that emit reactions today: `recon3d`, `metatlas` (Human-GEM), `rhea`,
`reactome`, `kegg`. KEGG is audited in the same pass. TCDB, Cellinker,
Mebocost and MRClinksDB publish transport pairs directly and stay as they are.
PathwayCommons publishes pairs too, and gets its own change (section 4b).

## 2. Reaction types (SBO)

| SBO | Meaning | Recon3D | Human-GEM | Rhea | Reactome (BioPAX) |
|---|---|---|---|---|---|
| SBO:0000176 | biochemical reaction | native `sbo` (4,562) | other subsystems | equation not transport | BiochemicalReaction |
| SBO:0000185 | transport | native (4,229) | subsystem "Transport reactions" (+ TCDB id) | `_is_transport(equation)` | Transport |
| 176 + 185 | transport with conversion | — | — | — | TransportWithBiochemicalReaction (both types) |
| SBO:0000627 / 0000628 / 0000632 | exchange / demand / sink | native (1,560 / 145 / 101) | "Exchange/demand reactions" | — | — |
| SBO:0000629 | biomass production | native (3) | biomass reaction | — | — |
| SBO:0000177 | non-covalent binding | — | — | — | ComplexAssembly |
| SBO:0000179 | degradation | — | — | — | Degradation (today a `biological_process`) |
| SBO:0000183 / 0000184 | transcription / translation | — | — | — | TemplateReaction |

Check every non-Recon3D id against the SBO ontology before use. Add `sbo` to
`omnipath_core.naming.Namespace`, and SBO to the reference library's ontology
hub so types get labels and a hierarchy (e.g. `is_a` transport).

## 3. Reaction identity

The resolver anchors reactions on Rhea master ids and MetaNetX
(`identity_runtime.REACTION_HUBS`). Both are compartment-free, so
compartment-specific model reactions collapse into one entity.

- Find exactly where the merge happens: the resolver's reaction library, the
  entity key, or both.
- An entity is identified by its resolved reaction plus a compartment
  signature: the sorted (participant, side, compartment) set, or the source's
  own compartment-specific id where it has one (BiGG/VMH ids are).
- Keep the cross-resource link as the reference key, as for proteins
  (`reference_entity_key`). "All glucose transports, any membrane" still
  groups; "this transport" stays distinct.
- This is a resolution-rule change: update the entity-resolution living doc.

## 4. Projection step

**Where:** the Parquet build, per resource, after resolution (in or next to
`ParquetWriter`, before `_write_tables`). The projected relations land in that
resource's `relation.parquet`, credited to that resource, so the explorer,
the Postgres loader, the subsets and the client all read them with no logic of
their own. A release-level step was rejected: provenance is per resource, and
resource builds are independent.

**Grain:** per source reaction row (relation evidence `row_id`), not per merged
entity. That is correct today and stays correct after the identity fix.

**Transport rule (v1):**

- Reaction typed as transport (SBO:0000185 or a descendant), or untyped and
  the compartment rule holds. The rule: a participant appears as input and
  output of the same row, in different compartments.
- Cargo = each such participant, chemical entities only in v1. ATP/ADP never
  qualify because they are not on both sides; a co-transported H+ or Na+ does,
  and that is a real co-transport.
- Transporter = each catalyst as stated (protein, complex, gene, gene-rule
  group, Reactome set); nothing is expanded (decided 2026-10-16). No catalyst means no relation;
  the reaction alone stays.
- Emitted: transporter `affects` cargo with `object_aspect_qualifier=transport`
  and `causal_mechanism_qualifier=relocalization`. Annotations
  `transport_from` / `transport_to`, and a link to the reaction entity and its
  source row.
- Provenance: evidence `dataset = projection:transport`, so a projected
  relation can be told apart from a directly stated one, and the explorer can
  say "from reaction X".

**Candidates later, not v1:** enzyme→substrate/product pairs for
non-transport reactions (MetaLinksDB's "a conversion that moves nothing has no
pair reading"); protein cargo (nuclear import, receptor internalisation).

## 4a. Reactome specifics

- **Type:** the parser keeps the BioPAX class instead of mapping five classes
  to `'reaction'` (`parsers/reactome.py:920`, and the controlled side at
  `:1436`). ComplexAssembly that dissolves a complex is dissociation, not
  binding: decide by whether the complex is an output or an input.
  TemplateReaction is transcription or translation, by product type.
- **Degradation** is a `biological_process` today. Becoming a typed
  `molecular_activity` changes what Reactome's degradation controls point at;
  decide in the migration.
- **Compartments:** the parser joins compartment names with `; `. It keeps the
  GO id of the BioPAX `CellularLocationVocabulary` instead.
- **Catalysts** stay where they are: a `catalyzes` relation (controller →
  reaction) in the `controls` dataset, with its own Reactome Catalysis id and
  evidence. The projection reads `catalyzes` as well as `enabled_by`.
- **Pairing by source reaction, not by entity.** 887 Reactome reactions
  resolve to a Rhea id, and 197 entities hold several Reactome reactions:
  16,022 rows → 15,533 entities. RHEA:21292 (PI3K, PIP2 → PIP3) holds 50, one
  per signalling context, each with its own catalyst complex. The
  compartment signature does not split them (same compartments, same
  chemistry), and that is right. But catalysis lives in another dataset, so
  the projection must join a catalyst to its reaction by the Reactome stable
  id of the controlled reaction. Joining by the merged entity would pair one
  context's catalyst with another context's cargo.
- **Expected transport:** chemical cargo from ≤416 reactions, 315 with a
  catalyst (169 complexes, 153 proteins, 87 generic sets). About 520 more
  reactions move a protein or complex (question 4).

## 4b. PathwayCommons transport pairs

PathwayCommons has no reaction entities to project from: its extended SIF is
pairs only. Its `controls-transport-of-chemical` type ("the protein controls
a reaction that changes the cellular location of the small molecule") is
mapped to aspect `localization` today: 3,942 lines → 3,695 relations. The
`MEDIATOR_IDS` column names the interactions behind each line, so the input
module can tell direct transport from the rest:

| Mediator | Lines | Reading | Model as |
|---|---|---|---|
| a catalysis of the transport (`reactome:Catalysis`, `humancyc:Catalysis`, `pathbank:TransportCatalysis`, `reconx:control_`) | 3,328 (84%) | the protein is the transporter, or a subunit of it; top cargos Ca²⁺, K⁺, Na⁺, Cl⁻ | transport, the same shape as the projected relations |
| regulation only (`reactome:Control`) | 154 | the protein regulates a transport | `localization`, as now |
| neither | 460 | the protein moves *with* the cargo in a complex (LDLRAP1 in LDL clearance, NUP50 in snRNP assembly), or indirect PID signalling (JAK2 → Ca²⁺) | `localization`, as now |

A line with both a catalysis and a regulation mediator counts as transport:
one direct path is enough. Keep the mediator ids as provenance annotations;
the Reactome `R-HSA-` ids link a pair to its reaction entity.

Three differences from the projected relations, none fixable from the SIF:

- **The subject is a member protein.** PathwayCommons splits a complex
  catalyst into its proteins, so each subunit becomes a transporter.
- **No from/to compartments.**
- **Overlap with direct inputs.** Most lines come from Reactome (2,469) and
  Recon X (1,036, an older Recon), which we also ingest directly with more
  detail. HumanCyc (710), PathBank (282), PID (123), PANTHER and INOH are
  sources we have no direct input for.

The protein→protein `controls-transport-of` (8,562) stays `localization`: its
cargo is a protein, which is question 5 below.

## 5. Downstream simplification

- Postgres: remove `_stage_transport_pairs` and the `transport_event_type`
  capability (`relational/pipeline.py`). Check that `_stage_reaction_hyperedges`
  still reads one reaction per row once the dataset split is gone.
- MetaLinksDB preset (`network_views/_definitions.py`): the reaction component
  filters on aspect `transport` like the rest. Rewrite its docstring numbers.
- COSMOS: replace `TRANSPORT_ENTITY_TYPE = "omnipath:unavailable_transport_type"`
  with the SBO transport type.
- Delete the `Transport` entity-type sections of the retired
  `network_views/sql/metalinksdb.sql`.
- Explorer: the transport filter works unchanged. Add a reaction-type facet
  for molecular activities.

## 6. Independent fixes (can land first)

- TCDB: restore the good `tcdb_substrates.tsv` on nicesrv (copy under
  `pypath-cache/_variants/6c565fddb6c89186/tcdb/`), rebuild TCDB.
- Build guard: fail a resource when a download is an HTML error page, or when
  a dataset yields 0 records from a non-empty file.
- Rhea: expand the per-reaction ChEBI list back to equation order (each
  stripped term → its first-appearance index), so every term gets its own
  ChEBI, role and side.
- Water (InChIKey XLYOFNOQVPJJNP) displays as "FA 18:5;OH": fix the label
  choice.

## Order of work

1. Independent fixes (section 6).
2. SBO namespace and ontology, `_reactions.py` helper, the reaction-record
   contract with tests.
3. Migrate `recon3d`, `metatlas`, `rhea`, `reactome` to one typed `reactions`
   dataset with GO compartments; audit `kegg`; PathwayCommons mediator split
   (section 4b).
4. Reaction identity with compartments; living doc.
5. Projection step with tests (fixtures: a Recon3D transport, a Rhea
   `X(out) = X(in)` symport, a Reactome Transport with a complex catalyst, a
   metabolic reaction that must yield nothing).
6. Rebuild the affected resources on nicesrv and compare with the baseline.
7. Downstream removals (section 5) and a Postgres reload on dev5.

## Baseline (release 2026.10.15.1)

| Resource | Transport relations now | Expected after |
|---|---|---|
| TCDB | 0 | substrate file restored: some thousands (9,124 TC→substrate lines, expanded to proteins) |
| Recon3D | 0 | ~3,350 catalyst–cargo pairs (per-row count today) |
| Human-GEM | 0 | ~4,430 |
| Rhea | 0 | non-zero once the parser is fixed (1,485 transport reactions) |
| Reactome | 0 | from ≤416 chemical-moving reactions, 315 with a catalyst |
| PathwayCommons | 0 (3,695 as `localization`) | ~3,100 transport relations from the 3,328 catalysis lines; ~600 stay `localization` |
| Mebocost, MRClinksDB, Cellinker, Macrophage | 1,467 | unchanged |

## Open questions (for colleagues)

1. Do reactions without a named transporter count as transport anywhere, or
   are they reaction-only? (About 2 in 5 in Recon3D and Human-GEM.)
2. Model pseudo-reactions (exchange, demand, sink, biomass): drop them, or keep
   them typed and hidden by default?
3. GPR genes: `enabled_by` (Biolink's range allows genes), or keep
   `associated_with` because a gene rule is not a physical catalyst?
4. Protein cargo (nuclear import, internalisation): a projected relation with
   aspect `localization`, or nothing?

Decided 2026-10-16:

- Compartments: `anatomical_context_qualifier` with GO cellular-component
  CURIEs; Rhea's in/out stays a membrane side.
- Transporter subject: as stated, no expansion.
- PathwayCommons: catalysis-mediated `controls-transport-of-chemical` lines
  become transport; the rest stay `localization` (section 4b).

## Implementation (2026-10-16)

**Where it lives**

| Piece | Code |
|---|---|
| Vocabulary: SBO types, GO compartments, compartment identity | `packages/core/omnipath_core/reactions.py`, new terms in `source_attributes.py` |
| Builders the reaction modules share | `pypath/pypath/inputs_v2/_reactions.py` |
| Recon3D, Human-GEM, Rhea, Reactome: one typed `reactions` dataset, GO compartment per participant, compartment set on the reaction | the four modules and their parsers |
| PathwayCommons: catalysis-mediated `controls-transport-of-chemical` lines as transport, mediator ids kept | `pypath/pypath/inputs_v2/pathwaycommons.py` |
| Compartments in a reaction's entity key, reference key unchanged | `omnipath_build/silver.py`, `writer.py` (`entity_map`), `extract/observations.py` |
| Transport projection, per source row, catalysts by row or by `omnipath:source_reaction` | `omnipath_build/projections.py`, called in `ParquetWriter.close()` after complexes |
| Postgres: `_stage_transport_pairs` removed, `transport_event_type` available | `relational/db/derived_tables.py`, `relational/pipeline.py` |
| COSMOS: transport by SBO type, else by movement | `cosmos/build.py`, `cosmos/sql/project_edges.sql` |
| Explorer labels for the new terms | `packages/web/src/lib/utils/annotation-presentation.ts` |
| Download and empty-dataset guards | `pypath/share/downloads.py`, `Dataset(allow_empty=…)`, `two_phase.py`; replay allows empty |
| Tests | `packages/build/tests/test_reaction_projection.py`, `test_two_phase.py` |

**Decisions taken while implementing**

- The SBO type does not gate the projection; the movement does. A reaction a
  source types as biochemical that still moves a cargo (623 Human-GEM rows,
  e.g. the carnitine shuttle) states a transport too.
- Compartment names stay beside the GO qualifier (`biopax:cellularLocation`),
  for display and for the Postgres party compartment, which reads them.
- The reaction's compartments on a relation endpoint (Reactome's controls)
  are its identity only; they are not copied onto the relation as endpoint
  annotations, where a qualifier is not allowed.
- Open questions 1–4 took the conservative default: no transporter → reaction
  only; pseudo-reactions kept, typed; GPR stays `associated_with` (the
  projection reads it); no protein cargo.
- Fixed on the way: Rhea's fallback UniProt pattern skipped every O/P/Q
  accession. Recon3D `metabolites` and CellChat `zebrafish_metabolite_genes`
  are empty by design and declare `allow_empty`.

**Not done, and why**

- No `sbo` identifier namespace and no SBO ontology in the reference library:
  types are annotation values; the explorer shows the CURIE. Follow-up.
- Explorer reaction-type facet: an entity-level facet through engine, API and
  web. Follow-up.
- Retired `network_views/sql/metalinksdb.sql` unchanged: it is a checksummed
  parity reference for `main`, executed by nothing.
- KEGG reactions stay untyped: KEGG states no compartments or type, and
  typing them all biochemical would add a claim.
- Water's "FA 18:5;OH" label: fixed separately the same day, see below.

**Results** (version 2026.10.16 built on nicesrv into
`/root/data/scratch/reaction-model-20261016/out`, not the served data root;
compared with release 2026.10.15.1)

| Resource | Transport relations before → after | Reaction entities before → after | Notes |
|---|---|---|---|
| TCDB | 0 → 19,908 | — | substrate file restored; error page kept as `tcdb_substrates.tsv.error.html` |
| Rhea | 0 → 27,294 | 18,624 → 18,624 | all transporters UniProt proteins, every organism; participant fallback on 135 polymer reactions (0.7%), no transport among them |
| Recon3D | 0 → 3,334 | 7,748 → 10,545 | 10,600 model reactions; types 4,534 biochemical, 4,216 transport, 1,547 exchange, 145 demand, 100 sink, 3 biomass |
| Human-GEM | 0 → 4,412 | 9,778 → 12,709 | 12,931 model reactions; types 6,984 biochemical, 4,075 transport, 1,652 exchange |
| PathwayCommons | 0 → 3,131 | — | the rest of `controls-transport-of-chemical` stays localization |
| Reactome | 0 → 608 | 15,533 → 15,682 | transporters 248 proteins, 216 complexes, 144 sets; all 6,733 catalyses still reach a reaction with participants; Reactome's BioPAX types nearly every conversion BiochemicalReaction, so transport comes from movement |

Pitfall met: Reactome's parser caches its parse (`_CACHE_VERSION`); a parser
change without a version bump reuses the old rows.

## Follow-up: lipid-name labels on non-lipids (2026-10-16)

Water was labelled "FA 18:5;OH" because ChEBI lists "HOH" as a synonym and
Goslin parses it as that hydroxy fatty acid; the label rule put any Goslin name
before ChEBI's. Glutamine's "FA 5:0;2NH2,5NH2;5oxo" is a true parse of its
systematic name ("diamino-oxopentanoic acid"): Goslin reads any substituted acid
as a fatty acyl, and an "-oate" anion as a wax ester without its alcohol
("hexanoate" → "WE 6:0/0:0", in SwissLipids too).

- **Label rule** (`omnipath_resolver.identity_runtime`, record policy v5): a
  fatty-acyl shorthand labels only from a SwissLipids or LIPID MAPS record; a
  `WE …/0:0` shorthand never. On the 2026-10-13 library, 0.46% of the 318,732
  entities with a lipid-name label change (sample of 5,000), all FA or WE:
  water, glutamine, perfluorooctanoic acid, abscisic acid, crotonate, resolvins.
  A few real fatty acids no lipid database holds fall back to their systematic
  name.
- **Hub index** (`omnipath_build.identity.goslin`): `WE …/0:0` is no claim in
  any hub; in ChEBI, HMDB and RefMet a parse of a name not written in lipid
  notation must have the record's composition (formula or InChI, hydrogen
  aside). Records with a lipid name: ChEBI 20,649 → 19,271, HMDB 168,512 →
  168,463, SwissLipids 777,126 → 776,910; LIPID MAPS and RefMet unchanged.
  Takes effect with the next hub-index and library build.
