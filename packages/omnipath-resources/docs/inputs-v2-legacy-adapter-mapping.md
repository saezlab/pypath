# Mapping legacy complexes, enzyme–substrate records and networks to inputs_v2

Design proposal, 2026-09-05. This document specifies adapters; it does not claim that they have been implemented or validated against live datasets. Scope and resource entry points are in the [inventory](inputs-v2-legacy-adapter-inventory.md).

## 1. Architecture and boundaries

Use three explicit output contracts, sharing a legacy-function runner and participant conversion:

```text
legacy function + pinned arguments
    → output normalizer → serializable parsed rows
    → EntityBuilder / RelationBuilder → silver extraction → existing build
```

Construct ordinary `Dataset(download=None, raw_parser=..., mapper=...)` objects and expose them through a module-level `Resource`. The [existing Dataset implementation](../pypath/inputs_v2/base.py) already supports this and [miRBase](../pypath/inputs_v2/mirbase.py) demonstrates a legacy-backed parser. No new dataset subclass is required. Resource-specific modules should contain metadata, the chosen function/arguments and reviewed semantic mappings; common execution code belongs in a shared adapter module. Proposed API names below are illustrative, not existing interfaces.

| Shared component | Responsibility |
| --- | --- |
| Legacy runner | Lazily resolve an allowlisted callable; bind configured input arguments; translate refresh/species options; report errors with resource/dataset/function context. |
| Complex normalizer | Accept dictionary values or an iterable of `Complex`; convert components, IDs, evidence and attributes into plain serializable rows. |
| Enzyme–substrate normalizer | Accept configured dictionary/namedtuple fields; produce one assertion per explicitly supported enzyme–substrate–site observation. |
| Network normalizer | Apply a supported `NetworkInput` contract and explicit semantic overrides; produce endpoint, predicate/qualifier and evidence fields. |
| Participant mapper | Map explicit identifiers/types/taxonomy and recursively convert actual complex participants. Preserve unknown or ambiguous identity instead of guessing. |

Do not auto-detect a biological schema from a function suffix or tuple length. Output-shape checks validate a declared contract; they do not choose it. A list of proteins can mean alternative identifiers, candidate enzymes, a family or a physical complex depending on the source.

### Invocation and cache contract

- Keep adapter execution options separate from legacy function arguments. Discovery passes `source` and `dataset`; `Dataset.raw` passes `force_refresh`. Consume those at the adapter boundary. Unknown legacy arguments should produce a clear error, not be silently filtered out.
- Normalize the requested taxon once, then explicitly translate to the legacy parameter (`organism`, `ncbi_tax_id`, numeric ID, species name). Human-only inputs must reject unsupported requests rather than label human rows with the requested species.
- Invoke each legacy function once per dataset execution. Iterate where the legacy API permits it; wrapping a function returning a materialized list does not make its download/parser streaming.
- `force_refresh=True` must refresh all applicable legacy download and intermediate caches. A Curl cache context alone does not invalidate bespoke pickle/SQL caches. Each registration must declare its refresh strategy; unsupported refresh should fail clearly rather than pretend to succeed. If using global Curl contexts, keep the context around full generator consumption and isolate concurrent executions.
- Convert namedtuples/sets/complexes to deterministic plain structures before parsed-cache/evidence serialization. Preserve sequence order, align member counts, sort unordered collections deterministically, and never serialize arbitrary legacy objects using `str()` as a substitute for their fields.
- Include function, normalized arguments, adapter version and semantic mapping version in cache/build identity. Two hu.MAP releases or differently filtered network datasets must not share an indistinguishable cache.
- Retain the original legacy-returned record alongside normalized values, or a lossless serializable projection of it. This is **legacy output**, not necessarily the original downloaded source row. Verify that the actual build persists this payload; returning extra dictionary keys alone is not proof of evidence retention.

### Shared schema conventions

Follow [BIOLINK_MODELING.md](../pypath/inputs_v2/BIOLINK_MODELING.md) and the current migrated modules. Biological types/predicates/qualifiers use the pinned native Biolink model; identifier systems use `Namespace`. Names and synonyms are identifiers (`Namespace.NAME`, `Namespace.SYNONYM`) under the current application convention. Place real accessions ahead of name aliases. Do not recreate the older convention of name/synonym annotations.

Taxonomy and publications use valid `NCBITaxon:` and `PMID:` CURIEs. Preserve PMC/DOI references using their actual identifier systems; arbitrary reference strings are not PMIDs. Resource metadata remains in `ResourceConfig`; primary-source and via-source provenance remain distinguishable. No universal “human”, “UniProt”, “phosphorylation” or positive-effect default is allowed across all inputs.

## 2. Complex adapter

### Accepted input and normalized row

The principal contract is `dict[str, intera.Complex]`; scConnect instead returns a set of `Complex`. Iterate dictionary **values**, not composition-string keys. ComplexPortal's optional `(complexes, details)` result requires an explicit selector/join; never infer that any tuple contains complex records.

The normalized row should include `members` (identifier, namespace, count and count provenance), `identifiers` (namespace/accession pairs), `name`, `taxon`, `publications`, `sources`, source-specific `attributes`, and the original legacy key when available. Keep member records together until projection to aligned builder columns.

| Legacy field | New representation | Rule / limitation |
| --- | --- | --- |
| Complex object | `EntityBuilder(entity_type=MacromolecularComplex)` | Only for an actual source complex. Do not promote a generic group to this class. |
| `components` keys | Nested `Protein` members with actual accession namespaces | UniProt is the conventional contract for these legacy complex functions, but confirm it per registration. Reuse the same participant conversion for complexes inside network rows. |
| Membership | `MembershipBuilder` and `Member` / `MembersFromList`, predicate `has_member` | Keep the complex intact; no automatic all-pairs interaction expansion. |
| `components` values | Membership annotation using native `slots.stoichiometry`, when supported by source evidence | A default count of 1 often means membership was listed once, not experimentally known 1:1 stoichiometry. Retain unqualified counts in raw payload; publish them as stoichiometry only under a reviewed source contract. |
| `ids: source → set[accession]` | `IdentifiersBuilder` with an explicit source-to-`Namespace` map | Source labels are not automatically namespaces. Unknown ID systems require registry work or raw retention. Preserve several IDs without inventing cross-resource equivalence beyond the returned object. |
| `name` | `Namespace.NAME` identifier | Do not mint a resource accession from a display name or composition string. |
| `ncbi_tax_id` | `slots.in_taxon` on the complex and members when justified | The legacy class defaults to human and assumes one organism. Check that the input actually propagated the requested species. Never treat that default as independent source evidence. |
| `references` | `slots.publications` | Normalize and validate individual references; preserve original forms in evidence. |
| `sources` | Resource/via provenance; `supporting_data_source` where appropriate | Distinguish the adapter's acquisition source from databases cited by an aggregator. |
| `attrs` | Explicit native/external attribute mappings or retained parsed payload | hu.MAP confidence and Compleat prediction status need source-specific meaning. Do not convert structured data into description strings. |
| `interactions` | Separate relations only after auditing the source's internal interaction contract | Complex co-membership alone does not establish a direct physical interaction. |

The current membership builder supports membership annotations, so stoichiometry does not require a new `Membership` field. Verify its integer value survives silver extraction and serving projection. Do not repeat a member entity N times to encode stoichiometry: deduplication could erase the count.

For complexes without source accessions (for example hu.MAP1), preserve structural membership and let the existing identity/resolution path handle accessionless complexes. Do not invent a `hu.MAP` accession from `COMPLEX:...`. Test that the extractor/resolver retains such entities; if it does not, that is an explicit implementation dependency before enabling the resource.

### Per-resource configuration

| Resources | Required choices |
| --- | --- |
| ComplexPortal | Pin species and `return_details`; audit loss from discarded non-UniProt participants, stripped isoforms and composition-based merges. The adapter cannot recover these from the resulting objects. |
| Compleat | Pin `predicted`; retain prediction/method evidence if it survives the function output; review source IDs and stoichiometry semantics. |
| Havugimana | Confirm identifier mapping and species; preserve source evidence actually returned. |
| hu.MAP1/2/3 | Separate versioned datasets; explicitly pin `humap_version` for 2/3 and `min_confidence`; do not assume source IDs exist. |
| PDB | Preserve PDB identifiers and supported component counts; check organism handling and ambiguity across structures. |
| SPIKE | Pin `min_confidence`; preserve complex identity and source-specific evidence. |
| scConnect | Accept set output; pin organism; audit whether annotation-derived grouping is a physical complex assertion. |
| HPMR / 3DComplex | HPMR documents empty/ambiguous complex output; 3DComplex's complex function raises `NotImplementedError`. Keep disabled until source-specific work resolves these. |

## 3. Enzyme–substrate adapter

### Row contract and expansion

Normalize `kinase`/`enzyme`, substrate identifiers, enzyme/substrate isoforms, `typ`/`ptm_type`, `resaa`, `resnum`, `instance`, `start`, `end`, references, supporting databases and source scores. Rename `kinase` internally to `enzyme`: DEPOD enzymes are phosphatases.

If the source supplies a list of candidate enzymes, expansion must preserve its candidate/inference status. Do not reinterpret that list as a physical complex or assert every candidate as experimentally established. Alternative identifiers for one substrate belong to that participant; separate substrate candidates require explicit expansion semantics.

| Legacy field | New representation | Rule / limitation |
| --- | --- | --- |
| Enzyme and substrate IDs | Subject and object entity builders in `RelationBuilder` | Use the registry's `id_type_enzyme` / `id_type_substrate`, including alternate field names. Keep RefSeq protein and nucleotide namespaces distinct. Resolve uncertain identity through the existing resolver; do not claim a nucleotide accession directly identifies a protein without an explicit translation. |
| Modification type | `affects` plus reviewed `object_aspect_qualifier`, effect direction and optional causal mechanism | Phosphorylation increases the phosphorylation aspect; dephosphorylation decreases that aspect. Neither implies increased/decreased protein activity. Other PTMs need an exact crosswalk. |
| `resaa`, `resnum` | Site-specific assertion context, retained together with substrate and isoform | Never attach the position as a general attribute of the substrate protein. The current `Relation` has no dedicated residue/site fields; see the staged site contract below. |
| `substrate_isoform`, `enzyme_isoform` | Explicit isoform identity when the namespace supports it; otherwise preserved assertion context | A site coordinate must remain tied to its source sequence. Do not silently map an isoform coordinate onto the canonical sequence. |
| `instance`, `start`, `end` | Site/motif evidence payload | A peptide motif is not the full `has_biological_sequence` of the protein. Preserve coordinate convention; do not silently shift offsets. |
| `references` | `slots.publications` on the relation | A publication count (`npmid` in MIMP) is not a publication accession. |
| `databases` | Supporting-source and via-source evidence | Preserve aggregator provenance, especially ProtMapper and MIMP. |
| `score`, `npmid`, other attributes | Reviewed source-specific quantitative attributes or raw evidence | Keep score scale/type and field origin; do not label a source score as a probability without evidence. |

Concrete proposed projection:

```text
enzyme E, substrate S, typ=dephosphorylation, site=S42, isoform=2
→ E —affects→ S (isoform 2)
  object_aspect_qualifier = phosphorylation
  object_direction_qualifier = decreased
  original site observation retained as evidence
```

`S42` here means residue serine at position 42, not a new protein identifier. This protein-level relation is a projection of a site observation; do not describe it as a fully modeled site graph.

### Site representation: implementation gate

The current silver `Relation` supports annotations but has no dedicated structured site slot. Initial mapping can emit the qualified protein relation **only if** the exact residue, position, isoform, motif and source record survive as separate evidence occurrences. Two sites with the same endpoints must not be lost through pair-level deduplication.

Exposing sites as serving attributes or first-class sequence features requires a separately validated schema mapping: select published terms/classes, encode the site-to-substrate/isoform link, and test extractor, storage, identity and API behavior. Do not invent a local `resnum` ontology predicate or misuse a genomic coordinate slot for protein positions. Until this is done, report “protein relations with retained site evidence”, not lossless site-level serving support. If the build cannot preserve those occurrences, block enzyme–substrate enablement rather than drop sites.

### Registered resources and exceptions

| Resource | Required overrides / audit |
| --- | --- |
| phosphoELM | UniProt endpoints; translate species argument explicitly; pin `ltp_only`; preserve experiment evidence and sites. |
| dbPTM | Enzyme gene symbol, substrate UniProt; exact `typ` crosswalk; multiple candidate enzymes and organism filtering. |
| HPRD | Enzyme gene symbol; use `substrate_refseqp` with RefSeq protein namespace instead of blindly reading `substrate`; human registry context. |
| Li2012 | Gene-symbol endpoints; phosphorylation may be a reviewed dataset constant; preserve source residue/position. |
| DEPOD | UniProt endpoints; dephosphorylation decreases phosphorylation, not necessarily activity. |
| PhosphoSite | Pin `raw=True`, species translation and `strict=True` for source-species observations. `raw=False` returns legacy objects. The old core processor defaults `strict=False`, so parity with that inferred output is not automatic. |
| PhosphoNetworks | Gene-symbol endpoints; preserve score and site; phosphorylation only as a reviewed dataset constant. |
| MIMP | Gene-symbol enzyme; substrate symbol and alternate `substrate_refseq` nucleotide ID require explicit interpretation. Kinase family expansion has already occurred upstream; preserve inference and publication-count evidence. |
| ProtMapper | UniProt endpoints; `interactions=False`; pin `only_evidences` / `only_literature`; preserve supporting databases and individual evidence. A phosphorylation default must be justified for the selected export. |
| KEA | UniProt endpoints; explicit phosphorylation records; preserve source PMID and site. |
| iPTMnet | Namedtuple aliases: `enzyme→enzyme`, `ptm_type→modification`; preserve both isoforms, site, score and references. |

The old [EnzymeSubstrateProcessor](../pypath/core/enz_sub.py) also maps IDs, validates sequences/sites, searches isoforms and may infer orthology. This adapter does not silently invoke that entire processor. Distinguish source observations from separately configured normalization/inference, recording validation status. Never mark a wrapped source row sequence-validated merely because the old full pipeline used to validate it.

## 4. Interaction / NetworkInput adapter

### Translate an explicit supported subset

Reuse [NetworkInput](../pypath/internals/input_formats.py) and [data_formats.py](../pypath/resources/data_formats.py) as extraction metadata. Do not invoke the legacy `Network` constructor simply to obtain normalized rows: it adds mapping, filtering, aggregation and potentially complex expansion beyond this adapter's contract.

| NetworkInput field | Adapter behavior / schema mapping |
| --- | --- |
| `input`, `input_args` | Resolve allowlisted callable and explicit arguments. Stale functions are configuration errors. Direct paths/URLs remain a different download/parser route. |
| `separator`, `header` | Needed for text rows only; materialized record outputs must not lose their first row because file-header rules were applied twice. |
| `id_col_a`, `id_col_b` | Extract endpoint values; validate required columns and distinguish missing endpoints from unsupported shapes. |
| `id_type_a/b` | Explicit namespace crosswalk, including supported per-row definitions. No universal UniProt conversion. |
| `entity_type_a/b` | Explicit Biolink type crosswalk. RNA, genes, proteins, chemicals and complexes retain distinct domains. Unknown types are reported. |
| `taxon_a/b`, `ncbi_tax_id`, `only_default_organism` | Implement documented constant/column/dictionary taxon forms and selection rules; preserve endpoint-specific species. Do not attach organism taxonomy to chemicals. |
| `is_directed` | Decode supported boolean/column/value/separator forms into the normalized row. Orientation and Biolink predicate symmetry need independent semantic review. |
| `sign` | Decode positive/negative values exactly. Translate causal effects to direction qualifiers only with a reviewed predicate/aspect mapping. Unknown and conflicting signs must remain explicit. |
| `references` (`refs` on object), `must_have_references` | Normalize configured references and apply explicit inclusion policy with counts. Do not silently inherit process-global settings. |
| `positive_filters`, `negative_filters` | Implement the selected legacy forms faithfully and count exclusions. Callable behavior needs explicit audit/tests; unsupported forms fail during configuration validation. |
| `resource` | Decode constant/source column and delimiters; preserve primary and supporting/via sources. |
| `extra_edge_attrs` | Review field-to-slot/qualifier/measurement mappings individually; retain unmapped structured fields in evidence. |
| `extra_node_attrs_a/b` | Apply only source-supported participant attributes. Context-specific values may belong to the assertion instead. |
| `interaction_type`, `data_model`, `dataset` | Hints for reviewed mapping selection and provenance; not automatically valid Biolink predicates or categories. |
| `expand_complexes` | Initial adapter preserves actual complex participants; a configuration expecting expansion is an explicit behavior difference that must be rejected or overridden in the registration. |
| `allow_loops`, `unique_fields` | Explicit loop policy and evidence-aware identity/deduplication; preserve source observations with different sites, signs or evidence. |
| `mark_source`, `mark_target` | Review role annotations; a legacy boolean marker is not automatically an ontology term. |
| `curl_args`, `huge` | Execution metadata; do not forward as biological annotations or silently assume the legacy callable accepts them. |

Validate configuration before downloading: every nontrivial field must be implemented, explicitly overridden with a documented reason, or rejected. Start with constant endpoint namespaces/types, positional/namedtuple rows, standard reference delimiters and documented sign/direction/filter forms. Complex callbacks and multi-return APIs remain opt-in. Do not copy suspected bugs or accidental global behavior just to claim compatibility.

### Biological mapping profiles

These profiles guide per-resource registrations; they are not automatic assignments based only on the resource name.

| Evidence / example resources | Proposed representation | Things not to infer |
| --- | --- | --- |
| Physical/binary PPI: BioGRID subsets, DIP, HINT, HuRI, InnateDB, MatrixDB | `interacts_with`; `directly_physically_interacts_with` only for explicitly direct physical evidence | A network record or co-complex experiment alone does not prove direct contact. |
| Causal signaling: ACSN, SPIKE, SignaLink, Wang, TRIP | `affects` with exact direction/aspect/mechanism crosswalks where supported | Positive sign does not establish a particular aspect; directedness alone does not establish causation. |
| TF regulation: DoRothEA, CollecTRI, TRRUST, HTRIdb, PAZAR | Reviewed regulation relation, commonly `affects` with expression aspect and reported direction | Do not interpret binding-only evidence as changed expression; keep complex TFs and target RNA/gene domains. |
| miRNA/RNA: miRTarBase, miRecords, TransmiR, ncRDeathDB, lncrnadb | RNA-aware participants and source-supported regulatory/association predicates | Do not label all RNA-related targets as proteins or impose a universal negative sign. |
| Ligand–receptor: CellTalkDB, LRdb, CellCall, talklr, scConnect, Ramilowski2015 | `interacts_with` with explicit molecular/complex participants | A ligand/receptor label does not by itself demonstrate agonism or a downstream effect. |
| Drug/chemical relations: DrugBank, DGIdb, DDInter, CancerDrugsDB, SLC | Source-specific interaction, association or qualified effect mapping; retain measurements | Drug–drug interaction, transporter substrate and protein binding are different assertions. |
| Functional associations: STRING links, CPDB subsets | Source-reviewed `associated_with` or more specific predicate | STRING functional links must not all become physical PPIs; keep physical subset separate. |
| Negative evidence: Negatome | Explicitly negated assertion only after storage/consumer support is tested | Never encode absence of interaction as an inhibitory effect or a positive `interacts_with` edge. |
| Aggregators: PathwayCommons, reaction resources, ProtMapper | Preserve source/via labels and source-supported relation meaning | Never treat co-membership or process-derived edges as automatically causal; avoid double-counting native migrated resources. |

For a symmetric predicate, do not force a directed interpretation because the source lists A before B. If a source's directedness cannot be represented by the selected predicate, choose a faithful reviewed predicate or retain the row for review. Effect direction is not endpoint reversal. Negative/mixed observations must not disappear when a pair is canonicalized.

## 5. Verification and implementation sequence

Implement only after the resource's configuration is explicit. Suggested initial fixtures: ComplexPortal/PDB plus scConnect (dictionary versus set, stoichiometry and species); KEA/DEPOD plus iPTMnet (opposite PTM effects and namedtuple aliases); BioGRID plus TRRUST/CollecTRI (undirected versus causal, and nested complex participants).

1. Add the common runner and serializable normalizers; keep all resource discovery/import paths free of downloads.
2. Add complex mapping with source ID crosswalks, supported stoichiometry, accessionless identity and species validation. Keep HPMR and the 3DComplex stub disabled.
3. Add enzyme–substrate mapping only after evidence persistence demonstrably retains distinct sites and isoforms. Document the protein-level projection and remaining serving-site limitation.
4. Add the validated NetworkInput subset and explicit biological profiles. Expand support per unsupported configuration form, not by silently defaulting.
5. Enable inventory resources incrementally after representative fixtures and available full snapshots pass the actual build/extractor path. Avoid duplicating already-native v2 resources.

Required checks should include:

- Synthetic contract tests without network calls: dict versus set complexes, multiple source IDs, unknown namespace, default/known stoichiometry, absent accession, missing member and nonhuman request.
- PTM tests: phosphorylation versus dephosphorylation, unknown modification, alternative/candidate enzymes, missing enzyme/site, isoform coordinates, nucleotide versus protein RefSeq IDs, and publication count versus PMID.
- Network tests: exact sign/direction parsing, unknown/conflicting signs, reference delimiters/policies, positive/negative filters, endpoint taxonomy, loops, complex preservation, unsupported callbacks and stale entry points.
- End-to-end `SilverExtractor` and evidence round trips: nested membership/counts; native types/namespaces/qualifiers; separate site/evidence occurrences; no accidental inferred activity or positive negated edge. Verify any downstream deduplication preserves these distinctions.
- For each available snapshot, report returned legacy rows, normalized rows, emitted observations, excluded rows by reason, unresolved identifiers and unmapped fields. Explain intentional differences from the legacy core pipeline; do not call them parity.
- Refresh and cache identity tests, including function-specific intermediate caches and differing species/version/filter arguments.

Full live-source availability is separate from adapter correctness. A broken endpoint needs maintenance; a stub needs implementation; an unsupported semantic field needs a mapping decision. None is repaired merely by adding a generic wrapper.
