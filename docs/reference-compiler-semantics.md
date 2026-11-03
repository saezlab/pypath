# Reference-library semantics (branch parquet-migration) — technical brief

Paths are relative to `packages/`. `R` = `resolver/omnipath_resolver`, `B` = `build/omnipath_build/reference`,
`RS` = `resolver/rust/reference/src/lib.rs`, `PC` = `resolver/src/precomputed.rs`.
Scale (docs/reference-resolver.md): 288,430,295 entities, 1,891,149,589 admitted lookup keys, 173 GiB LMDB.
Nothing here is built locally (no data/reference), so all claims are from code, not from data inspection.

Overall pipeline: hub parquet (`source_type, source_id, hub_id, taxonomy_id, backend`) -> `build_reference.py`
(assigned reference: records/edges/graph/components/decisions/members/entities/claims/gene-products/gene-resolution/
source-gene-resolution) -> `finalize_source_reference.py` (RaMP fix-ups) -> `direct_compact_index.build_direct_compact`
(entities -> enrich -> products -> identifiers -> publish -> candidate_limit -> gene_role_index) -> LMDB runtime
(`R/index.py FullRuntime`) + Rust kernels (`PC`, `RS`). Ordering of compile stages: `B/direct_compact_index.py:299-327`.

---------------------------------------------------------------------------------------------------------------
## 1. Lookup key

### 1.1 Byte layout  (`R/observations.py:48-54`; mirrored in Rust `RS:92-111`, SQL `B/full_index_identifiers.py:236`)

```
byte 0      0x01                      format version
byte 1      target   1 = chemical library, 2 = gene/protein library
byte 2      route    1 = identity route, 2 = gene-resolution route   (3 exists only in compile-time tables)
bytes 3-4   namespace code, big-endian u16 = CODES[ns] = 1-based index in NS tuple
byte 5      scope flag: 0 = unscoped, 1 = taxon-scoped
[bytes 6-9  taxon, big-endian u32, ONLY if flag==1; python: int(scope).to_bytes(4,"big")]
rest        identifier, UTF-8, exactly as normalised (no terminator)
```
`bytes([1, kind, route]) + CODES[ns].to_bytes(2,"big") + (b"\x01"+int(scope).to_bytes(4,"big") if scope else b"\x00") + value.encode()`.
Offset of identifier = 6 (unscoped) or 10 (scoped) (`R/index.py:85`). `scope` is a string; any non-empty string (including "0")
makes the flag 1 (so scope "0" would encode 4 zero bytes; entity taxa are never "0" so it can only miss).
Runtime validation: `len>=7, key[0]==1, key[1] in (1,2), key[2] in (1,2), key[5] in (0,1)`, code must be in manifest
`namespace_codes` (`R/index.py:74-88`). The manifest pins CODES; any reorder of NS breaks the index (`R/index.py:45`).

NS order (code = position+1) (`R/observations.py:9-41`):
1 inchikey, 2 pubchem, 3 chebi, 4 chembl, 5 hmdb, 6 lipidmaps, 7 swisslipids, 8 bigg, 9 metanetx, 10 uniprot,
11 uniprot-sec, 12 uniprot_entry, 13 entrez, 14 ensg, 15 enst, 16 ensp, 17 hgnc, 18 refseq, 19 refseq_protein,
20 genesymbol, 21 genesymbol-syn, 22 inchi, 23 kegg, 24 cas, 25 drugbank, 26 goslin, 27 refmet, 28 genbank, 29 ramp,
30 ramp_gene, 31 kegg_gene.

Storage envelope: shard = `md5(identifier)[:2]` (identifier only, NOT the full key; `R/index_storage.py:13`), physical LMDB key
`\0+key` (or `\1+sha256(key)` if len>500), value prefixed J/C (chunked 64 KiB)/L-V envelopes (`index_storage.py:17-65`).

### 1.2 target (byte 1)
Chosen per observation from the entity-type policy: `target = 1 if library=="chemical" else 2` (`observations.py:108`).
Policy by biolink type (`R/canonical/policy.py:58-203`): gene_protein = Gene + all `protein` descendants + all `transcript`
descendants; chemical = ChemicalEntity + SmallMolecule; cv_term/complex/generic: `library=None`, never matched
(`match.py:179` skips; they keep native identifiers).

### 1.3 route (byte 2) — decided PER VOTE (`observations.py:116`)
`route = 2 if target==2 and ns in {"entrez","ramp_gene","kegg_gene"} else 1`. Route 2 = "this identifier is a *gene* source
id; resolve to gene entity (independent of products)". Everything else (chemical, uniprot, ensg, hgnc, symbols, refseq, ...) is route 1.
Compile-time route 3 (`_product_entrez`, `_product_ensg`) is a vestigial product adjacency never served (see 2.7).
Runtime accepts only route 1/2 (`index.py:77`).

### 1.4 namespace code
Namespace after `normalize_ns` (glossary slug via `omnipath_core.naming.normalize_namespace`; remaps `bigg_metabolite->bigg`,
`uniprot_trembl->uniprot`; `R/canonical/identifiers.py:25-29`). A vote with a namespace not in CODES raises (`observations.py:113`).
Only namespaces in `policy.match_namespaces` become id votes; `policy.symbol_namespaces` (`genesymbol`, `genesymbol-syn`, gene_protein only)
become symbol votes ONLY if the observation has a taxon (`match.py:107-110`):
- gene_protein match ns: uniprot, uniprot-sec, uniprot_entry, entrez, ensg, enst, ensp, hgnc, refseq, refseq_protein, genbank, ramp_gene, kegg_gene (`policy.py:62-78`).
  `ramp` observed namespace is renamed `ramp_gene` when library is gene_protein (`match.py:83-84`).
- chemical match ns: inchikey, inchi, chebi, chembl, pubchem, hmdb, kegg, cas, lipidmaps, swisslipids, drugbank, bigg, metanetx, goslin, refmet, ramp (`policy.py:110-128`).
- Names, synonyms, smiles, formula etc. never vote (they populate `observed` only).

### 1.5 taxon-scope flag (byte 5) — decided PER OBSERVATION, not per vote (`observations.py:110,115`)
```
symbol_scoped = (target == 2) and any(v.kind == "symbol" for v in votes)
scope         = (v.taxon or "") if symbol_scoped else ""      # applies to EVERY vote of the observation
```
So: a gene/protein observation with any genesymbol/genesymbol-syn vote (needs taxon) looks up ALL its votes (uniprot, entrez,
ensg ...) under the observation's taxon scope; observations without a symbol vote look up everything unscoped and the observation
taxon is ignored. Chemicals are never scoped. The scoped key matches only candidates whose entity `taxon` equals the scope
(see 2.4). Symbols exist ONLY scoped in the index (no unscoped symbol keys are compiled) (`full_index_identifiers.py:226,274-304`).
Scope must be digits (`int(scope)`); an obs taxon like "Homo sapiens" would raise.

### 1.6 Which namespaces are "taxon-scoped"
Index-wise every non-symbol key exists both unscoped and per taxon (present entity taxa); symbols only per taxon. Runtime-wise
the scope is requested only under the rule 1.5 (so effectively: genesymbol, genesymbol-syn always scoped; everything else scoped
only when co-observed with a symbol).

### 1.7 How each vote / key is built (`match.py:73-138`, `identifiers.py:55-70`, `observations.py:104-146`)
Inputs per observation `obs`: `(namespace, identifier)` primary (primary=True), `obs.identifiers[*]` ({ns|type, id|identifier}),
for gene_protein also `molecular_form.sequence_identifiers[*]` and `molecular_form.isoform_identifier` (all primary=False).
Normalisation `normalize_identifier(ns,value)` -> list of (ns,id):
- `ensembl` ns + regex `ENS[A-Z]*([GTP])\d+(\.\d+)?` -> ensg/enst/ensp. ensg/enst/ensp: trailing `.\d+` stripped (`normalize_id`).
- chebi -> `CHEBI:<digits>`; hgnc -> `HGNC:<digits>`; hmdb -> `HMDB` + 7-digit zero-padded (strip leading zeros then zfill(7));
  inchikey -> upper, strip `InChIKey=`; chembl -> upper; everything else unchanged except trim. (SQL twin `normalize_id_sql`, `identifiers.py:73-96`;
  the build uses the SQL form on hub ids/values, so index keys and votes agree. NB SQL `hmdb` branch pads to >=7 via lpad.)
- **uniprot isoform**: if slug==uniprot and id matches `^([A-Z][0-9][A-Z0-9]{3}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)-\d+$` two pairs are emitted:
  `(uniprot,"P04637-2")` AND `(uniprot,"P04637")` -> two votes (intersected later). `-PRO_n` chain suffix does not match, so no base vote.
- **uniprot-sec**: votes as ordinary ns `uniprot-sec` (CODES 11). The index also stores secondary accessions under ns `uniprot` (see 2.3b),
  so an obs with `uniprot:Q9XXXX(secondary)` resolves; no special runtime rewrite.
- **refseq typing** (`match.py:93-101`): ns refseq/refseq_protein re-typed by prefix: `AP|NP|XP|YP|WP|ZP_` -> `refseq_protein`; `NM|NR|XM|XR_` -> `refseq`;
  `.version` stripped for those (version retained in molecular form only). Unrecognised prefix keeps its ns and keeps version.
- **genesymbol-syn**: ordinary symbol vote when observed as ns genesymbol-syn. In `FullRuntime.resolve` (`index.py:130-137`): if the same input has any
  `genesymbol` vote whose lookup returned candidates, that input's `genesymbol-syn` votes are dropped.
- **SMILES -> InChIKey** (`match.py:128-137`; `canonical/structures.py:30-70`): chemical only, only when `observed["inchikey"]` is empty and smiles observed.
  RDKit parse (no CXSMILES, no name parse); reject 0 atoms, dummy atoms (generic_structure), non-`InChI=1S/`, key whose chars 23:25 != "SA".
  Derived key is added as `Vote("inchikey", key, None, "id", primary=False)` and, like a full InChIKey, becomes a PrimaryAnchor (anchor field set,
  `observations.py:117`). Result cached in sqlite keyed on smiles and (policy version, rdkit, inchi versions).
- **goslin** (`observations.py:84-101`): chemical only, when `all votes are refmet` or there are none (`any(v.ns!="refmet")` returns early otherwise).
  Parse every observed `name`+`synonym` with pygoslin (`R/goslin.py:9-35`): status parsed requires level >= SPECIES, no adduct; value
  `"<level>:<lipid string>"` with sorted chains. Keep only parsed names of MAX specificity; each distinct goslin string -> `Vote("goslin", value, None, "id", False)`.
- **Bundle** `observation_bundle` (`observations.py:104-146`): set of `(ns,id,scope,route,anchor,primary)` sorted; `ordinal` = position in that sorted list;
  `anchor = v.id` iff `target==1 and ns=="inchikey" and INCHIKEY_RE (^[A-Z]{14}-[A-Z]{10}-[A-Z]$) fullmatch`, else "".
  Row = `input_id, ns, identifier, scope, anchor, target, route, ordinal, primary, lookup_key`. If `target==2` and the entity type is an RNA/transcript
  type (`RNA_ENTITY_TYPES`, descendants of biolink `transcript`), every row gets `gene_only=True` (forces gene-level interpretation of all votes;
  RNA never asserts a protein). `primary` is never read by any kernel.
- Observations whose library is not loaded are skipped (`match.py:179`); a vote ns missing from CODES raises (cannot happen: match namespaces are all in CODES); empty evidence -> NotFound/abstain.

---------------------------------------------------------------------------------------------------------------
## 2. Identifier-index value (what `FullRuntime.lookup(key)` returns)

```
{ "candidates": [ [num, entity_id, kind, anchor, quarantined, reviewed, gene_ids?], ... ],   # ordered by num
  "gene": bool, "products": bool }
```
(`R/index_codec.py:15-109`: candidate packed as `[num, eid, kind, anchor-or-None-if==eid, flags(bit0 q, bit1 reviewed, bit2 anchor==eid), (genes)]`,
top-level flags bit0 gene, bit1 products.) Miss => `candidates=[]`, `gene = target==2 and ns in {hgnc,ensg,genesymbol,genesymbol-syn}`,
`products = target==2 and ns in {hgnc,ensg}` (`index.py:96-100`) — irrelevant to decisions because an empty posting is skipped by both kernels.

### 2.1 Candidate fact fields
Built in SQL as `json_array(id, entity_id, kind, anchor, quarantined, reviewed, gene_ids)` (`full_index_identifiers.py:20`):
- `num`: numeric entity id (see section 5).
- `entity_id`: `inchikey:<27>`, `uniprot:<ACC>`, `entrez:<N>`, or a native record id `<hub>:<local_id>` (e.g. `pubchem:123`, `ramp_gene:RAMP_G_x`, `goslin:...`).
- `kind`: 1 chemical, 2 protein, 3 gene. Set from entity_id prefix: `inchikey:`->chemical, `uniprot:`->protein, else domain default
  (chemical domain->chemical; gene_protein domain->gene) (`build_reference.py:608`). So native `ramp_gene:*` / anchorless records are kind 3 but are not `entrez:*`.
- `anchor`: entity `primary_anchor` = `min(anchor) where anchor_count=1` over members (`build_reference.py:610`), i.e. `inchikey:..`/`uniprot:..`/NULL. Anchored entities: anchor==entity_id.
- `quarantined`: `bool_or(decision='multiple_anchor_claims')` over member records: a record claiming >1 distinct full anchor (chemical) is its own entity, flagged.
- `reviewed`: proteins only: `records-uniprot.reviewed` joined on anchor; true iff some uniprot_entry claim of the record has `split_part(entry,'_',1) <> local_id`
  (Swiss-Prot style names; TrEMBL entry names start with the accession) (`build_reference.py:452,410`; `full_index.py:165-175`). Chemicals/genes false. **Not used by any current decision** (Policy disables reviewed preference).
- `gene_ids`: list of NCBI gene ids as strings. `entrez:N` entity -> `[N]`; protein -> `list(DISTINCT entrez_id)` from `gene_products` grouped by `protein_entity_id`; others `[]`
  (`full_index.py:182,192`). (The entity *record's* `gene_ids` differs: `[]` for entrez entities, see 3.)

### 2.2 Admission of an assertion into the identifier index (`full_index_identifiers.py:208-226`)
Raw assertion rows (`full_index_assertions.py`) are `(target, route, namespace, identifier, entity_id, tag)` then LEFT JOINed to entity metadata. Admitted iff
1. tag != 'chemical_fallback' OR no row with same (target,route,namespace,identifier) tagged 'native' exists (native beats fallback);
2. the entity exists in the entity table (`id IS NOT NULL`);
3. `kind == target` OR (`target==2 AND kind==3`) (gene entities allowed in the gene/protein library; chemical lib: kind 1 only).
Quarantined entities ARE admitted (decided by kernel). Candidate rows are DISTINCT over (target,route,ns,identifier,id,...,taxon,gene_ids), list ordered by `id`.
Rows tagged `product` are excluded from the identifier index.

### 2.3 Raw assertion rows by (target, route) (`full_index_assertions.py:57-159`)
Source relations are projections of the assigned reference (see section 7). Notation: `claims-X` = `claims-X/identifier_claims.parquet (entity_id,namespace,identifier)`
(all roles except `attribute` source types; role itself unused by compiler); `identity` = `{domain}-entities/identity_identifiers.parquet (namespace,identifier,entity_id)`.

**Target 1, route 1 (chemical)**
- from `identity` where namespace in (inchikey, inchi, chebi, chembl, pubchem, hmdb, kegg, cas, lipidmaps, swisslipids, drugbank, bigg, metanetx, goslin, refmet, ramp), tag `native`.
  In practice identity namespaces are: the hub names (chebi,pubchem,chembl,hmdb,lipidmaps,swisslipids,bigg,metanetx,refmet,ramp) with the hub-normalised local id; `inchikey` rows
  `substr(entity_id,10)` for every entity whose id starts `inchikey:`; `goslin` rows from `goslin-identifiers/claims` joined to members (`build_reference.py:625,640`).
  `inchi`, `kegg`, `cas`, `drugbank` never appear in identity => unreachable as native keys.
- from `claims-{hub}` for hub in chebi,pubchem,chembl,hmdb,lipidmaps,swisslipids,bigg,metanetx,refmet,ramp,(kegg if present): namespaces `cas, drugbank, kegg`
  (+ `chebi` only from the chebi hub, + `bigg` only from the bigg hub) (`projections.py:47-57`). tag `chemical_fallback` for kegg/bigg, `regular` for cas/drugbank/chebi.
  Other claim namespaces (cross-refs like hmdb ids inside chebi claims) are NOT lookup keys (only entity aliases).
- Kernel consequence: a posting for `inchikey` is exactly one entity (`inchikey:<key>`); `goslin` postings may be many (limit-10 rule).

**Target 2, route 1 (gene/protein)** (protein_ns whitelist for identity: uniprot, uniprot-sec, uniprot_entry, ensg, enst, ensp, hgnc, refseq, refseq_protein, genbank, genesymbol, genesymbol-syn; `full_index_assertions.py:58`)
 a. `identity` rows in that whitelist: i.e. `('uniprot', local_id) -> entity of that uniprot record` (tag regular). (`entrez` and `ramp_gene` identity rows are excluded here.)
 b. `claims-uniprot` (= `gene_protein-uniprot-forward`) rows, ns in exact set {uniprot_entry, ensg, ensp, enst, hgnc, refseq, refseq_protein, genbank, genesymbol, genesymbol-syn} as-is, plus expansions:
   - `genesymbol` rows are ALSO emitted as `genesymbol-syn` (regular): an exact symbol answers a synonym query.
   - `genesymbol-syn` rows are emitted as BOTH `genesymbol` and `genesymbol-syn` with tag `symbol_synonym` (weaker; see 2.5).
   - `uniprot-sec` rows are emitted as BOTH `uniprot` and `uniprot-sec` (so secondary ACs resolve under `uniprot`). (Neither `uniprot` nor `uniprot-sec` from claims is emitted as-is; `uniprot` comes via (a).)
   - `refseq_protein` matching `(AP|NP|XP|YP|WP|ZP)_[0-9]+` and `genbank` matching `[A-Z]{3}[0-9]{5,}` also emitted with `.version` stripped (versioned original kept by the as-is rule).
   - `refseq` matching `(NM|NR|XM|XR)_[0-9]+` (after version strip) emitted version-stripped (also from `claims-entrez`).
   Entity of these rows is the *protein* entity of the uniprot-hub record (`uniprot:ACC`). NOTE: the claims `entrez`/uniprot `ensg` etc. attach to proteins; they reach genes only by collapse (2.6).
 c. `claims-entrez` (= `gene_protein-entrez-forward`) rows with ns in (ensg, enst, refseq) -> entity `entrez:N` (gene). (`ensp`/`refseq_protein` of gene2ensembl deliberately NOT asserted: they are product ids.)
 d. **Isoform projection** (`full_index_assertions.py:26-38,164-167`): any asserted entity_id matching `uniprot:[A-Z0-9]+-[0-9]+` whose parent (`-N` stripped) has a `records-uniprot` row with `anchor_count=1`
    is replaced by the parent entity_id (non-product tags). Isoform entities without a primary parent stay as-is (their anchor contains '-', hence never a "product" in the kernel).

**Target 2, route 2 (gene resolution; candidates are GENES, never proteins)**
 - `gene-resolution` (`entrez_identifiers.parquet`): `('entrez', N) -> entity entrez:N` for (i) every identity row with namespace entrez (an entrez hub record) [admission `gene_identity`],
   (ii) every `gene_products.entrez_id` even without an entrez hub record [`explicit_gene_product`] (`build_reference.py:90-99`). Posted as key(2,2,'entrez',scope,N).
 - `source-gene-resolution` (`identifiers.parquet`): `('ramp_gene', RAMP_G_x) -> entrez:...` (admission `explicit_source_gene`) or `-> native ramp_gene entity` (`source_gene_only`):
   ramp_gene hub rows with source_type in (uniprot, entrez, ensg, ensp, enst, hgnc, refseq_protein, genbank) are matched to uniprot claims with equal value and ns
   (or ns=uniprot vs claim ns uniprot-sec) -> protein entity -> gene via `gene_products` (taxon compatible: ramp taxon in ('','0') or equal to product taxon); entrez source ids map directly to
   `gene-resolution` entities; `source_gene_only` adds ramp_gene identity rows whose entity is not `uniprot:*` and not otherwise mapped (`build_reference.py:728-746`).
 - There is no route-2 builder for `kegg_gene` (code reserved; always misses).
 - "Gene->protein" routing (Rust `Role::GeneToProtein`, `allow_gene_to_protein`) is DISABLED in current policy (`RS:174-208`; `PC:99` uses `Policy::omnipath()` with all gene->protein flags false). In the molecular kernel route==2 only means "this vote is a gene vote".

**Route 3 / product rows** (`full_index_assertions.py:142-159`): `gene-products-forward` (`_product_entrez`, protein entity) and `claims-uniprot` ensg (`_product_ensg`) tagged `product`; compiled to a temp LMDB
(`product_partition`, `full_index_identifiers.py:78-101`) and only read by `Products.get/projected`, which `identifier_partition` never calls (only `Products.gene` is used, which reads the entity store).
=> vestigial; the `products` flag is always `false` in compiled values (`full_index_identifiers.py:243`; symbols/ensg-fallback/collapse also set False).

### 2.4 Scoped copies (`full_index_identifiers.py:224-230`)
For non-symbol namespaces: key written once unscoped (all candidates) and once per `taxon` that occurs among its candidates, with only candidates of that taxon.
`taxon` here = candidate entity's taxon (min non-'0' taxon of its member records; NULL if none; finalize replaces protein taxon by authoritative UniProt taxon for RaMP-touched proteins,
`finalize_source_reference.py:100-111`). Candidates with NULL taxon appear only in the unscoped list. A requested scope with no candidate of that taxon = key absent = miss.
Symbols (`genesymbol`, `genesymbol-syn`): only taxon-scoped keys, only for candidates with non-empty taxon, grouped by (ns, identifier, taxon); `symbol_synonym`-tagged rows are dropped
if any non-synonym row exists for same (ns, identifier, taxon) among admitted rows (exact beats synonym) (`:274-279`).

### 2.5 `gene` flag (compile) (`full_index_identifiers.py:242`)
`gene = target==2 AND (route==2 OR namespace IN ('hgnc','ensg','enst','refseq'))`; symbols always `gene=True`; ensg fallback `gene=True`. Gene flag FALSE for uniprot, uniprot-sec, uniprot_entry,
ensp, refseq_protein, genbank, and all chemical keys. Runtime: `gene_flag OR vote.gene_only` -> `is_gene` in kernel (`index.py:145`; `PC:184`).

### 2.6 Gene collapse — how "gene level" candidates are formed (`collapse_gene_posting`, `full_index_identifiers.py:344-368`)
Applies to keys with target 2 and (route==2 or ns in hgnc,ensg,enst,refseq code set) and to all symbol keys (`:263-268, :291-300`). For a posting with `gene=True`:
- if ANY candidate is quarantined -> posting left unchanged;
- each candidate must be either a gene entity (`kind==3` and id startswith `entrez:`) or have non-empty `gene_ids` (candidate[6]); otherwise the whole posting is left unchanged;
- genes = union of those gene ids; every gene must exist in the entity store as kind 3 and non-quarantined, else posting left unchanged;
- result replaces candidates by the gene facts `[num, 'entrez:N', 3, anchor(null), False, reviewed, ['N']]` sorted by num, `gene=True, products=False`.
So hgnc/ensg/symbol/route-2 postings are normally gene entities; the candidate-limit counts genes, not proteins (a gene with 40 protein entries stays 1 candidate). Product-ish keys
(uniprot, uniprot_entry, ensp, refseq_protein, genbank) are never collapsed: they stay protein facts with their `gene_ids`.
Dead/redundant: an "ensg fallback" for noncoding genes from `claims-entrez` (`:305-325`) only fires for ensg ids lacking any target2/route1 ensg row, which (2.3c) already provides => expected to be ~0.

### 2.7 Gene-role component override (`R/gene_role_index.py:94-102`, `B/gene_role_index.py`)
At runtime, before the base index, for `key[1:3]==(2,1)` (target 2, route 1) and ns in {ensg, hgnc, genesymbol, genesymbol-syn, enst, refseq}: if the component has the key, ITS value
(`{gene:true, products:false, candidates:[gene facts]}`) replaces the base value — even an empty posting (retained on purpose when >10 distinct genes: `B/gene_role_index.py:232-241`).
Absent => base lookup. Component assertions (`B/gene_role_index.py:37-85`):
 1. `claims-entrez` rows with ns in (ensg, enst, hgnc, genesymbol, genesymbol-syn) or refseq `(NM|NR|XM|XR)_\d+(\.\d+)?` attached to NON-quarantined kind='gene' entities (entities.parquet), taxon = gene taxon.
 2. `claims-uniprot` rows with ns in (ensg, hgnc, genesymbol, genesymbol-syn) of admitted (kind=protein, non-quarantined) proteins that link to EXACTLY one entrez gene
    (`gene_products` grouped by protein having `count(distinct entrez_id)=1` and at most one distinct non-0 taxon), the gene must exist, and gene/link/protein taxa must agree (NULL/''/'0' = wildcard);
    gene taxon coalesced as gene.taxon, link.taxon, protein.taxon.
 ensg/enst/refseq identifiers have `.version` stripped. Taxon `'0'`/'' -> NULL.
 Symbols expand to both ns genesymbol and genesymbol-syn; within (ns, identifier, taxon) exact assertions (is_synonym=false) hide synonyms (`:211-213`).
 Keys: unscoped for non-symbols; taxon-scoped for any assertion with a taxon (incl. non-symbols); no unscoped symbol keys.
 Candidates = distinct gene entity ids; if > 10 -> `candidates=[]` (archived as "rejected_gene_keys"), else `[num,'entrez:N',3,anchor,False,reviewed,['N']]` sorted by num.
 Component records restore aliases only if the gene survives in the alias' own key (3.5).

### 2.8 Why `products` is meaningless today
`products` is passed into the Rust 6-tuple (`index.py:146`) but `resolve_molecular_batch` ignores it (`PC:170 ... _`), and the chemical path only uses it via `Role::GeneIdentity{multiple_products}`
which is unreachable for chemicals (gene flag false). An on-demand implementation can drop it.

---------------------------------------------------------------------------------------------------------------
## 3. Entity-index value (`FullRuntime.record(entity_id)` -> record dict)

`{entity_id, kind(1/2/3), anchor, taxon, label, identifiers, gene_ids}` (the stored object is `{record, meta}`; meta = `{id, entity_id, kind, anchor, quarantined, reviewed}`; runtime returns only `record`
after the gene-role merge, `index.py:102-108`). `match.py:213-230` consumes `identifiers` (as `[ns,id]` pairs), `taxon or None`, `label` (None if label == identifier part), `kind`
(1 -> reference_library chemical, else gene_protein), `gene_ids` (protein -> `entrez:`+id as `protein_gene_candidates`). `writer.py:608-660` reads `identifiers` for every referenced node_id.
Missing entity -> ValueError.

### 3.1 Which entities exist
`{chemical,gene_protein}-entities/entities.parquet (entity_id, kind, primary_anchor, taxon, quarantined)`:
group of all member records' entity ids (`build_reference.py:594-640`) + extra rows: gene_protein: `entrez:N` for every `gene_products.entrez_id` not otherwise present (kind gene, taxon=min(p.taxon), not quarantined);
chemical: every anchor in `multiple_anchor_claims` not already an entity (kind chemical, anchor = itself). Entity assignment summary (`build_reference.py:479-585`):
anchored record (anchor_count=1) -> its anchor; entrez hub records -> `entrez:N` always; anchorless components -> attached to the single boundary anchor (if exactly one distinct and none of the boundary
claims is multi-anchor), else root record id if no anchor, else own record id ("ambiguous_native"); multi-anchor records -> own record id, quarantined. Gene<->product edges (`semantics='gene_product'`) never join entities.
Anchors: chemical = full InChIKey `^[A-Z]{14}-[A-Z]{10}-[A-Z]$` excluding MOSFIJXAXDLOML-UHFFFAOYSA-N / -UHFFFAOYNA-N; protein = uniprot record where `source_type='uniprot' AND normalized_id=local_id`.
`taxon` = `min(taxon) FILTER (taxon<>'0')` of member records (record taxon = min non-''/'0' hub taxonomy_id else '0').

### 3.2 label (SQL, `full_index.py:194-225`; Python twin for corrected records `B/identifier_admission.py:6-47`; display twin `R/canonical/label.py`)
From the entity's alias set (3.3): eligible namespaces and priority
- chemical: `name`(0), `chebi`(1), `pubchem`(2); value not matching `[A-Z]{14}-[A-Z]{10}-[A-Z0-9]` and `1<=len(trim)<=120`;
- gene_protein: `genesymbol`(0), `name`(1); `1<=len(trim)<=120`. (`genesymbol-syn` is never a label.)
Select `arg_min` over `(priority, long_flag, len, value)`, `long_flag` = genesymbol len>15, name len>80 (chemical) / >120 (gene_protein), chebi/pubchem have len=0 / never long.
Formatting: chebi -> prefix `CHEBI:` if missing; pubchem -> `CID:<id>`; else value. Fallback `substr(entity_id after first ':')`.
Gene-role merge recomputes label for `entrez:` entities: min over merged `genesymbol` values (1..120 chars) by `(len>15, len, value)` (`R/gene_role_index.py:116-120`).
`match.py:224`: `label = row["label"] if != identifier else None`, then `assign_label` fallback by obs policy (label.py: chemical name sorted by (len>80,len,s) etc.).

### 3.3 identifiers (alias list)
`list([namespace, identifier] ORDER BY namespace, identifier)` of DISTINCT union (`full_index.py:103-115`, `direct_compact_index.py:114-126`) of every projection named `{domain}-*-forward` or `{domain}-*-names` plus `gene-resolution-forward` (gene_protein):
- chemical: `chemical-identity-forward` (identity_identifiers: hub-native ids incl. pubchem, inchikey, goslin), `chemical-{hub}-forward` for hub in chebi,chembl,hmdb,lipidmaps,swisslipids,bigg,metanetx,refmet,ramp (ALL claim namespaces/roles of
  that hub's records except attribute types and smiles/inchi/formula/names; `pubchem`/`kegg` claim tables are excluded from aliases), `chemical-{hub}-names` (name, synonym, lipid_shorthand, systematic_name from raw hub rows for those hubs).
- gene_protein: `gene_protein-identity-forward`, `gene_protein-{uniprot,entrez,ramp_gene}-forward` (all claims), `gene_protein-ramp_gene-names` (names), `gene-resolution-forward` (`entrez` id of entrez:N entities), and
  `gene_protein-gene-forward` (gene entity `entrez:N` receives `genesymbol,genesymbol-syn,hgnc,ensg` claims of each linked protein ONLY IF that protein links to exactly one gene: `projections.py:90-104`).
=> The alias list is a SUPERSET of the lookup keys (includes names, cross-refs, `entrez` ids on proteins, `ensp`, `uniprot_entry`, isoform ids, ...). The writer later filters by `policy.alias_namespaces`
   (`policy.py:83-100,134-153`; gene_protein: uniprot, uniprot-sec, uniprot_entry, entrez, ensg, enst, ensp, hgnc, refseq, refseq_protein, genbank, ramp_gene, kegg_gene, genesymbol, genesymbol-syn;
   chemical: inchikey, chebi, chembl, pubchem, hmdb, kegg, cas, lipidmaps, swisslipids, drugbank, bigg, metanetx, goslin, refmet, ramp, name) and drops form-specific ids unless they are the entity's own
   canonical id (`writer.py:140-148, 626-655`: ensp, enst, `*_sequence_sha256`, uniprot isoforms/chains, refseq NP/XP/NM/... accessions).
Admission rules applying to the list itself: none beyond the candidate-limit rewrite (4.3). Identifiers are NOT filtered by "was it indexed".

### 3.4 taxon, anchor, gene_ids
`taxon`: 3.1 (string or null; `""` treated as none by Python). `anchor`: `primary_anchor`. `gene_ids` (record): protein -> sorted distinct `entrez_id` from `gene_products`; others `[]` (note: for `entrez:` entities the record's gene_ids is `[]`,
the candidate fact's is `[N]`).

### 3.5 gene-role component merge (`R/gene_role_index.py:104-121`)
Only for ids starting `entrez:`: load `{"identifiers":[[ns,id],...]}` from the component `records` store; merged identifiers = sorted set-union with the base record; label recomputed from merged `genesymbol`s. Aliases in the component are the
(ensg, hgnc, genesymbol, genesymbol-syn) assertions of that gene, restored only if the gene appears in the supported posting of that very alias key (scope=taxon for symbols else ''), `B/gene_role_index.py:273-283`.

### 3.6 gene_products (`build_reference.py:69-87`) — the only gene<->protein link
Union of (a) uniprot-hub edges with `semantics='gene_product'` and `source_anchor_count=1`: `(entrez_id=target_id, protein_entity_id=source anchor, taxon=s.taxon, reviewed=s.reviewed)` unless both taxa known and differ
(target record may be absent); (b) entrez-hub edges to uniprot with `target_anchor_count=1`: `(entrez_id=s.local_id, protein_entity_id=t.anchor, taxon=t.taxon, reviewed=t.reviewed)` requiring compatible taxa. `projection_eligible` column is unused downstream.

---------------------------------------------------------------------------------------------------------------
## 4. Candidate limit and other filters

### 4.1 The >10 rule
- Policy `max-10-candidates-per-scoped-lookup-v1` (`candidate_limit.py:22`). Evaluated on the FINAL compiled value of each physical key: `len(value["candidates"]) > 10` (audit in `CompactWriter.put_objects`,
  `compact_index.py:57-60`, `audit_limit=10`; archive asserts `len(ids) > 10`, `candidate_limit.py:97`). So it is applied AFTER gene collapse (counts genes for gene keys; proteins for product keys; entities for chemicals) and per (target, route, ns, identifier, SCOPE):
  unscoped and each taxon list judged independently (a symbol can be dropped globally and kept in taxon 9606).
- Action: `apply_limit` (`:295-409`) deletes the whole key from the identifiers LMDB (never truncates), archives key + all candidate ids + gene/products flags to `ambiguous.parquet` (schema `candidate_limit.py:23-37`; audit-only, not read at runtime),
  then rewrites entity records to remove aliases that correspond to rejected keys (`Rules.excluded`, `:264-288`): for each alias (ns,id) of an entity whose kind-target (`1 if kind==1 else 2`) has a rejected key with any route: excluded if the entity's taxon is among rejected scopes;
  or if the unscoped (`None`) key is rejected AND the entity has no usable taxon-scoped key (looked up in identifiers shards via `lookup_key(target,route,ns,str(taxon),ident)`). `clean_record` drops them and, if anything was dropped, recomputes the label (3.2 twin, fallback entity_id local part).
  Semantics for matching: a removed key is simply absent => behaves as "no evidence" (NOT a conflict).
- Gene-role component uses the same cap on distinct supported genes, but retains an empty posting (blocks fallback to base) and drops the restored aliases.
- Archived scale: 218,173 of ~1.89 B keys (0.0115%), 6.7 M candidate links.

### 4.2 Everything else that removes keys or candidates
- compile: chemical_fallback (kegg, bigg claims) suppressed where a native key exists (2.2.1); candidates whose entity not in entity table; kind != target (chemical lib only chemicals; gene/protein lib only protein+gene);
  `symbol_synonym` suppressed by exact within taxon; symbols require entity taxon (so taxon-less proteins have no symbol keys); scoped copies only for entities with taxon; isoform assertions re-pointed to primary parent;
  versioned refseq/genbank forms are asserted both with and without version only for syntactically valid accessions; version stripped for ensg/enst/ensp (hub emit + normalizer) and refseq RNA.
- gene collapse refuses (keeps protein facts) if any candidate quarantined / unlinked.
- goslin: only the most specific parsed level per hub record is admitted as a claim (`goslin_identifiers.py:178-184`); virtual `goslin:*` record joins components so conflicting anchors show.
- hubs only: UniProt hub idmapping types kept (`hubs/sources/uniprot.py:30-42`): UniProtKB-ID, Gene_Name, Gene_Synonym, GeneID, Ensembl, Ensembl_TRS, Ensembl_PRO, HGNC, RefSeq, RefSeq_NT, EMBL-CDS + uniprot-sec; RefSeq/genbank sequence mappings added by `prepare_sequence_mappings.py`.
- runtime: `genesymbol-syn` votes dropped when same input's `genesymbol` has candidates (`index.py:130-137`); molecular kernel errors if a posting > 100,000 ids or batch > 5,000,000 postings / 1,000,000 candidates; chemical kernel marks > 100,000 ids Incomplete, >10,000 final candidates Incomplete (all unreachable after the 10-cap).
- runtime flags: 64 MiB per stored value budget; missing partitions fail (no fallback).

---------------------------------------------------------------------------------------------------------------
## 5. Numeric entity ids (`num` / `id`)

Assignment (`full_index.py:171`): per entity partition `part = md5(entity_id)[:2]` (hex, 00..ff) and per domain (chemical=0, gene_protein=1):
```
id = part_int * 2**40 + domain_index * 2**39 + row_number() OVER (ORDER BY entity_id)      -- row_number starts at 1
```
Hence ids are a pure function of (entity_id, set of entity_ids in same md5-prefix partition and domain), stable only within one build; unique across the universe; used only as (a) the posting list ordering/identity for the Rust kernel
(strictly increasing, `PC:33`, `RS:330`), (b) join key in compiled facts. The gene-role component copies `meta["id"]` from the base so ids must match the same base generation (`B/gene_role_index.py:149-157`; manifest sha pins this).
Leaves the runtime? NO. `FullRuntime.resolve` maps ids back to `entity_id` strings in Rust (`PC:100-120`, molecular returns strings) and returns `entities=sorted(ids)` / strings; `match.targets` builds `Node` from `eid`;
`writer.py` only uses `node_id = entity_id` (e.g. `writer.py:281,336,412,611-626`); `candidate_count` (len of surviving candidates / distinct ids) is returned but unused. Grep of `packages/` shows no consumer of `num`.
Only ambiguous.parquet / manifest carry entity_id strings, not nums. An on-demand design can use entity_id strings (or any local dictionary) for intersection; sorting by entity_id string reproduces the externally visible ordering
(outputs are sorted by entity_id, `index.py:181`, `match.py:212`).

---------------------------------------------------------------------------------------------------------------
## 6. How the decision is made

Inputs per observation: (target; ordered votes with `lookup_key, anchor, route, ordinal, gene_flag, products_flag`) plus, per candidate id in any returned posting: `(entity_id, kind, anchor, quarantined, reviewed, gene_ids[])`.
`FullRuntime.resolve` (`index.py:116-199`) batches: unique keys -> `lookup()` -> `postings=(key,[num...])`, `metadata[num]=candidate tuple`; drops syn votes; splits queries by target; chemical -> `resolve_precomputed_batch`, biological -> `resolve_molecular_batch`;
then fetches `record()` for accepted ids and (molecular) the chosen protein.

### 6.1 Chemical (`PC:24-124` -> `RS:283-624`, Policy::omnipath: same-connectivity multiple OFF, require_primary_anchor OFF, gene->protein OFF)
1. Role per vote: `anchor != ""` -> PrimaryAnchor(InChIKey syntax-checked); else `Identity` (route 2/gene never occur for chemicals).
2. If ANY vote is a PrimaryAnchor, ONLY primary-anchor votes participate; all others are `IgnoredByPolicy` (`RS:209-217`). Primary anchors come from explicit/derived full InChIKeys only (1.7).
3. Stop conditions in order: bad input; >1 DISTINCT primary anchors -> ConflictingEvidence; any participating posting over limit -> Incomplete; a primary anchor with no posting -> UnresolvedPrimaryAnchor;
   no participating non-empty posting -> NotFound (empty/missing postings of NON-anchor votes are silently skipped = no veto).
4. Candidates = intersection of all non-empty participating postings (ascending size). Empty -> ConflictingEvidence. >10,000 -> Incomplete.
5. Metadata checks: candidate kind must be chemical else hard error; ANY candidate quarantined -> QuarantinedCandidate; claimed anchors must equal candidate anchors (else hard error: the inchikey posting must be exactly one entity, `RS:339`).
6. `len==1` -> Unique (accepted); `>1` -> Ambiguous (not accepted; MultipleSameConnectivity disabled). Reviewed preference only for gene roles (n/a).
Accepted entity = that one `inchikey:*` / native-record entity. Not-accepted -> match.py keeps native id, but if observed has exactly one valid full InChIKey it becomes the identity (`match.py:309-323`).

### 6.2 Gene/protein (`PC:141-293`, `resolve_molecular_batch`)
Per vote (skip if posting empty): `is_gene = gene_flag(or gene_only) OR route==2`. For each candidate m in posting: record id into `count`, `quarantined |= m.quarantined` (GLOBAL across all votes/candidates).
- gene contribution: if m is `entrez:*` gene (kind 3) -> {m.entity_id}; else `all_linked &= gene_ids non-empty` and genes += `entrez:`+g for g in gene_ids.
  After the posting: if `all_linked` and genes non-empty -> push as one gene set (a vote with any unlinked candidate gives no gene evidence).
- product contribution: only if NOT is_gene: candidates with kind==2 and anchor `uniprot:*` without '-' (primary accession, not isoform/chain) form a product set (pushed if non-empty).
Combine: `genes` = intersection of all gene sets; `proteins` = intersection of all product sets. (Source-GeneID `anchor="entrez:N"` field handling at `PC:171-173,231-235` is supported by the kernel and tests but the Python bundle never sets it for target 2 today; entrez ids vote through route 2 instead.)
- `product_conflict` = product sets exist but intersection empty; `gene_conflict` = gene sets exist but intersection empty.
- `protein` = the single entity id if `|proteins|==1 && !product_conflict && !quarantined`.
- `status` = conflict if either conflict; resolved if |genes|==1; ambiguous if >1; missing if 0.
- `gene_candidates` = union of all gene sets if gene_conflict else the intersection.
- accepted entities = `[the gene]` if `!quarantined && !product_conflict && status==resolved`; else `[protein]` if protein; else `[]`.
- outcome = QuarantinedCandidate if quarantined; ConflictingEvidence if conflict; Ambiguous if accepted empty; else Unique. Output row: `(input, outcome, accepted, |count|, protein, gene_candidates, status)`.
Consequences (docs/resolution.md, tests `build/tests/test_molecular_resolution.py`): gene ids/symbols/ensg/hgnc/enst/refseq-RNA/entrez/ramp_gene NEVER select a protein (even with one product); RNA types (gene_only) never assert protein;
protein identity requires a product vote (uniprot, uniprot-sec, uniprot_entry, ensp, refseq_protein, genbank — all collapsing to a primary accession); a protein with ambiguous gene links is still returned as the entity (`status=ambiguous`);
a known conflicting gene (e.g. protein TP53 + GeneID 672) -> `conflict`, entity falls back to the protein; unknown gene ids contribute nothing; any quarantined candidate anywhere -> abstain.
`reviewed` and `candidate_count` do not influence the outcome. `Match` then: gene resolved -> entity `entrez:N` (library gene_protein, protein_* fields describe the chosen product with its `gene_ids`),
`gene_mapping_status` and `gene_candidates` recorded; else the product entity or the retained asserted/native identity (`match.py:200-257`, `writer.py:48-79`).
`_retain_asserted_protein` (`match.py:382-444`): unmatched but explicitly reported single UniProt/ENSP/RefSeq-protein accession is kept as `protein_*` with `protein_mapping_status=reported`.

---------------------------------------------------------------------------------------------------------------
## 7. Compiler input tables (what an on-demand implementation must be able to query)

### 7.1 Hub parquet (frozen under `inputs/{hub}.parquet`; contract `hubs/schema.py:11-17`)
Columns: `source_type, source_id, hub_id, taxonomy_id, backend` (all strings; `file_row_number` read by claim stages).
Hubs: chemical `chebi, pubchem, chembl, hmdb, lipidmaps, swisslipids, bigg, metanetx, refmet, ramp`; gene_protein `uniprot, entrez, ramp_gene`
(`build_reference.py:22-34`; `uniprot` hub additionally carries secondary ACs `uniprot-sec`, `refseq_protein`/`genbank` from idmapping_selected, `prepare_sequence_mappings.py`; `entrez` hub from gene2ensembl carrying ensg/enst/ensp/refseq/refseq_protein).
Identity row for each hub record: `source_type=<hub>, source_id=hub_id=hub_id`. Chemical hub source_types include name/synonym/lipid_shorthand/systematic_name/formula/smiles/inchi/inchikey + xrefs.

### 7.2 Assigned-reference tables (output of `build_reference.py`; read by compilers via `projections.specification`, `projections.py:25-135`)
Schema per table (from the SQL in `build_reference.py`):
- `records-{hub}/records.parquet`: `hub, local_id, record_id('hub:local_id'), anchor_count, anchor, taxon, reviewed`  (compiler reads `records-uniprot`: anchor,taxon,anchor_count,reviewed).
- `records-{hub}/edges.parquet`: `source_record, target_hub, target_id, target_record, assertion_id, semantics ('gene_product'|'identity_xref')`.
- `records-{hub}/multiple_anchor_claims.parquet`: `record_id, anchor`.
- `{domain}-graph/{edges,endpoints,vertices,assessed_edges,keyless_edges}.parquet` (assessed_edges adds `source_anchor, source_anchor_count, target_anchor, target_anchor_count, target_exists`) — compile-time only; `gene_products` is derived from them.
- `{domain}-decisions/{vertices,boundaries,boundary_summary,assignments,edge_assessments}.parquet`, `{domain}-components/components.parquet` — entity assignment (not read by index compiler).
- `members-{hub}/members.parquet`: records cols + `entity_id, decision, component`  (read for `{hub}-names` join: `m.local_id = norm(hub_id)`).
- `{domain}-entities/entities.parquet`: `entity_id, kind('chemical'|'protein'|'gene'), primary_anchor, taxon, quarantined`.
- `{domain}-entities/identity_identifiers.parquet`: `namespace(hub|'inchikey'|'entrez'|'goslin'), identifier, entity_id, record_id`  (chemical: RaMP rows replaced by `finalize_source_reference`).
- `claims-{hub}/identifier_claims.parquet`: `record_id, entity_id, namespace, identifier, assertion_id, role` — all hub rows except source types smiles/inchi/formula/name/synonym/lipid_shorthand/systematic_name;
  role in attribute|native_identity|anchor|cross_reference|unverified_alias (not consumed by the compiler; used by `claims-entrez` (record_id) too).
- `goslin-identifiers/claims.parquet`: `hub, hub_id, record_id, name_type, source_name, name, goslin, level` (hubs swisslipids, lipidmaps, hmdb, chebi, refmet; most specific parsed level per record).
- `gene-products/gene_products.parquet`: `entrez_id, protein_entity_id, taxon, reviewed, projection_eligible`.
- `gene-resolution/entrez_identifiers.parquet`: `namespace='entrez', identifier, entity_id, admission, record_id`.
- `source-gene-resolution/identifiers.parquet`: `namespace='ramp_gene', identifier, entity_id, admission, record_id`.
- `inputs/{hub}.parquet` raw hub rows (for names: `source_type in (name, synonym, lipid_shorthand, systematic_name)`).
- `manifest.json` (`fingerprint`, `sources`, `stages`) and finalize's `source-chemical-mappings/{provenance,qc}.parquet` (audit).

### 7.3 Compiler projections actually read (name -> source -> key -> columns) (`projections.py`)
```
{chem,gp}-entities            entities.parquet                       entity_id
{chem,gp}-identity            identity_identifiers.parquet           identifier
{chem,gp}-identity-forward    identity_identifiers.parquet           entity_id   (entity_id,namespace,identifier)
claims-{hub}  (chebi,pubchem,chembl,hmdb,lipidmaps,swisslipids,bigg,metanetx,refmet,ramp,kegg)  identifier  (entity_id,namespace,identifier) WHERE ns in (kegg,cas,drugbank) OR (chebi/bigg own ns)
claims-uniprot                claims-uniprot/identifier_claims       sequence_key(identifier; refseq_protein/genbank version stripped)  (entity_id,namespace,identifier)
uniprot-ensg-forward          claims-uniprot WHERE ns='ensg'          entity_id
claims-entrez                 claims-entrez                           identifier (record_id,namespace,identifier) WHERE ns='ensg'
gene-resolution, gene-resolution-forward   entrez_identifiers
source-gene-resolution        source-gene-resolution/identifiers
records-uniprot               anchor,taxon,anchor_count,reviewed
gene-products (key entrez_id), gene-products-forward (key protein_entity_id)
gene_protein-gene-forward     gene_products JOIN claims-uniprot (ns in genesymbol,genesymbol-syn,hgnc,ensg) for proteins with exactly one gene
{domain}-{hub}-forward        claims-{hub}  (entity_id,namespace,identifier)   hubs: chemical minus pubchem,kegg; gene_protein uniprot,entrez,ramp_gene
{domain}-{hub}-names          inputs/{hub} JOIN members-{hub}        hubs except uniprot, entrez (names/synonyms/lipid_shorthand/systematic_name)
```
Gene-role builder additionally reads `gene_protein-entities/entities.parquet, claims-uniprot, claims-entrez, gene-products` directly (`B/gene_role_index.py:326-336`).

### 7.4 Minimal logical model an on-demand implementation must answer
Given a vote `(target, route, ns, identifier, scope)`:
1. posting = candidates from the assertion relations of 2.3 for that exact (ns, identifier) [+ expansions, isoform parent projection, version-stripped duplicates], filtered by 2.2 admission, restricted to entity taxon==scope if scoped (symbols require scope), exact-over-synonym for symbols, gene-collapse (2.6), then the >10 rule (4.1) and gene-role override (2.7).
2. candidate facts: kind, anchor, quarantined, gene_ids (reviewed optional).
3. entity record: label (3.2), alias list (3.3), taxon, anchor, gene_ids.
4. decision per 6.

---------------------------------------------------------------------------------------------------------------
## 8. Quirks / dead code worth knowing before reimplementing

- Scope applies to the whole observation when any symbol vote exists (1.5) — easy to miss; changes which entities uniprot/entrez ids can match (taxon-filtered).
- Missing/removed keys never veto (chemical non-anchor; all gene/protein votes); only a found-but-different posting vetoes. Hence >10-removed identifiers silently disappear from evidence.
- Quarantine is contagious in the molecular kernel (any quarantined candidate in any posting of any vote), but only on the final intersection for chemicals.
- Chemical anchor override: a full InChIKey (explicit or SMILES-derived) discards ALL other chemical votes.
- `products` flag, `reviewed` bit, `candidate_count`, `primary` column, `projection_eligible`, route-3/`product` LMDB, ensg fallback, `kegg_gene`, `inchi`/`cas`/`drugbank`/`kegg` native keys, kernel `anchor="entrez:*"` on votes: unused or unreachable in current output.
- `replay_resources.VSCHEMA` (`replay_resources.py:61-70`) has no `gene_only`: persisted replay vote tables lose the RNA gene-only flag (only live `observation_bundle` rows carry it).
- Runtime miss defaults for `gene`/`products` (`index.py:96-100`) differ from compiled gene namespace set (enst/refseq) but are harmless.
- Gene-role override covers only target 2 route 1 and ns {ensg, hgnc, genesymbol, genesymbol-syn, enst, refseq}; hgnc/ensg/etc. otherwise fall back to the collapsed base posting; empty override postings are authoritative.
- Scoped taxon in key is the CANDIDATE entity taxon, not a property of the identifier claim; entities with NULL taxon (e.g. RaMP/anchorless) are visible only unscoped.
- Isoform handling is spread across: votes (isoform + base vote), assertions (isoform entity re-pointed to parent), kernel (anchors with '-' never products).
