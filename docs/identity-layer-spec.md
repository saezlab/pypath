# Identity layer and on-demand resolution: implementation spec (phase 1)

Status: implemented on `parquet-migration` (8 October 2026), including reactions
from phase 2; miRNA resolution and the omnipath-utils pair export are not done.
The current behaviour is described in
[entity resolution](../core_documentation/resolution.md). Written 6 October 2026
as the working spec for branch `resolution-v2`. Rules come from the
proposal "Entity resolution: rules and implementation"; current semantics are mapped
in the compiler brief (section numbers below refer to it as "brief §").

Phase 1 covers chemicals and genes/proteins from the existing hub Parquet files
(`/root/data/hubs/20260909` on nicesrv). Phase 2 adds miRNA, reactions and the
omnipath-utils pair export. Nothing in phase 1 may change the published Parquet
contract (`omnipath_core/schema.py`).

## 1. What replaces what

| Today | Phase 1 |
| --- | --- |
| `build_reference.py` graph, components, decisions, claims | (a) **hub index**, per hub, independent of the rules; (b) **identity decisions**, small, dependent on the rules |
| LMDB compiler (`direct_compact_index`, `full_index*`, `candidate_limit`, `gene_role_index`) | none |
| `FullRuntime` (LMDB) | `omnipath_resolver/identity_runtime.py: IdentityRuntime`, same `resolve()` contract |
| entity records precomputed for 304M entities | built on first use, cached per identity fingerprint |

Design principle: nothing is materialized for the whole universe that follows trivially from a
record. 295,964,650 hub records exist; only about 0.9M (0.31%) need a real decision. Every other
record's entity is its own anchor (`inchikey:K`, `uniprot:X`) or its gene id (`entrez:N`). Hub
data changes rebuild only that hub's index; rule changes rebuild only the decisions.

The Rust kernels, `LibraryMatcher`, `observations.py` key encoding and the writer stay.
`LibraryMatcher` picks the runtime from the manifest `format`.

## 2. Rules to implement (phase 1)

Identity (build time):
1. Anchors: chemical = valid full InChIKey (same regex and empty-key exclusions as today);
   protein = UniProt record whose own primary accession it is (same as today);
   gene = every NCBI Gene hub record is the gene entity `entrez:N` (same as today);
   **lipid name** = for a chemical record with no InChIKey, its most specific Goslin name at
   species level or finer becomes its anchor, entity id `goslin:<level>:<name>` (same string
   format as today's goslin ids). Class/category level is never an anchor.
2. A record claiming two or more anchors of the kind (two InChIKeys, or two lipid names at its
   most specific level) is **quarantined**: its own entity, flagged, never connects anything.
3. **Attach (one step):** a record without an anchor joins an entity when its own identity
   cross-references (edges with `semantics='identity_xref'`, either direction, never
   `gene_product`) reach records whose anchors are all the same single anchor and none of those
   records is quarantined. Chains are never followed.
4. **Group:** records that have no anchor, were not attached, are not quarantined and have no
   identity edge to any anchored record are connected through their identity edges. A connected
   group becomes one entity only if it holds at most one record per hub; otherwise every record
   stays its own entity. Group entity id = record of the most preferred hub
   (chebi, lipidmaps, swisslipids, hmdb, chembl, pubchem, kegg, metanetx, bigg, refmet, ramp,
   ramp_gene), lowest `local_id` within that hub.
5. A full-structure Goslin name resolves to an InChIKey only when every record carrying that name
   and an InChIKey agrees on exactly one InChIKey; otherwise the name is not a lookup key for a
   structure entity.
6. Gene–product links (`gene_products`): same derivation as today (brief §3.6).

Lookup (runtime):
7. No candidate cutoff. Postings keep every candidate; more than one surviving candidate means
   the kernel abstains (Ambiguous). Keep the kernel's hard safety limits.
8. Precedence within one posting, applied before the kernel: native identity rows beat
   fallback rows (today's `chemical_fallback`); a primary-accession identity row beats
   `uniprot-sec` claim rows under ns `uniprot`; exact symbols beat `symbol_synonym` rows within
   one taxon. So a secondary accession resolves only when it maps to exactly one entry.
9. Taxon scope exactly as today (brief §1.5–1.6): scoped only when the observation has a symbol
   vote; symbols only scoped; scope = candidate entity taxon.
10. Gene-level namespaces (hgnc, ensg, enst, refseq RNA, genesymbol, genesymbol-syn, route 2)
    post gene entities, built like today's gene-role component (brief §2.7) but without the
    10-gene cap. Product namespaces post protein entities with their `gene_ids`.
11. **A SMILES-derived InChIKey is not a primary anchor.** It votes as an ordinary identity vote
    (`anchor=""`). A stated InChIKey keeps the decide-alone priority.

Entity records and labels:
12. Record = `{entity_id, kind, anchor, taxon, label, identifiers, gene_ids}` exactly as
    `FullRuntime.record` returns (brief §3), built from member records' hub rows.
13. Labels, first available, deterministic, ties to shorter then alphabetical:
    gene: NCBI Gene symbol of the entrez record › `entrez:N` local part;
    protein: primary gene name of the UniProt entry (`genesymbol` claim, not `genesymbol-syn`) ›
    entry name › accession;
    chemical: ChEBI name › HMDB name › ChEMBL name › PubChem name › any other hub `name` ›
    IUPAC/systematic name; skip values matching the InChIKey regex and values over 80 chars;
    lipid (entity with a Goslin anchor, or a structure whose records carry a Goslin name):
    most specific Goslin shorthand;
    otherwise the entity id local part. Never a database ID formatted as a label
    (`CHEBI:…`, `CID:…`) while a name exists.

## 3. On-disk layout

### 3a. Hub index (per hub; rebuilt only when that hub's Parquet changes)

`/root/data/hubindex/<hub>/<sha12 of hub parquet>/`, built by
`omnipath-build build-hub-index --hub NAME --hubs-dir DIR --output-root DIR`.

```
manifest.json   hub, input {sha256, bytes}, code_sha256 (hash of the hub-index code),
                counts, timings
records.parquet one row per hub record, sorted by local_id, ~64k-row groups:
                local_id, record_id ('<hub>:<local_id>'), taxon (NULL when unknown),
                anchor ('inchikey:K' | 'uniprot:X' | NULL; set only when anchor_count = 1),
                anchor_count, reviewed
by_id/part=XX/  lookup rows, XX = md5(identifier)[:2], sorted by (ns, identifier, local_id):
                ns, identifier, local_id, tag, taxon, anchor, anchor_count
                (record columns denormalized from records.parquet; they change only with the hub)
by_record/part=XX/ every normalized hub row incl. names, XX = md5(record_id)[:2],
                sorted by local_id: local_id, source_type, value
```

- Normalization exactly as today (`normalize_id_sql`, version stripping, refseq typing; brief §1.7).
- `tag`: `native` (the record's own id), `claim` (a cross-reference or alias), `secondary`
  (a `uniprot-sec` value repeated under ns `uniprot`), `symbol_synonym` (a `genesymbol-syn` value
  repeated under ns `genesymbol`), `version_stripped` (refseq/genbank without version). These
  are the data-level expansions of brief §2.3; which rows become candidates is decided at runtime.
- Names, synonyms, smiles, inchi, formula are only in `by_record`, never in `by_id`.
- kv store (`kv/`, `omnipath-hub-kv-v3`): key `part + ns + 0x00 + identifier`; for ns `genesymbol`
  and `genesymbol-syn` (symbols exist only taxon-scoped) `part + ns + 0x00 + identifier + 0x00 +
  taxon`, taxon empty for a record without one (`0`, `None` and blanks count as none). A scoped
  lookup reads the keys of its taxon and the empty one, never the head of another species'
  record. Format v2 has no taxon in the key; the runtime reads both (v2: every species, then
  the record taxon filter).
- Chemical hubs also get `goslin` rows (ns `goslin`, identifier `<level>:<name>`, the record's
  most specific parsed name at species level or finer) in both `by_id` and `by_record`, parsed
  with the persistent Goslin cache.

### 3b. Identity decisions (rule-dependent, small)

`/root/data/identity/<fingerprint>/`, built by
`omnipath-build build-identity --hub-index-root DIR --output-dir DIR`.
fingerprint = sha256 of (each hub index manifest sha) + rules code sha.

```
manifest.json            format "omnipath-identity-v2", fingerprint, hub_indexes {hub: dir},
                         rules_sha256, counts, timings
exceptions.parquet       sorted by record_id: record_id, entity_id, decision, quarantined
                         every record whose entity is NOT its own anchor / entrez id:
                         lipid_name, attached, grouped, structureless, ambiguous_native,
                         quarantined, ramp_gene source-gene mappings (decision
                         explicit_source_gene / source_gene_only, brief §2.3)
exception_members.parquet sorted by entity_id: entity_id, record_id (same rows, other order)
entities_extra.parquet   sorted by entity_id: entity_id, kind, taxon, quarantined, preferred_record
                         only for entity ids not derivable from a record (goslin:…, group and
                         structureless record ids, quarantined record ids)
gene_products_by_protein.parquet  sorted by protein_entity_id: protein_entity_id, entrez_id, taxon
gene_products_by_gene.parquet     the same rows sorted by entrez_id
lipid_structures.parquet goslin full-structure name -> inchikey, only names with exactly one
                         InChIKey across records (rule 5)
```

Entity of a record = its exception if present, else its `anchor` (anchor_count = 1), else
`entrez:<local_id>` for an entrez record, else its record_id. Kind from the entity id prefix
(`inchikey:` chemical, `uniprot:` protein, `entrez:` gene) or `entities_extra`.

## 4. Runtime contract

`IdentityRuntime(path, cache_dir=None)` opens a v2 identity directory and the hub indexes it names.

- `lookup_many(keys)` per batch: decode keys (brief §1.1); read `by_id` rows for the keys from
  the hubs of the key's domain (target 1: chemical hubs; target 2: uniprot, entrez, ramp_gene),
  only the partitions the identifiers hash to; map records to entities (rule above, joining
  `exceptions` for the batch's records); apply admission (which hub/ns/tag rows are lookup keys,
  brief §2.2–2.3), isoform → parent projection, precedence (rule 8), scope (rule 9, entity
  taxon = defining record's taxon or `entities_extra`/gene_products taxon), gene-level postings
  (rule 10, via `gene_products_by_protein`: proteins linking to exactly one gene, taxa
  compatible); return `{"candidates": [...], "gene": bool, "products": False}`.
- candidate fact = `[num, entity_id, kind(1/2/3), anchor|None, quarantined, reviewed, gene_ids]`;
  `num` = sha256(entity_id)[:8] big-endian >> 1; sorted by num; collision in a batch raises.
- `record_many(entity_ids)`: members = for an anchored entity the records whose `by_id` row
  (ns `inchikey`/`uniprot`, identifier K/X) has that anchor and no conflicting exception, plus
  `exception_members`; for `entrez:N` the entrez record plus the genesymbol, genesymbol-syn,
  hgnc and ensg rows of proteins linked to exactly that one gene; rows from `by_record`;
  label per rule 13. Cached in `<cache_dir>/<fingerprint>/entities.sqlite` (WAL).
- `resolve(queries, votes)` = `FullRuntime.resolve` with `lookup_many`/`record_many`.

## 5. Validation

- Gold set: `packages/build/tests/identity/` fixtures: tiny hub Parquet files with one case per
  rule above, expected entity per observation.
- Regression on nicesrv: extract observations once per resource (pre-resolution votes), resolve
  them with `FullRuntime` on the Oct 4 library and with `IdentityRuntime`; report counts per
  outcome and every differing observation, classified by the rule that explains it.

## 6. nicesrv conventions

- Code: `/root/projects/resolution/<checkout>`; push branches to the private bare repo
  `ssh://nicesrv/root/git/omnipath-build.git` (never to GitHub without asking).
- Data: hubs `/root/data/hubs/20260909`, pypath cache `PYPATH_DOWNLOAD_DATADIR=/root/data/pypath-cache`,
  hub indexes `/root/data/hubindex/`, identity decisions `/root/data/identity/`, regression output `/root/data/regression/`.
- Old library (oracle): `/root/projects/omnipath-releases/20261004-molecular-forms/reference/library`.
- 8 CPUs, 15 GiB RAM, explore.omnipathdb.org runs on the same host: run heavy jobs with
  `systemd-run --scope -p MemoryMax=9G nice -n 10 …`, DuckDB `memory_limit` ≤ 7GB, ≤ 6 threads.

## 7. Postponed work

Decided 6 October 2026: these stay out of the first rebuild.

- **miRNA resolution** (proposal section 1): miRBase precursor (`MI…`) and mature (`MIMAT…`)
  anchors, official miRBase names as votes, a gene reference from the precursor's NCBI Gene.
  Today only the miRBase resource emits miRNAs (48,885 mature, 38,589 precursor), already as
  accessions.
- **Identifier tables from omnipath-utils**: utils lacks RaMP compound and gene ids, RefMet's
  and BiGG's own exports, lipid shorthand/systematic names, SEED and MetaCyc ids, and formulas;
  `uniprot-sec`, Entrez–Ensembl coverage across organisms and per-record provenance need a real
  export to confirm. Switch hub by hub once utils covers a hub completely.

## 8. Reactions (implemented 6 October 2026)

- Hubs `rhea` (one record per master reaction: its directional ids as `rhea` claims, Rhea's own
  cross-references to KEGG, Reactome, MetaCyc, EcoCyc and M-CSA, EC numbers as attributes, the
  master's `DEFINITION` as `name`) and `metanetx_reaction` (one record per MNXR with its Rhea,
  KEGG, MetaCyc, BiGG, SEED, VMH and SABIO-RK ids; Rhea ids written as their masters).
- A Rhea record is anchored on its master (`rhea:<id>`). A MetaNetX reaction attaches one step to
  the single Rhea master it names; without one it is its own entity
  (`metanetx_reaction:<MNXR>`); naming several it is ambiguous and evidence for each.
- Library `reaction` for `molecular_activity`, key target 3. For one id, Rhea's own
  cross-reference beats MetaNetX's (native over fallback). Reactome ids lose prefix and version.
- Decisions use the existing chemical kernel path (intersection, quarantine, unique/ambiguous):
  the runtime passes reaction queries and entities to it as chemical ones; no structure anchors.
- Label: Rhea's equation, else the source's name.
