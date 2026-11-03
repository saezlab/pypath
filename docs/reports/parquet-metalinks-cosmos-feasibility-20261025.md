# Building MetalinksDB and COSMOS from Parquet: feasibility

2026-10-25

**Short answer: feasible.** Both products take seconds to compute with DuckDB over a release's Parquet files.

The work is porting what they read from PostgreSQL's derived layers, mainly the interaction classification and COSMOS's finishing steps. Compute is not the obstacle.

## How they are built today

Both are built after a release is loaded into PostgreSQL.

The dev5 load of 2026.10.19-b took **3.0 h** and about 210 GB of PostgreSQL, from 16 GB of Parquet. The parts that matter here:

| Step | Time |
|---|---|
| Interaction projection: the record table, and reaction headers merged across sources | 26 min |
| Other derived layers: labels, classes, metabolic domain, bitmaps | 59 min |
| COSMOS product | 79 s |
| MetalinksDB | 0 s (a registry row only) |

**MetalinksDB** is a recipe stored in `network_registry`, not a table. A consumer applies it to `interaction_fact_resource`. Nothing in our Parquet API, engine or client executes presets.

**COSMOS** is the `cosmos_edge` table. It reads:

- `interaction` and `interaction_party`: reaction stars, merged across sources on their chemistry;
- direction and compartment annotations;
- the exclusion metabolites product;
- for labels, published identifiers plus the utils PostgreSQL `id_mapping`, which lives on another server.

## What the Parquet release already has

`relation.parquet` has, per statement:

- endpoint keys, already resolved and shared across resources;
- endpoint types, predicate and taxon;
- the effect qualifiers: aspect, direction and mechanism;
- the contributing sources.

`relation_evidence` records the dataset, which is what ChEMBL's mechanism-of-action flag comes from. `entity_annotation` holds reaction directions and compartments. `entity_identifier` holds published aliases.

Transport pairs are already projected per resource at build time (`omnipath_build/projections.py`).

## Prototypes

Both ran on nicesrv over the 2026.10.19-b pins plus MetaChatDB, using DuckDB with 2 threads and 3–4 GB. The scripts are in `/root/data/scratch/build-20261025/`: `metalinks_parquet2.py` and `cosmos_parquet.py`.

**MetalinksDB.** The preset's own composition, read from `METALINKSDB`, ran as one query in **2.8 s**.

- **Components:** ChEMBL mechanism-of-action pairs; the curated sources; the Rhea, Recon3D and MetAtlas transport pairs.
- **Gate:** one end is a chemical, the partner type is allowed, and the predicate is allowed.
- **Collapse:** on the two endpoints.
- **Result:** 891,120 source records became 827,378 pairs, 11,288 of them from more than one source. All 12 contributing sources delivered rows.

One substitution was needed. The reaction component filters on `interaction_class IN (transport, ligand_receptor)`. In Parquet, `interaction_class` only says whether a relation is directed. The prototype used the transport aspect qualifier instead.

**COSMOS core.** The reaction stars of KEGG, Rhea, MetAtlas and Recon3D, merged across sources on their reactant and product sets, then binarised around each enzyme. It took **1.5 s**.

| | PostgreSQL (dev5) | Parquet prototype |
|---|---|---|
| Reactions (headers) | 34,411 | 34,492 |
| Reactions without a catalyst | 19,631 | 20,099 |
| Gene nodes (enzyme × reaction) | 385,901 | 358,794 |
| Edges | 2,671,581 | 2,536,590 (forward only) |
| Metabolite nodes | 28,434 | 14,881 (no compartments) |

The prototype leaves out what explains most of the gaps:

- compartment-specific metabolite nodes;
- reverse halves of reversible reactions;
- edges for reactions without a known enzyme;
- connectors;
- exclusion metabolites;
- label translation.

## What would have to be ported

1. **Interaction classification.** PostgreSQL's classes (`relational/classify/interaction_class.*` plus the class-evidence step of the derive) don't exist in Parquet. Either port the rules to a DuckDB step at release level, or rewrite the recipe on qualifiers, as the prototype did.
2. **MetalinksDB attributes.** The `intercell` layer (role and location terms from resource annotations), references, evidence and license scope. These are joins on tables Parquet already has.
3. **COSMOS finishing steps:**
   - reverse halves, using the entity direction annotation;
   - compartment nodes, using the entity compartment annotation;
   - orphan edges for reactions without a known enzyme;
   - connectors;
   - label translation.

   `project_edges.sql` and `project_connectors.sql` are about 430 lines of PostgreSQL SQL. Most of it maps directly to DuckDB: arrays, `bool_or` and window functions. The rest is psycopg plumbing: unlogged staging tables, `ANALYZE`, advisory locks.
4. **Exclusion metabolites**, which COSMOS reads. This is about 900 lines of Python and SQL over the ChEBI hierarchy and structure keys. All of its inputs are in Parquet.
5. **Label translation.** Published aliases cover most labels. The utils `id_mapping` sits in a separate PostgreSQL either way: keep calling it, or use the identity library's hub indexes, which are already Parquet.
6. **Where the output lives.** Today a release has only per-resource versions, so this would be a new release-level product. For example, `releases/<version>/products/{metalinksdb,cosmos_edge}.parquet`, pinned by digest in the release manifest.

   MetalinksDB is fast enough to be a query in the engine (an explorer preset, or a prebuilt client query) instead of a file. COSMOS fits a file better.

## Rough effort

- **MetalinksDB: small.** One query, the classification decision, the attribute joins, and a parity check against dev5.
- **COSMOS: medium.** The edge SQL, connectors, labels and exclusion metabolites. The 26-minute interaction projection doesn't need porting in full, only the reaction headers.

## Open questions

- **Which side keeps the products.** dev5 consumers read them from PostgreSQL. Keeping both implementations conflicts with having no redundant code, so pick one owner per product.
- **Merge semantics.** COSMOS's PostgreSQL headers are keyed on chemistry, and the prototype agrees within 0.2%. `docs/parquet-subsets-milestone.md` says the reverse ("does not reproduce the legacy merge"), and that doc looks outdated.
- **Parity.** Before switching, compare full outputs pair for pair against dev5. The prototypes only compared counts.

## Side finding (fixed the same day)

Release 2026.10.25 served TCDB (2026.10.7) and Rhea (2026.10.9) builds that predate the transport projection, so they had no transport pairs.

Release **2026.10.25.1** fixes this. It pins the newest complete build of all 173 resources and is now the site's default. 42 resources moved to newer builds, and none lost entities or relations. TCDB went from 0 to 44,721 relations.

It still differs from what dev5 has loaded (2026.10.19-b): 17 resources are newer 2026.10.16 builds on the current identity library.
