# Glossary

The terms we use in code, documentation and discussion. The key terms come
first; read those before the other pages. When two words could mean the same
thing, this page says which one we use.

## Key terms

### Resource
One upstream data provider as built by us, for example SIGNOR or ChEBI. A
resource can contain several datasets.

### Resource version
One immutable, numbered directory of Parquets for a resource, for example
`resources/signor/1.0.0/`. Independent of the provider's own release number.

### Release
A JSON manifest pinning one version of each included resource, plus the digest
of each build manifest, for example `2026.09`. Immutable once published.

### Latest
The API default: the highest numeric version of each resource, including newly
added resources. Changes as builds are published.

### Entity
One row of the `entity` table: a typed thing such as a protein, gene, chemical,
complex or ontology class. Identified by its [entity key](#entity-key).

### Entity type
The Biolink type reported by the source, for example `protein` or `gene`.
Resolution never changes it.

### Reference entity key
The shared grouping reference of an entity, normally `entrez:<GeneID>`. A
protein, a gene and an RNA of the same gene share it while keeping separate
entity keys. Also called *gene reference*.

### Molecular form
The occurrence-level description of what exactly was observed: product
pointer, isoform, sequence identifiers, modifications and variants. Missing
fields mean unspecified.

### Statement
A subject–predicate–object claim with its qualifiers and evidence, stored as one
row of the `relation` table. Its `statement_kind` is `relation` or `ontology`.

### Evidence
One source occurrence supporting an entity or statement:
`{source, dataset, row_id, upstream_id, annotations, molecular forms}`.
Evidence is a multiset, so identical occurrences are kept.

### Reference library
What matching reads during resource builds: the hub indexes plus the identity
library (`omnipath-identity-v2`) of identity decisions. Built offline from hubs;
also called the identity layer.

### Anchor
An identifier that defines identity on its own: the full InChIKey for
chemicals, the Goslin name for lipids without a structure, the primary UniProt
accession for proteins, the NCBI GeneID for genes and the Rhea master reaction
for reactions. Reference entities are built around anchors.

### Abstain
Resolution's choice to leave an observation unresolved because its evidence is
missing, conflicting or ambiguous. Deliberate behaviour, not a failure.

### Parquet contract
The published tables of a resource version and their schemas
(`omnipath_core/schema.py`), currently `serving_schema_version = 5`. The only
interface between build, serving and PostgreSQL.

## Other terms

### Ambiguous
More than one candidate entity survives resolution, so matching abstains.

### Annotation
A `{term, value, quantity, source, dataset}` attribute attached to an entity,
statement or evidence occurrence. Terms are Biolink slots or ontology CURIEs.
On statements, `scope` says whether it describes the relation, subject or object.

### Base checkpoint
The point in a PostgreSQL load where all base tables, constraints and indexes
are committed and validated. Later phases resume from here.

### Biolink
The [Biolink Model](https://biolink.github.io/biolink-model/) supplies our entity
types, predicates, qualifiers and annotation slots.

### Build manifest
`build_manifest.json` in each resource version: versions, row counts, sizes and
SHA-256 checksums of the three Parquets, plus build provenance.

### Candidate
An entity in the reference library that a lookup key points to. Matching
intersects the candidates of all votes.

### Component
A connected group of anchorless hub records without a cross-reference to an
anchored record, linked by cross-references. It becomes one entity only if it
holds at most one record per hub; otherwise every record stays its own entity.

### Connectivity
The first 14-character block of an InChIKey, which encodes the molecular
skeleton without stereochemistry or charge. Used to group chemicals at query time.

### Entity key
SHA-256 of `(entity_type, namespace, identifier)`. Two rows with the same
identifier but different types have different keys.

### Evidence payload
The raw source record as JSON, stored in `evidence_payloads.parquet` and linked
to evidence by `row_id`.

### Fallback
The result when resolution cannot establish a reference: the observation keeps
its own native identifier (or reported product), with annotations explaining why.
Also called a *native identity*.

### Gene reference
See [reference entity key](#reference-entity-key).

### Gene mapping status
`resolved`, `ambiguous`, `conflict` or `missing`: the outcome of gene-level
resolution for one occurrence. Anything but `resolved` is stored as the
`omnipath:gene_mapping_status` annotation.

### Generation
One complete, immutable identity library, named by the fingerprint of its hub
indexes and rules code. Each resource build pins one and records it.

### Hub
A normalized export of a reference database (UniProt, NCBI Gene, ChEBI,
PubChem, …): rows of *record carries identifier*, used to build the reference library.

### inputs_v2
The pypath module that parses source resources into typed entities and
relations. It is the input of every resource build.

### Native identity
See [fallback](#fallback).

### Occurrence
A single appearance of an entity or statement in a source row. Each occurrence
becomes one [evidence](#evidence) item.

### Phase
A durable, separately committed step of a PostgreSQL load (the base, a
derivation or a product), recorded in `parquet_phase`.

### Policy
The per-entity-type matching rules: which library is used, which namespaces
vote, which need a taxon and how labels are chosen.

### Predicate
The Biolink relationship of a statement, for example `affects` or `physically_interacts_with`.

### Product
A gene product with its own reusable entity row: a primary UniProt protein or a
specifically reported transcript. Referenced from molecular forms.

### Qualifier
A Biolink annotation that changes what a statement means (for example direction
or aspect). Qualifiers are part of the relation key.

### Quarantined
A hub record that claims more than one anchor. It stays its own entity and
never joins others; in matching, its identifier votes for every anchor it
claims, so the other identifiers decide.

### Relation
In Parquet, a row of the `relation` table (a [statement](#statement)). In
PostgreSQL, `relation` is main's graph triple (subject, predicate, object); several
qualified statements can share one triple.

### Serving schema version
The version of the Parquet record contract (currently 5). Separate from the
build-manifest `schema_version` (1).

### Subset
A curated product built inside a loaded PostgreSQL schema: MetSigDB, network
views or COSMOS.

### Taxon
The NCBI Taxonomy ID of an entity or statement, for example `9606` for human.
Empty when assertions conflict.

### Vote
A normalized identifier from an observation that is allowed to take part in
matching. Names never vote; gene symbols vote only with a taxon.
