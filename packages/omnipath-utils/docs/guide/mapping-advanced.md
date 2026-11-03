# Advanced ID Translation

This page covers the internal workings of the translation system for
power users and developers. For the basic API reference, see
[ID Translation](mapping.md).

## Controlling translation behavior

### Raw mode: bypass all special handling

Every translation function accepts `raw=True`. In raw mode, the
translator performs a single direct table lookup -- no case fallbacks,
no chain translation, no UniProt cleanup. This is useful when you need
maximum speed and know your identifiers are already in the exact form
stored in the mapping table.

=== "Python"

    ```python
    from omnipath_utils.mapping import map_name, translate

    # Raw: exact match only
    map_name('TP53', 'genesymbol', 'uniprot', raw=True)
    # {'P04637'}

    map_name('tp53', 'genesymbol', 'uniprot', raw=True)
    # set()  -- lowercase, no fallback

    # Batch raw translation
    translate(
        ['TP53', 'EGFR', 'tp53'],
        'genesymbol', 'uniprot',
        raw=True,
    )
    # {'TP53': {'P04637'}, 'EGFR': {'P00533'}, 'tp53': set()}
    ```

=== "REST API"

    ```bash
    curl "https://omnipathdb.org/mapping/translate?\
    identifiers=TP53,tp53&\
    id_type=genesymbol&\
    target_id_type=uniprot&\
    raw=true"
    ```

    The response `meta` object includes `"raw": true` to confirm the
    parameter was applied.

=== "DataFrame"

    ```python
    from omnipath_utils.mapping import translate_column

    result = translate_column(
        df, 'gene', 'genesymbol', 'uniprot',
        raw=True,
    )
    ```

What raw mode skips:

| Feature | Normal mode | Raw mode |
|---------|-------------|----------|
| Gene symbol case fallbacks (UPPER, Capitalize) | Yes | No |
| Gene symbol synonym lookup | Yes | No |
| Chain translation (source -> uniprot -> target) | Yes | No |
| UniProt cleanup (secondary -> primary, TrEMBL -> SwissProt) | Yes | No |
| RefSeq version stripping | Yes | No |
| Ensembl version stripping | Yes | No |
| miRNA reciprocal fallback | Yes | No |
| CURIE prefix stripping | Yes | No |
| Reverse table scan | Yes | No |

### Explicit backend selection

By default, the mapper auto-selects a backend based on which data
sources support the requested ID type pair. You can force a specific
backend with the `backend=` parameter.

=== "Python"

    ```python
    from omnipath_utils.mapping import map_name, translate

    # Force BioMart for Ensembl-centric lookups
    map_name('TP53', 'genesymbol', 'ensg', backend='biomart')

    # Force UniProt REST for protein IDs
    translate(
        ['P04637', 'P00533'],
        'uniprot', 'genesymbol',
        backend='uniprot',
    )
    ```

=== "REST API"

    ```bash
    curl "https://omnipathdb.org/mapping/translate?\
    identifiers=TP53&\
    id_type=genesymbol&\
    target_id_type=uniprot&\
    backend=biomart"
    ```

    !!! note
        The REST API uses the database backend for lookups. The `backend`
        parameter is recorded in the response metadata but does not
        change the query behavior in database mode. It is primarily
        meaningful for the Python API.

Available backends:

| Backend | Data source | When to force it |
|---------|-------------|------------------|
| `uniprot` | UniProt REST API | Default for protein IDs; most comprehensive protein cross-refs |
| `uniprot_ftp` | UniProt FTP idmapping files | Bulk download; faster for large-scale builds |
| `uploadlists` | UniProt ID Mapping batch service | Specific, targeted ID sets |
| `biomart` | Ensembl BioMart | When you need fresh Ensembl data or Ensembl-specific ID types |
| `mirbase` | miRBase | miRNA names and accessions |
| `unichem` | UniChem (EMBL-EBI) | Chemical compound cross-references |
| `ramp` | RaMP-DB | Metabolite cross-references and synonym mappings |
| `hmdb` | HMDB | HMDB-centric metabolite mappings |
| `metanetx` | MNXref | Pairwise metabolite ID translation via MetaNetX bridge (bigg↔chebi, kegg↔chebi, hmdb↔chebi, lipidmaps↔chebi, swisslipids↔chebi, metanetx↔*) |
| `bigg` | BiGG Models | BiGG metabolite mappings (bigg↔chebi, bigg↔hmdb, bigg↔kegg, bigg↔metanetx) |

When `backend=` is specified, the mapper skips its cached table and
forces a reload from the requested backend. This is useful when the
auto-selected backend returned incomplete data and you want to try a
different source.

### Querying SwissProt, TrEMBL, and synonyms directly

The `uniprot` target type runs the full cleanup pipeline (secondary AC
resolution, TrEMBL-to-SwissProt preference, proteome filtering, format
validation). To bypass this and query specific subsets of UniProt, use
the dedicated target types:

```python
from omnipath_utils.mapping import map_name

# Default: prefers SwissProt, runs full cleanup
map_name('TP53', 'genesymbol', 'uniprot')
# {'P04637'}  -- SwissProt entry

# Only reviewed (SwissProt) entries, no cleanup
map_name('TP53', 'genesymbol', 'swissprot')
# {'P04637'}

# Only unreviewed (TrEMBL) entries, no cleanup
map_name('TP53', 'genesymbol', 'trembl')
# TrEMBL entries for TP53, if any exist

# Gene symbol synonyms
map_name('p53', 'genesymbol-syn', 'uniprot')
# Looks up 'p53' as a synonym; may find TP53's UniProt AC
```

Target type behavior:

| Target type | Source data | Cleanup pipeline | Typical use |
|-------------|-----------|-----------------|-------------|
| `uniprot` | SwissProt + TrEMBL | Yes (full 4-step) | Default; production use |
| `swissprot` | SwissProt only | No | When you specifically need reviewed entries |
| `trembl` | TrEMBL only | No | When you need unreviewed entries |
| `genesymbol-syn` | Gene symbol synonym table | No (source type) | Looking up old/alternative gene names |

## The translation pipeline in detail

### Step-by-step walkthrough

Here is what happens when you call:

```python
map_name('tp53', 'genesymbol', 'uniprot')
```

**Step 1: Alias resolution.**
`IdTypeRegistry.resolve()` normalizes the type names. `'genesymbol'`
is already canonical. Variants like `'GeneSymbol'`, `'gene_symbol'`,
or `'genesymbol_syn'` would be resolved to their canonical forms.

**Step 2: Direct lookup.**
The mapper looks up `'tp53'` in the `genesymbol -> uniprot` table.
Tables are case-sensitive, so `'tp53'` is not found (the table has
`'TP53'`). Result: miss.

**Step 3: Gene symbol fallbacks.**
Since `id_type` is `'genesymbol'`, the fallback chain runs:

- `'tp53'.upper()` = `'TP53'` -- lookup finds `{'P04637', 'A0A024R1R8', ...}`.
  Match found; remaining fallbacks are skipped.

If UPPER had failed, the system would try:

- `'tp53'.capitalize()` = `'Tp53'` (for rodent symbols)
- Synonym table lookup for `'tp53'` and `'TP53'`
- Append "1": `'tp531'` (non-strict mode only)

**Step 4: UniProt cleanup.**
Since the target is `'uniprot'`, the cleanup pipeline runs on the
raw result set `{'P04637', 'A0A024R1R8', ...}`:

1. **Secondary -> primary:** Each AC is checked against the
   `uniprot-sec -> uniprot-pri` table. If it is a secondary AC,
   it is replaced with the current primary AC.

2. **TrEMBL -> SwissProt:** `A0A024R1R8` is not in the SwissProt
   reference list. The pipeline looks up its gene symbol via the
   `trembl -> genesymbol` table, gets `'TP53'`, then finds the
   SwissProt entry for `'TP53'`: `'P04637'`. The TrEMBL AC is replaced.

3. **Proteome filter:** The result set is intersected with the human
   proteome reference list. `P04637` is present; stale or cross-organism
   ACs are removed.

4. **Format validation:** Each AC is checked against the UniProt AC
   regex. Invalid strings are discarded.

**Step 5: Result.**
`{'P04637'}` -- just the SwissProt entry.

### Vectorization

The `translate()` and `translate_column()` functions use
`translate_core()` internally, which implements a two-pass strategy:

1. **First pass (batch):** All identifiers are looked up in the
   mapping table in a single sweep. This is an O(1) dict lookup per ID.
   Identifiers that get a direct hit are collected immediately.

2. **Second pass (per-ID fallback):** Only identifiers that missed
   in the first pass go through the full `map_name()` pipeline with
   all its fallback strategies.

This means: if your mapping table covers 99% of your identifiers,
only 1% go through the slower per-ID fallback path. For typical gene
symbol-to-UniProt translations with properly cased input, the first
pass handles nearly everything.

```python
# Efficient: translate_core handles the batch
from omnipath_utils.mapping import translate
result = translate(gene_list, 'genesymbol', 'uniprot')

# Less efficient: each call goes through the full pipeline
from omnipath_utils.mapping import map_name
result = {g: map_name(g, 'genesymbol', 'uniprot') for g in gene_list}
```

When the target type is `uniprot`, cleanup runs on the batch results
from the first pass as well (not just the fallback results), ensuring
consistent output regardless of which pass produced the hit.

### Memory mode vs Database mode

The Python API uses **memory mode**: mapping tables are downloaded from
upstream APIs (UniProt, BioMart, etc.), stored as `dict[str, set[str]]`
in memory, and cached as pickle files on disk. This is fast for
interactive use and small to medium batch jobs.

The REST API uses **database mode**: mapping tables are pre-loaded into
PostgreSQL during a build step. Translation queries run as SQL queries
against indexed tables. This handles concurrent requests and very large
ID sets efficiently.

Both modes share the same fallback logic and cleanup pipeline.
The difference is the data access layer:

| Aspect | Memory mode | Database mode |
|--------|-------------|---------------|
| Data loading | On-demand from upstream APIs | Pre-built into PostgreSQL |
| Lookup | Python dict (O(1)) | SQL query (indexed) |
| Fallbacks | `map_name()` per-ID pipeline | Batch SQL queries per fallback step |
| Cleanup | Same pipeline | Same pipeline |
| Caching | Pickle files, 5-min memory expiry | Database is the cache |
| Concurrency | Single-process | Multi-process via connection pool |

## REST API fallback chain

The REST API implements the same fallback strategies as the Python API,
but uses batch SQL queries instead of per-ID lookups.

### How REST fallbacks work

1. **Direct DB query:** All input IDs are looked up in a single SQL
   query against the `id_mapping` table.

2. **Gene symbol fallbacks (batch):** For unmapped IDs where the source
   type is `genesymbol`:
    - Uppercase batch: all unmapped IDs are uppercased and queried
    - Capitalize batch: remaining unmapped IDs are capitalized and queried
    - Synonym batch: remaining unmapped IDs are queried against
      `genesymbol-syn`
    - Uppercase synonym batch: remaining unmapped IDs uppercased and
      queried against `genesymbol-syn`

3. **Chain translation (batch):** For unmapped IDs where neither side
   is `uniprot`: a batch intermediate query (`source -> uniprot`),
   then a batch final query (`uniprot -> target`).

4. **UniProt cleanup:** Same pipeline as the Python API (secondary ->
   primary, TrEMBL -> SwissProt, proteome filter, format validation).

When `raw=true` is passed, steps 2-4 are skipped entirely.

### REST query parameters reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `identifiers` | string | required | Comma-separated IDs (GET) or JSON list (POST) |
| `id_type` | string | required | Source ID type |
| `target_id_type` | string | required | Target ID type |
| `ncbi_tax_id` | int | 9606 | Organism NCBI Taxonomy ID |
| `raw` | bool | false | Skip all special-case handling |
| `backend` | string | null | Force specific backend (metadata only in DB mode) |

**GET examples:**

```bash
# Basic translation
curl "https://omnipathdb.org/mapping/translate?\
identifiers=TP53,EGFR&\
id_type=genesymbol&\
target_id_type=uniprot"

# Raw mode -- direct lookup only
curl "https://omnipathdb.org/mapping/translate?\
identifiers=TP53&\
id_type=genesymbol&\
target_id_type=uniprot&\
raw=true"

# Specify organism (mouse)
curl "https://omnipathdb.org/mapping/translate?\
identifiers=Trp53&\
id_type=genesymbol&\
target_id_type=uniprot&\
ncbi_tax_id=10090"
```

**POST example:**

```bash
curl -X POST "https://omnipathdb.org/mapping/translate" \
     -H "Content-Type: application/json" \
     -d '{
       "identifiers": ["TP53", "EGFR", "BRCA1"],
       "id_type": "genesymbol",
       "target_id_type": "uniprot",
       "raw": false,
       "backend": null
     }'
```

## Backend details

### Auto-selection algorithm

When no backend is specified, the mapper's `_find_backends()` method
builds a candidate list:

1. **Column-based backends** (`uniprot`, `uniprot_ftp`, `biomart`):
   For each backend, it checks whether both the source and target ID
   types have a column defined in `id_types.yaml` under that backend's
   key. If both columns exist, the backend is added to the candidate list.

2. **Custom backends** (`mirbase`, `unichem`, `ramp`, `hmdb`, `metanetx`, `bigg`): These
   are always appended to the candidate list. They perform their own
   internal support check in their `read()` method and return an empty
   dict if the ID type pair is not supported.

The mapper tries each candidate in order. The first backend that
returns non-empty data wins.

### Backend-specific behavior

#### uniprot

- **Data:** UniProt AC, gene symbols, Entrez, HGNC, RefSeq, PDB, and
  all cross-references stored in UniProt entries.
- **Access:** `pypath.inputs.uniprot` (preferred) or direct HTTP to the
  UniProt REST API.
- **Organism-specific:** Yes. Queries are filtered by NCBI taxonomy ID.
- **Caching:** Pickle file per `(source, target, organism)` triple.
- **When to force:** Default for most protein-related translations.

#### uniprot_ftp

- **Data:** Same scope as `uniprot`, but downloaded from UniProt FTP
  idmapping files (bulk, pre-computed).
- **Access:** `pypath.inputs.uniprot_ftp` or direct HTTP to UniProt FTP.
- **Organism-specific:** Yes, but only 12 model organisms have
  pre-computed files.
- **Caching:** Same as `uniprot`.
- **When to force:** Faster for full-proteome builds; useful when the
  REST API is slow or rate-limited.

#### uploadlists

- **Data:** Same scope as UniProt REST.
- **Access:** Direct HTTP to the UniProt ID Mapping batch service
  (submit job, poll, collect results).
- **Organism-specific:** Yes.
- **When to force:** When you have a specific, bounded set of IDs and
  want the most up-to-date cross-references.

#### biomart

- **Data:** Ensembl gene/transcript/protein IDs, gene symbols, Entrez.
- **Access:** `pypath.inputs.biomart` or direct HTTP to the Ensembl
  BioMart XML service.
- **Organism-specific:** Yes.
- **Caching:** Pickle file per triple.
- **When to force:** When translating Ensembl IDs or when you need
  fresh Ensembl-specific data.

#### mirbase

- **Data:** Precursor names (`mir-pre`), mature names (`mir-mat-name`),
  miRBase accessions (`mirbase`).
- **Access:** `pypath.inputs.mirbase` (requires pypath).
- **Organism-specific:** Yes.
- **When to force:** miRNA translations.

#### unichem

- **Data:** Cross-references between chemical databases (ChEMBL, ChEBI,
  DrugBank, PubChem, KEGG, etc.).
- **Access:** `pypath.inputs.unichem` (requires pypath).
- **Organism-specific:** No (chemicals are universal).
- **When to force:** Chemical compound ID translation.

#### ramp

- **Data:** Metabolite cross-references and synonym mappings from
  RaMP-DB.
- **Access:** `pypath.inputs.ramp` (requires pypath).
- **Organism-specific:** No.
- **When to force:** Metabolite ID translation, especially when you
  need common-name-to-database-ID resolution.

#### hmdb

- **Data:** HMDB, PubChem, ChEBI, DrugBank, KEGG compound.
- **Access:** `pypath.inputs.hmdb` (requires pypath).
- **Organism-specific:** Human-derived data, but identifiers are universal.
- **When to force:** HMDB-centric metabolite lookups.

#### metanetx

- **Data:** Pairwise metabolite ID cross-references from MNXref
  `chem_xref.tsv` (3.4M cross-reference entries).
- **Access:** `pypath.inputs.metanetx.metanetx_mapping()` (requires pypath).
- **Organism-specific:** No (chemicals are universal).
- **Supported pairs:** bigg↔chebi, kegg↔chebi, hmdb↔chebi,
  lipidmaps↔chebi, swisslipids↔chebi, and all metanetx↔*
  combinations.
- **Coverage:** 82K hmdb→chebi, 45K kegg→chebi, 23K lipidmaps→chebi,
  11K bigg→chebi via MetaNetX bridge.
- **When to force:** Metabolite ID translation, especially for
  lipidmaps→chebi and swisslipids→chebi which are not available in
  other backends.

#### bigg

- **Data:** BiGG Models universal metabolite TSV (9,090 universal
  metabolites across 85+ models).
- **Access:** `pypath.inputs.bigg.bigg_metabolite_mapping()` (requires pypath).
- **Organism-specific:** No (chemicals are universal).
- **Supported pairs:** bigg↔chebi, bigg↔hmdb, bigg↔kegg,
  bigg↔metanetx.
- **Coverage:** 2,145 BiGG metabolites with ChEBI (10,319 pairs
  including ChEBI ontology hierarchy).
- **When to force:** BiGG metabolite ID translation. Combined with
  the MetaNetX backend, gives maximum BiGG→ChEBI coverage.

### Writing a new backend

To add a new data source:

1. Create a module in `omnipath_utils/mapping/backends/` (e.g.
   `_mybackend.py`).

2. Subclass `MappingBackend`:

    ```python
    from omnipath_utils.mapping.backends._base import MappingBackend
    from omnipath_utils.mapping.backends import register


    class MyBackend(MappingBackend):
        name = "mybackend"
        yaml_key = "mybackend"  # key in id_types.yaml

        def _read_via_pypath(
            self, id_type, target_id_type, ncbi_tax_id,
            *, src_col, tgt_col, **kwargs,
        ):
            import pypath.inputs.mybackend as mymod
            # ... fetch and return dict[str, set[str]]
            raise ImportError  # if pypath not available

        def _read_direct(
            self, id_type, target_id_type, ncbi_tax_id,
            *, src_col, tgt_col, **kwargs,
        ):
            # Direct HTTP implementation
            # Return dict[str, set[str]]
            return {}

    register("mybackend", MyBackend)
    ```

3. Add column definitions to `id_types.yaml` under the `mybackend` key
   for each supported ID type.

4. The backend will be automatically discovered by `_find_backends()`
   if it is a column-based backend, or by the custom backends list if
   you add it to `_CUSTOM_BACKENDS` in `Mapper`.

## HMDB identifier normalisation

HMDB identifiers have two historical formats: the old 5-digit format
(`HMDB00001`) and the current 7-digit format (`HMDB0000001`). All
translation APIs (Python and REST) automatically normalise the old format
to the 7-digit form. This normalisation is applied transparently at input
time -- both formats are accepted, and results always use the 7-digit
format.

The normalisation is applied by the mapper before any backend lookup, so
it works consistently across all backends that handle HMDB identifiers
(`hmdb`, `metanetx`, `bigg`, `ramp`, `unichem`).

## Database tables for special cases

During a database build, several auxiliary tables are created beyond
the main `id_mapping` table:

| Table content | Source -> Target | Purpose |
|---------------|-----------------|---------|
| Gene symbol -> SwissProt | `genesymbol -> swissprot` | TrEMBL-to-SwissProt cleanup step |
| TrEMBL -> gene symbol | `trembl -> genesymbol` | Reverse lookup in TrEMBL-to-SwissProt cleanup |
| Gene symbol synonyms -> UniProt | `genesymbol-syn -> uniprot` | Synonym fallback |
| SwissProt reference list | reflist | Proteome membership check (SwissProt) |
| TrEMBL reference list | reflist | Proteome membership check (TrEMBL) |
| HGNC IDs | `hgnc -> uniprot` | Additional protein ID coverage |
| RefSeq protein IDs | `refseqp -> uniprot` | Additional protein ID coverage |
| Secondary -> primary UniProt | `uniprot-sec -> uniprot-pri` | Obsolete AC resolution |

These tables enable the cleanup pipeline and fallback strategies to run
entirely within the database, without loading data into memory.

## Caching and performance

### Memory mode caching

Mapping tables are cached as pickle files in
`~/.cache/omnipath_utils/mapping/`. Each unique combination of
`(id_type, target_id_type, ncbi_tax_id, backend)` produces a
deterministic filename via MD5 hash:

```
mapping_genesymbol__uniprot__9606__a1b2c3d4e5f6.pickle
```

In-memory tables auto-expire after 5 minutes (300 seconds) of
inactivity. The `lifetime` parameter on the `Mapper` constructor
controls this. Expired tables are removed on the next
`remove_expired()` call.

To clear the cache:

```bash
rm -rf ~/.cache/omnipath_utils/mapping/
```

To change the cache directory:

```python
from omnipath_utils.mapping._mapper import Mapper
mapper = Mapper(cachedir='/custom/cache/path')
```

### Database mode performance

The database schema uses indexes on `(source_type_id, target_type_id,
ncbi_tax_id, source_id)` for fast lookups. The `COPY` command is used
during builds for fast bulk inserts. Partitioning by organism is
available for very large deployments.

### Translation performance tips

**Use batch functions for multiple IDs.** `translate()` and
`translate_column()` use `translate_core()` which performs a vectorized
first pass. Calling `map_name()` in a loop forces every ID through the
full fallback pipeline.

```python
# Good: vectorized first pass, fallback only for misses
result = translate(gene_list, 'genesymbol', 'uniprot')

# Bad: every ID goes through full pipeline
result = {g: map_name(g, 'genesymbol', 'uniprot') for g in gene_list}
```

**Use `raw=True` when you do not need fallbacks.** If your identifiers
are already in canonical form (e.g. uppercase gene symbols, primary
UniProt ACs), raw mode avoids all overhead:

```python
result = translate(clean_gene_list, 'genesymbol', 'uniprot', raw=True)
```

**Pre-fetch the translation table for repeated lookups.** If you are
translating IDs across multiple DataFrames or in a loop, fetch the
table once and reuse it:

```python
from omnipath_utils.mapping import translation_table

table = translation_table('genesymbol', 'uniprot')
# Use table directly as a dict
for gene in genes:
    uniprots = table.get(gene, set())
```

**For DataFrames, use `translate_column()`.** It handles deduplication,
batch lookup, and optional row expansion in a single call:

```python
from omnipath_utils.mapping import translate_column

df = translate_column(df, 'gene', 'genesymbol', 'uniprot')
```

## Identify and All-Mappings endpoints

### GET /mapping/identify

Given one or more identifiers, search all mapping tables to find which
ID types contain them. Useful when the type of an identifier is unknown.

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `identifiers` | string | yes | Comma-separated identifiers |
| `ncbi_tax_id` | int | no | NCBI Taxonomy ID (default: 9606) |

```bash
curl "https://omnipathdb.org/mapping/identify?\
identifiers=P04637,HMDB0000001"
```

**Response:**

```json
{
    "results": {
        "P04637": [
            {"id_type": "uniprot", "role": "source", "mappings_count": 5},
            {"id_type": "uniprot", "role": "target", "mappings_count": 1}
        ],
        "HMDB0000001": [
            {"id_type": "hmdb", "role": "source", "mappings_count": 3}
        ]
    },
    "meta": {
        "ncbi_tax_id": 9606,
        "total_input": 2
    }
}
```

Each match includes:

- **id_type** -- the canonical ID type name where the identifier was found.
- **role** -- `"source"` or `"target"`, indicating whether the identifier
  appears as a source or target in mapping rows.
- **mappings_count** -- the number of distinct mapped partners.

### GET /mapping/all

Given identifiers and their type, return all known mappings to every
other target type in a single request.

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `identifiers` | string | yes | Comma-separated identifiers |
| `id_type` | string | yes | Source ID type |
| `ncbi_tax_id` | int | no | NCBI Taxonomy ID (default: 9606) |

```bash
curl "https://omnipathdb.org/mapping/all?\
identifiers=P04637&\
id_type=uniprot"
```

**Response:**

```json
{
    "results": {
        "P04637": {
            "genesymbol": ["TP53"],
            "entrez": ["7157"],
            "ensg": ["ENSG00000141510"],
            "hgnc": ["11998"]
        }
    },
    "meta": {
        "id_type": "uniprot",
        "ncbi_tax_id": 9606,
        "total_input": 1
    }
}
```

### Python API

Both functions are available from `omnipath_utils.mapping`:

```python
from omnipath_utils.mapping import identify, all_mappings

# Identify unknown identifiers
identify(["P04637", "HMDB0000001"])

# Get all mappings
all_mappings(["P04637"], "uniprot")
```

And from the client:

```python
from omnipath_client.utils import identify, all_mappings

identify(["P04637", "HMDB0000001"])
all_mappings(["P04637"], "uniprot")
```

Both require database mode (PostgreSQL) on the server side.
