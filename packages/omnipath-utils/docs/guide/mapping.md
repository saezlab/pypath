# ID Translation

!!! tip "Using omnipath-client"

    All the translation functions below are also available through the
    [omnipath-client](https://saezlab.github.io/omnipath-client/) package,
    which queries the web service and requires no local setup:

    ```python
    from omnipath_client.utils import map_name, translate_column

    map_name("TP53", "genesymbol", "uniprot")  # {"P04637"}

    # Translate a DataFrame column (pandas, polars, or pyarrow)
    translate_column(df, "gene", "genesymbol", "uniprot")
    ```

    ```bash
    pip install omnipath-client
    ```

## Overview

omnipath-utils translates between 97 biological identifier types --
gene symbols, UniProt accessions, Ensembl IDs, Entrez gene IDs, small
molecule identifiers, miRNA names, and more. Data comes from UniProt,
Ensembl BioMart, miRBase, HMDB, RaMP, UniChem, MetaNetX, and BiGG.

Biological ID mapping is inherently one-to-many. A gene symbol can
correspond to multiple UniProt accessions (reviewed and unreviewed entries,
isoforms of different genes sharing a symbol). A UniProt accession can map
to multiple Ensembl gene IDs when gene models differ between databases.
omnipath-utils returns `set[str]` to reflect this reality.

The mapper also handles the messy details: outdated secondary UniProt
accessions, versioned Ensembl and RefSeq identifiers, case-variant gene
symbols, CURIE prefixes, and confused miRNA precursor/mature forms. When
no direct mapping table exists, it chains through UniProt as an
intermediate (e.g. Entrez &rarr; UniProt &rarr; Ensembl).

## Python API

### Core functions

All functions are available from `omnipath_utils.mapping`:

```python
from omnipath_utils.mapping import (
    map_name,
    map_names,
    map_name0,
    translate,
    translation_table,
    id_types,
)
```

#### `map_name` -- translate a single identifier

```python
def map_name(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    raw: bool = False,
    backend: str | None = None,
) -> set[str]
```

Returns all target identifiers matching the input. Empty set if no
mapping is found.

```python
map_name('TP53', 'genesymbol', 'uniprot')
# {'P04637'}

map_name('P04637', 'uniprot', 'genesymbol')
# {'TP53'}

map_name('TP53', 'genesymbol', 'ensg')
# {'ENSG00000141510'}

map_name('HMDB0000122', 'hmdb', 'chebi')
# {'15903'}
```

#### `map_names` -- translate multiple, return union

```python
def map_names(
    names: Iterable[str],
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    raw: bool = False,
    backend: str | None = None,
) -> set[str]
```

Translates each input identifier individually and returns the union of all
results. Useful when you need a flat set of targets and do not need to
know which input produced which output.

```python
map_names(['TP53', 'EGFR', 'BRCA1'], 'genesymbol', 'uniprot')
# {'P04637', 'P00533', 'P38398'}
```

#### `map_name0` -- translate to a single result

```python
def map_name0(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    raw: bool = False,
    backend: str | None = None,
) -> str | None
```

Convenience wrapper that picks one result from the set. Returns `None` if
no mapping exists. If the mapping is ambiguous (multiple targets), the
choice is arbitrary.

```python
map_name0('TP53', 'genesymbol', 'uniprot')
# 'P04637'

map_name0('NONEXISTENT', 'genesymbol', 'uniprot')
# None
```

#### `translate` -- batch translate with per-input results

```python
def translate(
    identifiers: Iterable[str],
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    raw: bool = False,
    backend: str | None = None,
) -> dict[str, set[str]]
```

Returns a dict mapping each input to its set of targets. Inputs that
could not be translated map to an empty set.

```python
translate(['TP53', 'EGFR', 'FAKE'], 'genesymbol', 'uniprot')
# {'TP53': {'P04637'}, 'EGFR': {'P00533'}, 'FAKE': set()}
```

!!! note
    `translate` uses vectorized table lookup for the first pass, then
    falls back to per-ID `map_name` (with full special-case handling)
    for any identifiers that miss in the table. Use `raw=True` to
    restrict to table lookup only with no fallbacks.

#### Full UniProt across all organisms -- `full_uniprot`

The curated tables cover the common organisms. In **database-backed mode** a
comprehensive table (`id_mapping_ftp`) built from the complete UniProt idmapping —
**all organisms, ~1 billion rows** — is also available. The `full_uniprot` argument
(on `map_name`, `map_names`, `translate`, and the `/mapping/translate` endpoint)
controls how it is used:

| value | behaviour |
|-------|-----------|
| `'fallback'` (default) | curated table first; the full table is consulted **only for identifiers the curated table did not resolve** |
| `'never'` | curated only |
| `'both'` | query both and merge (deduplicated) |
| `'only'` | the full table only |

```python
# EGFR in chimpanzee (9598): absent from the curated tables, resolved from the full table
map_name('EGFR', 'genesymbol', 'uniprot', ncbi_tax_id=9598)                     # resolves
map_name('EGFR', 'genesymbol', 'uniprot', ncbi_tax_id=9598, full_uniprot='never')  # empty
```

The `id_mapping_ftp` row count is reported by the `/health` endpoint. The default
(`fallback`) keeps the common case fast on the small curated table and only touches
the billion-row table when needed.

#### `translation_table` -- full mapping table

```python
def translation_table(
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
) -> dict[str, set[str]]
```

Returns the entire mapping table as a dict. This is the raw table --
every source identifier known to the backend, mapped to all its targets.

```python
table = translation_table('genesymbol', 'uniprot')
table['TP53']
# {'P04637'}
len(table)
# ~20000 for human
```

#### `id_types` -- list all supported types

```python
def id_types() -> list[str]
```

Returns canonical names of all 97 supported ID types.

```python
id_types()
# ['uniprot', 'swissprot', 'trembl', 'genesymbol', 'genesymbol-syn',
#  'entrez', 'ensg', 'ensp', 'enst', 'refseqp', 'hgnc', 'hmdb',
#  'chebi', 'pubchem', 'drugbank', 'mirbase', 'mir-pre', ...]
```

### Parameters

All translation functions accept these parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `name` / `names` / `identifiers` | `str` / `Iterable[str]` | required | The identifier(s) to translate |
| `id_type` | `str` | required | Source ID type (e.g. `genesymbol`, `uniprot`, `ensg`, `hmdb`) |
| `target_id_type` | `str` | required | Target ID type |
| `ncbi_tax_id` | `int \| None` | `9606` (human) | NCBI Taxonomy ID for the organism |
| `raw` | `bool` | `False` | Skip all special-case handling (direct table lookup only) |
| `backend` | `str \| None` | `None` | Force a specific backend (e.g. `uniprot`, `biomart`) |

For details on these parameters, see [Advanced Translation](mapping-advanced.md).

The `Mapper.map_name` method (accessed via the singleton) accepts two
additional parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `strict` | `bool` | `False` | Skip fuzzy gene symbol fallbacks (case variants, synonym lookup, "1" suffix) |
| `uniprot_cleanup_flag` | `bool` | `True` | When target is `uniprot`, run the cleanup pipeline (secondary &rarr; primary, SwissProt preference, proteome filter) |

These are not exposed in the module-level convenience functions, but
can be accessed directly:

```python
from omnipath_utils.mapping._mapper import Mapper

Mapper.get().map_name(
    'TP53', 'genesymbol', 'uniprot',
    strict=True,
    uniprot_cleanup_flag=False,
)
```

### One-to-many results

Translation results are always sets because biological ID mapping is
inherently one-to-many:

- A gene symbol may map to multiple UniProt accessions. For example,
  `HBB` maps to the main hemoglobin beta chain (P68871) plus
  potentially unreviewed TrEMBL entries.
- A single Ensembl gene may correspond to multiple UniProt entries
  if the gene has been split or merged across databases.
- Small molecule databases assign different identifiers to the same
  compound, or the same identifier to stereoisomers.

```python
map_name('HBB', 'genesymbol', 'uniprot')
# Could return {'P68871'} or {'P68871', 'A0A0S2Z4L3', ...}
# depending on cleanup settings and available data
```

When you need exactly one result, use `map_name0` -- but be aware that
the choice among multiple candidates is arbitrary.

## REST API

The web service exposes translation via HTTP endpoints. These use the
database backend (PostgreSQL) rather than in-memory tables.

### GET /mapping/translate

Translate a comma-separated list of identifiers.

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `identifiers` | string | yes | Comma-separated identifiers |
| `id_type` | string | yes | Source ID type |
| `target_id_type` | string | yes | Target ID type |
| `ncbi_tax_id` | int | no | NCBI Taxonomy ID (default: 9606) |
| `raw` | bool | no | Skip special-case handling (default: false) |
| `backend` | string | no | Force specific backend (default: null) |

```bash
curl "https://omnipathdb.org/mapping/translate?\
identifiers=TP53,EGFR,BRCA1&\
id_type=genesymbol&\
target_id_type=uniprot"
```

### POST /mapping/translate

For large ID lists (hundreds or thousands of identifiers), use the POST
endpoint with a JSON body.

**JSON body:**

```json
{
    "identifiers": ["TP53", "EGFR", "BRCA1", "..."],
    "id_type": "genesymbol",
    "target_id_type": "uniprot",
    "ncbi_tax_id": 9606,
    "raw": false,
    "backend": null
}
```

```bash
curl -X POST "https://omnipathdb.org/mapping/translate" \
     -H "Content-Type: application/json" \
     -d '{
       "identifiers": ["TP53", "EGFR", "BRCA1"],
       "id_type": "genesymbol",
       "target_id_type": "uniprot",
       "ncbi_tax_id": 9606
     }'
```

### GET /mapping/id-types

Returns all supported ID types with metadata.

```bash
curl "https://omnipathdb.org/mapping/id-types"
```

### Response format

Both GET and POST `/mapping/translate` return the same JSON structure:

```json
{
    "results": {
        "TP53": ["P04637"],
        "EGFR": ["P00533"],
        "BRCA1": ["P38398"]
    },
    "unmapped": ["NONEXISTENT"],
    "meta": {
        "id_type": "genesymbol",
        "target_id_type": "uniprot",
        "ncbi_tax_id": 9606,
        "total_input": 4,
        "total_mapped": 3,
        "raw": false,
        "backend": null
    }
}
```

- **results** -- dict mapping each successfully translated input to a
  sorted list of target identifiers.
- **unmapped** -- list of input identifiers that could not be translated.
- **meta** -- request parameters and summary counts.

The `/mapping/id-types` endpoint returns a list of objects:

```json
[
    {
        "name": "uniprot",
        "label": "UniProt AC",
        "entity_type": "protein",
        "curie_prefix": "uniprot"
    },
    ...
]
```

## UniProt behavior

This is the most important section of this document. The UniProt cleanup
pipeline runs automatically whenever the target ID type is `uniprot`, and
it substantially affects results.

### SwissProt vs TrEMBL

UniProt has two sections: **SwissProt** (reviewed, manually curated) and
**TrEMBL** (unreviewed, computationally predicted). For human, SwissProt
contains ~20,400 entries while TrEMBL adds ~200,000 more. Most
bioinformatics workflows want SwissProt entries.

By default, when `target_id_type='uniprot'`, the cleanup pipeline runs
after every successful translation step. The pipeline has four stages:

**Step 1: Secondary &rarr; primary AC translation.** Some resources store
obsolete secondary UniProt accessions. The cleanup maps these to their
current primary AC using the `uniprot-sec` &rarr; `uniprot-pri` table.
If no secondary mapping exists, the AC is assumed to already be primary.

**Step 2: TrEMBL &rarr; SwissProt preference.** For each result AC, the
pipeline checks whether it is in the SwissProt reference list. If it is,
the AC is kept. If it is a TrEMBL entry, the pipeline looks up its gene
symbol (via the `trembl` &rarr; `genesymbol` table), then finds the
SwissProt entry for that symbol. If a SwissProt entry exists, it replaces
the TrEMBL AC. If no SwissProt entry is found for that gene, the TrEMBL
AC is kept.

**Step 3: Organism proteome filter.** The result set is intersected with
the organism's full proteome (all UniProt ACs for that NCBI taxonomy ID).
This removes stale or misassigned ACs. If the filter would discard all
results (e.g. due to an incomplete proteome list), the unfiltered set is
returned.

**Step 4: Format validation.** Each AC is checked against the UniProt AC
regex pattern (`^[OPQ][0-9][A-Z0-9]{3}[0-9]$` or the extended 10-character
format). Invalid strings are discarded.

### Controlling UniProt behavior

Five patterns cover common use cases:

```python
# Default: full cleanup, prefers SwissProt
map_name('TP53', 'genesymbol', 'uniprot')
# {'P04637'}  -- P04637 is the SwissProt entry

# Explicitly request only SwissProt (reviewed) entries
map_name('TP53', 'genesymbol', 'swissprot')
# {'P04637'}

# Explicitly request only TrEMBL (unreviewed) entries
map_name('TP53', 'genesymbol', 'trembl')
# Unreviewed entries only; may return empty set if gene
# has no TrEMBL entries

# Disable cleanup: get raw results from the mapping table
Mapper.get().map_name('TP53', 'genesymbol', 'uniprot',
                       uniprot_cleanup_flag=False)
# May include TrEMBL entries, secondary ACs, entries from
# other organisms

# Same type: cleanup still runs
map_name('Q9Y4K3', 'uniprot', 'uniprot')
# Translates secondary -> primary if Q9Y4K3 is a secondary AC
```

The three target ID types and their behavior:

| Target type | Backend filter | Cleanup pipeline | Result |
|------------|----------------|-----------------|--------|
| `uniprot` | Both SwissProt + TrEMBL | Yes (secondary &rarr; primary, TrEMBL &rarr; SwissProt, proteome filter, format check) | Prefers SwissProt, keeps TrEMBL only when no SwissProt exists |
| `swissprot` | SwissProt only (`reviewed=True`) | No | Only reviewed entries |
| `trembl` | TrEMBL only (`reviewed=False`) | No | Only unreviewed entries |

### UniProt &rarr; gene symbol

When mapping a UniProt AC to `genesymbol`, the system first checks the
SwissProt gene name table. If the AC is not found there, it tries the
TrEMBL table. If neither has it, the secondary &rarr; primary chain is
attempted: the AC is looked up in `uniprot-sec` &rarr; `uniprot-pri`,
and the resulting primary AC is looked up again.

## Translation pipeline

When you call `map_name('TP53', 'genesymbol', 'uniprot')`, the mapper
runs through an ordered sequence of strategies until one produces results.
Here is the full pipeline:

### 1. Alias resolution

ID type names are normalized via `IdTypeRegistry.resolve()`. Aliases and
variant spellings are mapped to canonical names:

- `genesymbol_syn` &rarr; `genesymbol-syn`
- `GeneSymbol` &rarr; `genesymbol`
- `gene_symbol` &rarr; `genesymbol`
- `ensembl_gene_id` &rarr; `ensg`

### 2. Same-type shortcut

If `id_type == target_id_type`, the input is returned as-is. Exception:
if the target is `uniprot` and cleanup is enabled, the cleanup pipeline
still runs (to resolve secondary ACs and filter the proteome).

### 3. Direct table lookup

The mapper looks for a loaded or loadable mapping table for the exact
`(source, target, organism)` triple. If the table exists and contains
the input, the result is returned.

### 4. Gene symbol fallbacks

Only triggered when `id_type` is `genesymbol` or `genesymbol-syn`. The
system tries progressively looser matches:

**(a) UPPER case.** Tries `name.upper()`. Human gene symbols are
uppercase (`TP53`), but input may be mixed case (`Tp53`).

**(b) Capitalized.** Tries `name.capitalize()` (first letter upper, rest
lower). Rodent gene symbols follow this convention (`Trp53` for mouse).

**(c) Synonym table.** Looks up the name in the `genesymbol-syn` table.
Gene symbols change over time; the synonym table maps old names to current
ones. Both exact and uppercase variants are tried.

**(d) Append "1".** Tries `name + "1"`. Some gene families have members
where the "1" suffix is optional in common usage (e.g. `ACTA` vs `ACTA1`).
**Skipped in strict mode.**

### 5. RefSeq version handling

Only triggered when `id_type` starts with `refseq`. RefSeq accessions
include a version suffix (e.g. `NM_000546.6`). If the exact ID is not
found:

- Strips the version suffix and tries the base ID (`NM_000546`)
- In non-strict mode, iterates common version numbers 1--19
  (`NM_000546.1`, `NM_000546.2`, ...) until a match is found

### 6. Ensembl version stripping

Only triggered when `id_type` starts with `ens` and the input contains a
dot. Strips the version suffix:

`ENSG00000141510.18` &rarr; `ENSG00000141510`

### 7. miRNA reciprocal fallback

Only triggered when `id_type` starts with `mir-`. Data sources often
confuse mature and precursor miRNA forms. If a direct lookup for
`mir-mat-name` (mature name, e.g. `hsa-miR-21-5p`) fails, the system
tries it as `mir-name` (precursor name), maps to a miRBase accession as
intermediate, then maps to the target. The reverse direction works the
same way.

### 8. CURIE prefix stripping

Only triggered when the input contains `:`. Strips the prefix and retries:

`CHEBI:15903` &rarr; `15903`

### 9. Chain translation

Only triggered when neither `id_type` nor `target_id_type` is `uniprot`.
The system chains through UniProt as an intermediate:

`entrez` &rarr; `uniprot` &rarr; `ensg`

Each leg of the chain runs through the full `map_name` pipeline
(including all fallback strategies).

### 10. Reverse lookup

If no forward table exists, the mapper checks for a reverse table
(`target` &rarr; `source`) and scans all values to find entries containing
the input. This is a linear scan and slower than a direct lookup, but it
avoids the need to maintain separate reverse tables.

### 11. UniProt cleanup

Applied after **any** successful step if the target is `uniprot` and
`uniprot_cleanup_flag` is True. See the [UniProt behavior](#uniprot-behavior)
section for the full cleanup pipeline.

### Strict mode

When `strict=True`, the following fallbacks are skipped:

- Gene symbol step 4d (append "1")
- RefSeq version iteration (steps beyond stripping the version suffix)

Strict mode is useful when you need exact matches and want to avoid false
positives from fuzzy matching.

## Backends

### How backends are selected

Backend selection is automatic. For each `(source, target)` pair, the
mapper checks which backends have column definitions for both types in
`id_types.yaml`. The column-based backends (`uniprot`, `uniprot_ftp`,
`biomart`) are checked first. Custom backends (`mirbase`, `unichem`,
`ramp`, `hmdb`) are always appended to the candidate list; they perform
their own support check internally.

The first backend that successfully returns data wins. If a backend
fails (network error, missing data), the next one is tried.

### Available backends

| Backend | Data source | Coverage | Organism-specific | Data access |
|---------|-------------|----------|-------------------|-------------|
| `uniprot` | UniProt REST API | UniProt AC, gene symbols, Entrez, HGNC, RefSeq, PDB, and all cross-references in UniProt | Yes | pypath.inputs.uniprot &rarr; direct HTTP |
| `uniprot_ftp` | UniProt FTP idmapping files | Same as `uniprot`, but bulk download per organism | Yes (12 model organisms) | pypath.inputs.uniprot_ftp &rarr; direct HTTP |
| `uploadlists` | UniProt ID Mapping batch service | Same scope as UniProt, but for targeted ID sets | Yes | Direct HTTP (submit/poll/collect) |
| `biomart` | Ensembl BioMart | Ensembl gene/transcript/protein IDs, gene symbols, Entrez | Yes | pypath.inputs.biomart &rarr; direct HTTP |
| `mirbase` | miRBase | Precursor names, mature names, miRBase accessions | Yes | pypath.inputs.mirbase |
| `unichem` | UniChem (EMBL-EBI) | Cross-references between chemical databases (ChEMBL, ChEBI, DrugBank, PubChem, etc.) | No | pypath.inputs.unichem |
| `ramp` | RaMP-DB | Metabolite cross-references plus synonym mappings | No | pypath.inputs.ramp |
| `hmdb` | HMDB | HMDB, PubChem, ChEBI, DrugBank, KEGG compound | Human only | pypath.inputs.hmdb |
| `metanetx` | MNXref `chem_xref.tsv` (3.4M cross-reference entries) | Pairwise metabolite ID translation via MetaNetX bridge: bigg↔chebi, kegg↔chebi, hmdb↔chebi, lipidmaps↔chebi, swisslipids↔chebi, and all metanetx↔* combinations | No | pypath.inputs.metanetx |
| `bigg` | BiGG Models universal metabolite TSV (9,090 universal metabolites across 85+ models) | bigg↔chebi, bigg↔hmdb, bigg↔kegg, bigg↔metanetx | No | pypath.inputs.bigg |

### Pypath integration

Most backends try to use `pypath.inputs` first. If pypath is not
installed, the `uniprot` and `biomart` backends fall back to direct HTTP
requests against the upstream APIs. The `mirbase`, `unichem`, `ramp`, and
`hmdb` backends require pypath (they raise `ImportError` if it is
unavailable). The `uploadlists` backend always uses direct HTTP.

### Using a specific backend (developer info)

Backends can be called directly, bypassing the mapper's automatic
selection. This is useful for debugging or when you need raw data
from a specific source.

```python
from omnipath_utils.mapping.backends import get_backend

# Load a UniChem mapping table
b = get_backend('unichem')
data = b.read('chembl', 'chebi', 0)
# data: {'CHEMBL25': {'15365'}, 'CHEMBL612': {'17303'}, ...}

# Load an Ensembl BioMart table
b = get_backend('biomart')
data = b.read('ensg', 'genesymbol', 9606)
# data: {'ENSG00000141510': {'TP53'}, ...}
```

The `read()` method returns `dict[str, set[str]]`. The third argument
is `ncbi_tax_id`; pass `0` for organism-independent backends.

## Small molecule identifiers

Small molecule mappings are provided by five backends:

- **UniChem** -- cross-references between chemical databases maintained by
  EMBL-EBI. Covers ChEMBL, ChEBI, DrugBank, PubChem, KEGG, and others.
- **RaMP** -- the RaMP-DB multi-source metabolite harmonization database.
  Provides both primary ID cross-references and synonym mappings (common
  names to database IDs).
- **HMDB** -- the Human Metabolome Database. Maps between HMDB, PubChem,
  ChEBI, DrugBank, and KEGG compound identifiers.
- **MetaNetX** -- the MNXref namespace reconciliation database. Provides
  pairwise metabolite ID translation via MetaNetX as a bridge identifier.
  Covers 82K hmdb→chebi, 45K kegg→chebi, 23K lipidmaps→chebi, and
  11K bigg→chebi mappings. Supported pairs include bigg↔chebi,
  kegg↔chebi, hmdb↔chebi, lipidmaps↔chebi, swisslipids↔chebi, and
  all metanetx↔* combinations.
- **BiGG** -- the BiGG Models database of genome-scale metabolic network
  reconstructions. Provides bigg↔chebi, bigg↔hmdb, bigg↔kegg, and
  bigg↔metanetx mappings from 9,090 universal metabolites across 85+
  models. Coverage includes 2,145 BiGG metabolites with ChEBI (10,319
  pairs including ChEBI ontology hierarchy). Combined with MetaNetX,
  gives maximum BiGG→ChEBI coverage.

Small molecule identifiers are not organism-specific. Backends receive
`ncbi_tax_id=0` (or ignore the parameter). HMDB data is human-derived
but the identifiers themselves are universal.

```python
map_name('HMDB0000122', 'hmdb', 'chebi')
# {'15903'}

map_name('CHEMBL25', 'chembl', 'drugbank')
# {'DB00945'}  -- aspirin

map_name('15903', 'chebi', 'pubchem')
# {'5793'}

# ChEBI to HMDB
map_name('15422', 'chebi', 'hmdb')

# PubChem to ChEBI
map_name('5957', 'pubchem', 'chebi')

# HMDB to KEGG
map_name('HMDB0000001', 'hmdb', 'kegg')
```

### HMDB identifier normalisation

HMDB identifiers have two historical formats: the old 5-digit format
(`HMDB00001`) and the current 7-digit format (`HMDB0000001`). The mapper
automatically normalises the old format to 7-digit in all translation
APIs (Python and REST). This is applied transparently -- you can pass
either format as input, and results always use the 7-digit format.

```python
# Both formats work; results always use 7-digit
map_name('HMDB00001', 'hmdb', 'chebi')
# {'16044'}

map_name('HMDB0000001', 'hmdb', 'chebi')
# {'16044'}
```

## Identifying unknown identifiers

When you have identifiers but do not know their type, use the `identify`
function to search all mapping tables:

```python
from omnipath_utils.mapping import identify

identify(['P04637', 'HMDB0000001'])
# {'P04637': [{'id_type': 'uniprot', 'role': 'source', 'mappings_count': 5}, ...],
#  'HMDB0000001': [{'id_type': 'hmdb', 'role': 'source', 'mappings_count': 3}, ...]}
```

Each result entry includes:

- **id_type** -- the ID type where the identifier was found.
- **role** -- whether the identifier appears as a `source` or `target`
  in mapping tables.
- **mappings_count** -- how many distinct mappings exist from/to that
  identifier.

This requires database mode (PostgreSQL).

### REST API

```bash
curl "https://omnipathdb.org/mapping/identify?\
identifiers=P04637,HMDB0000001"
```

## Get all mappings for an identifier

To retrieve all known mappings for an identifier to every other type,
use `all_mappings`:

```python
from omnipath_utils.mapping import all_mappings

all_mappings(['P04637'], 'uniprot')
# {'P04637': {'genesymbol': ['TP53'], 'entrez': ['7157'], ...}}
```

This returns a nested dict: `{identifier: {target_type: [target_ids]}}`.

This requires database mode (PostgreSQL).

### REST API

```bash
curl "https://omnipathdb.org/mapping/all?\
identifiers=P04637&\
id_type=uniprot"
```

## miRNA identifiers

miRNA translation uses the miRBase backend. Three ID types are supported:

| ID type | Description | Example |
|---------|-------------|---------|
| `mir-pre` | Precursor miRNA name | `hsa-mir-21` |
| `mir-mat-name` | Mature miRNA name | `hsa-miR-21-5p` |
| `mirbase` | miRBase accession | `MI0000077` (precursor), `MIMAT0000076` (mature) |

Data sources frequently confuse precursor and mature forms. The
reciprocal fallback (pipeline step 7) handles this: if you look up a
mature name but the table stores it as a precursor (or vice versa), the
mapper swaps the assumed type, chains through the miRBase accession, and
reaches the target.

```python
map_name('hsa-miR-21-5p', 'mir-mat-name', 'mirbase')
# {'MIMAT0000076'}

map_name('MI0000077', 'mirbase', 'mir-mat-name')
# {'hsa-miR-21-5p', 'hsa-miR-21-3p'}
```

miRNA mappings are organism-specific. Pass `ncbi_tax_id` for non-human
organisms:

```python
map_name('mmu-miR-21a-5p', 'mir-mat-name', 'mirbase', ncbi_tax_id=10090)
```

## Caching

Mapping tables are cached as pickle files in
`~/.cache/omnipath_utils/mapping/`. Each unique combination of `(id_type,
target_id_type, ncbi_tax_id, backend)` produces a deterministic cache
filename based on an MD5 hash.

In-memory tables auto-expire after 5 minutes (300 seconds) of inactivity.
The `lifetime` parameter on the `Mapper` constructor controls this.

The cache directory is configurable:

```python
from omnipath_utils.mapping._mapper import Mapper

mapper = Mapper(cachedir='/path/to/cache')
```

To clear the cache, delete the pickle files:

```bash
rm -rf ~/.cache/omnipath_utils/mapping/
```

## Database mode

When PostgreSQL is available (deployment scenario), the REST API queries
the database directly via SQL rather than loading mapping tables into
memory. The database stores pre-computed mapping tables in a normalized
schema (`id_mapping` table with `source_type_id`, `target_type_id`,
`source_id`, `target_id`, `ncbi_tax_id`).

The Python API uses the in-memory mode by default. The `Mapper` singleton
manages table loading, caching, and expiry. For the deployed web service,
translation queries hit PostgreSQL through SQLAlchemy, bypassing the
in-memory machinery entirely.

See the [Database Build](../deployment/database.md) and
[Web Service](../deployment/server.md) pages for deployment details.
