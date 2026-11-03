# Reference Lists

Reference lists are complete sets of identifiers for an organism.
They are useful for filtering, validation, and distinguishing reviewed
from unreviewed proteins.

## Available lists

| Name | Content |
|------|---------|
| `swissprot` | All reviewed UniProt accessions |
| `trembl` | All unreviewed UniProt accessions |
| `uniprot` | All UniProt accessions (SwissProt + TrEMBL) |

## Usage

### Check membership

```python
from omnipath_utils.reflists import is_swissprot, is_trembl

is_swissprot('P04637')  # True  (TP53 is reviewed)
is_trembl('P04637')     # False
```

### Get full lists

```python
from omnipath_utils.reflists import all_swissprots, all_trembls, all_uniprots

swissprots = all_swissprots()          # ~20,000 human reviewed ACs
trembls = all_trembls()                # ~200,000 human unreviewed ACs
all_acs = all_uniprots()               # union of both
```

### Other organisms

All functions accept `ncbi_tax_id`:

```python
from omnipath_utils.reflists import all_swissprots

mouse_sp = all_swissprots(ncbi_tax_id=10090)
```

### Generic access

```python
from omnipath_utils.reflists import get_reflist

ids = get_reflist('swissprot', ncbi_tax_id=9606)
```

## Data source

Reference lists are loaded from the UniProt REST API. When pypath is
available, it uses the pypath UniProt query interface; otherwise, it
falls back to direct HTTP requests to `rest.uniprot.org`.

Results are cached in memory after first load.
