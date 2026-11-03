# Taxonomy

omnipath-utils resolves organism identifiers across naming systems used
by different databases.

## `ensure_*` functions

Each function accepts any supported organism representation and returns
the requested name form:

```python
from omnipath_utils.taxonomy import (
    ensure_ncbi_tax_id,
    ensure_common_name,
    ensure_latin_name,
    ensure_ensembl_name,
    ensure_kegg_code,
    ensure_mirbase_name,
    ensure_oma_code,
)
```

### Input forms

All `ensure_*` functions accept any of these as input:

| Form | Example |
|------|---------|
| NCBI Taxonomy ID (int) | `9606` |
| NCBI Taxonomy ID (str) | `'9606'` |
| Common name | `'human'` |
| Latin name | `'Homo sapiens'` |
| Short latin | `'hsapiens'` |
| Ensembl name | `'hsapiens'` |
| KEGG code | `'hsa'` |
| miRBase code | `'hsa'` |
| OMA code | `'HUMAN'` |
| UniProt code | `'HUMAN'` |

### Examples

```python
ensure_ncbi_tax_id('human')     # 9606
ensure_ncbi_tax_id('hsapiens')  # 9606
ensure_ncbi_tax_id('hsa')       # 9606
ensure_ncbi_tax_id(9606)         # 9606

ensure_common_name(10090)        # 'mouse'
ensure_latin_name(9606)          # 'Homo sapiens'
ensure_ensembl_name(9606)        # 'hsapiens'
ensure_kegg_code(9606)           # 'hsa'
```

Returns `None` if the organism is not recognized.

## List all organisms

```python
from omnipath_utils.taxonomy import all_organisms

orgs = all_organisms()
# dict: {9606: {'common_name': 'human', ...}, 10090: {...}, ...}
```

## Supported organisms

omnipath-utils ships with 22 organisms pre-configured in `organisms.yaml`,
including all major model organisms:

| NCBI Tax ID | Common name | KEGG | Ensembl |
|------------|-------------|------|---------|
| 9606 | human | hsa | hsapiens |
| 10090 | mouse | mmu | mmusculus |
| 10116 | rat | rno | rnorvegicus |
| 9913 | cow | bta | btaurus |
| 9031 | chicken | gga | ggallus |
| 7955 | zebrafish | dre | drerio |
| 7227 | fruitfly | dme | dmelanogaster |
| 6239 | nematode | cel | celegans |
| 4932 | yeast | sce | scerevisiae |
| 3702 | thale cress | ath | athaliana |

And more. Use `all_organisms()` for the complete list.
