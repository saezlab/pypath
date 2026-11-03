# Quick Start

## Translate identifiers

```python
from omnipath_utils.mapping import map_name, map_names, translate

# Single gene symbol to UniProt
map_name('TP53', 'genesymbol', 'uniprot')
# {'P04637'}

# Multiple identifiers at once
map_names(['TP53', 'EGFR'], 'genesymbol', 'uniprot')
# {'P04637', 'P00533'}

# Batch translate with per-input results
translate(['TP53', 'EGFR'], 'genesymbol', 'uniprot')
# {'TP53': {'P04637'}, 'EGFR': {'P00533'}}

# Get a single result (or None)
from omnipath_utils.mapping import map_name0
map_name0('TP53', 'genesymbol', 'uniprot')
# 'P04637'
```

## Resolve organism names

```python
from omnipath_utils.taxonomy import (
    ensure_ncbi_tax_id,
    ensure_common_name,
    ensure_kegg_code,
)

# Any name form to NCBI Taxonomy ID
ensure_ncbi_tax_id('human')     # 9606
ensure_ncbi_tax_id('hsapiens')  # 9606
ensure_ncbi_tax_id('hsa')       # 9606

# NCBI Taxonomy ID to other forms
ensure_common_name(10090)       # 'mouse'
ensure_kegg_code(9606)          # 'hsa'
```

## Check reference lists

```python
from omnipath_utils.reflists import is_swissprot, all_swissprots

# Check if a UniProt AC is reviewed
is_swissprot('P04637')  # True (TP53 is reviewed)

# Get all reviewed human UniProt ACs
swissprots = all_swissprots()
len(swissprots)  # ~20,000
```
