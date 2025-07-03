PyPath inputs Module — Quick Reference

(for fixing & maintaining pypath/inputs/)

⸻

1. What it is
	•	150 + data sources, ~500 download/parse functions.
	•	Each sub-module = one resource: pypath.inputs.{resource}
	•	Function naming: {resource}_{suffix}

Suffix	Meaning
_interactions	PPI networks
_enz_sub	Enzyme-substrate
_complexes	Complex membership
_annotations	Feature/GO etc.
_raw	Near-raw dump
_mapping	ID x-refs


⸻

2. Standard data flow
	1.	Download via pypath.share.curl.Curl
	•	Cache: ~/.pypath/cache/ (hash-named)
	•	Retries: 3 · Time-outs configurable
	2.	Parse file → Python structures
	3.	Map IDs to UniProt / NCBI TaxID (pypath.utils.mapping)
	4.	Assemble objects (interactions, complexes, annotations)

⸻

3. Frequent breakpoints & fixes

Symptom	Checklist	Quick fix
Download error	Net, SSL cert, corrupt cache	

with curl.cache_delete_on(): data = func()
with curl.debug_on(): data = func()
``` |
| New file layout | Compare fresh vs cached, adjust columns / delimiters | Update parser, re-test |
| Wrong IDs | Old mapping, species mix-up |  
```python
mapping.map_name('GENE_ID','src','uniprot',ncbi_tax_id=9606)
``` |

---

#### 4. Test toolkit  

```bash
# Run tests
python input_module_maintenance/test_input_modules.py --module signor
python input_module_maintenance/test_input_modules.py --batch-test interactions

Report fields: SUCCESS / record-count / sample / time.
Typical failures: network, format drift, mapping.

Manual smoke-test:

from pypath.inputs import signor
ints = signor.signor_interactions()
print(len(ints), ints[0])


⸻

5. Maintenance workflow (9-step cheat-sheet)
	1.	Locate failing function (logs / status file)
	2.	Open URL in browser → verify accessibility & format
	3.	Clear cache, re-download, diff vs old
	4.	Patch parser & error handling
	5.	Re-run automated + manual tests
	6.	Commit docstring comment on change
	7.	Update module_test_status.md (✓ count | ✗ reason)
	8.	Confirm downstream DB build passes
	9.	Push & PR

⸻

6. Key repo paths
	•	pypath/inputs/{resource}.py – individual modules
	•	pypath/share/settings.py – globals
	•	pypath/resources/urls.py – source URLs
	•	pypath/resources/data_formats.py – format specs
	•	input_module_maintenance/module_test_status.md – live health log
	•	input_module_maintenance/modules.json – module inventory