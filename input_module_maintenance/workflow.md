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
# Test single module (RECOMMENDED - process one at a time)
python input_module_maintenance/test_input_modules.py --module signor

# Test specific function
python input_module_maintenance/test_input_modules.py --function signor.signor_interactions

# List available modules
python input_module_maintenance/test_input_modules.py --list-modules

Report fields: SUCCESS / record-count / sample / time.
Typical failures: network, format drift, mapping.

Manual smoke-test for debugging:

from pypath.inputs import signor
ints = signor.signor_interactions()
print(len(ints), ints[0])


⸻

5. Maintenance workflow (parallel testing approach)

### 5a. Systematic Testing Phase (RECOMMENDED)
	1.	Use 10 parallel subagents to test modules efficiently:
		• Launch 10 Task agents simultaneously with module names
		• Each agent tests one module and reports back results
		• Significantly faster than one-by-one testing
		• Update module_test_status.md with all results at once
	2.	Document patterns immediately:
		• ✅ Working modules (no action needed)
		• ⚠️ Parameter-dependent functions (expected behavior)
		• ❌ Failed modules (need investigation/fixes)
		• ⏳ Timeout modules (large datasets or network issues)
	3.	Commit batch documentation after each round of 10

### 5b. Fix Phase (after systematic testing)
	1.	Prioritize fixes based on patterns:
		• Infrastructure issues (missing methods, imports)
		• Simple typos (like _cons → _const)
		• URL/API changes
		• Format parsing errors
	2.	Fix ONE issue at a time:
		• Test the specific function manually first
		• Apply minimal fix
		• Re-test to confirm
		• Update module_test_status.md: ✅ FIXED
		• Commit with clear message
	3.	Skip non-critical issues:
		• 403 Forbidden (external site blocks)
		• Functions requiring parameters (document as ⚠️)
		• Large dataset timeouts (document as ⏳)

### 5c. Legacy Individual Testing (when needed)
	1.	Test ONE module at a time using test script
	2.	If test fails, investigate the specific error
	3.	For download/format errors:
		• Check URL in browser
		• Clear cache with curl.cache_delete_on()
		• Compare column names/format changes
	4.	Fix the module code and re-test

⸻

6. Key repo paths
	•	pypath/inputs/{resource}.py – individual modules
	•	pypath/share/settings.py – globals
	•	pypath/resources/urls.py – source URLs
	•	pypath/resources/data_formats.py – format specs
	•	input_module_maintenance/module_test_status.md – live health log
	•	input_module_maintenance/modules.json – module inventory