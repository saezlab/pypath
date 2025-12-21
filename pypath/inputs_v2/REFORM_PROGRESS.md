# inputs_v2 Refactor Progress

## Achieved
- Introduced `Dataset`/`Resource` interface and migrated inputs_v2 modules to dataset bricks.
- Removed legacy wrapper functions; module-level datasets are the API (e.g., `guidetopharma.targets`).
- Folded `Map` into `Column` and added `FieldConfig` to centralize per-file processing.
- Replaced `Map` usage across inputs_v2 modules with `Column(..., extract=..., map=...)`.
- Consolidated processing helpers into a single local block for:
  - `pypath/pypath/inputs_v2/guidetopharma.py`
  - `pypath/pypath/inputs_v2/intact.py` 
- Rolled out the “processing block + FieldConfig + f(...) everywhere” pattern across remaining inputs_v2 modules:
  - `pypath/pypath/inputs_v2/bindingdb.py`
  - `pypath/pypath/inputs_v2/hmdb.py`
  - `pypath/pypath/inputs_v2/lipidmaps.py`
  - `pypath/pypath/inputs_v2/signor.py`
  - `pypath/pypath/inputs_v2/swisslipids.py`
  - `pypath/pypath/inputs_v2/uniprot.py`
- Updated `pypath/pypath/inputs_v2/RESOURCE_CONFIG_GUIDE.md` to document `FieldConfig` usage and remove `Map` references.
- Switched remaining schemas to `f(...)` and removed direct `Column(...)` calls from inputs_v2 resource modules.
- Standardized on `FieldConfig` factories across inputs_v2 modules (including ontologies).

## Final Cleanup
- Done.
