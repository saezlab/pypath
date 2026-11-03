# FoodOn Mapping Plan (No Cross-Resource Dependency)

## Goal
Create a FoodOn mapping artifact via a standalone script, and have `foodb` and `phenol_explorer` inputs_v2 modules **optionally** enrich their records only if the mapping file is present. This avoids runtime dependencies between resources.

## Design
- **Mapping artifact is committed to git** in `pypath/pypath/inputs_v2/enrichment/`.
- **Inputs modules do not download FoodOn**; they only read the mapping file if it exists.
- **Regeneration is manual/occasional** via a script.

## Proposed File Locations
- Script: `scripts/build_foodon_mappings.py` (or `pypath/pypath/inputs_v2/scripts/build_foodon_mappings.py`)
- Mapping file: `pypath/pypath/inputs_v2/enrichment/foodon_mappings.csv`

## Mapping File Schema (proposed)
Single file with a `source` column:
- `source` ("foodb" | "phenol_explorer")
- `source_id` (resource ID; e.g., FooDB public_id or Phenol-Explorer id)
- `source_name`
- `foodon_id`
- `foodon_name`
- `method` ("ncbi_taxonomy" | "exact_name" | "scientific_name")
- `ncbi_taxon` (if applicable)

## Script Behavior
1. Load FoodB foods + Phenol-Explorer foods using existing inputs_v2 parsers.
2. Download and parse FoodOn terms (NCBI taxon, label, synonyms).
3. Build mappings with priority:
   - NCBI taxonomy ID
   - Exact name
   - Scientific name
4. Write combined CSV mapping file into `pypath/pypath/inputs_v2/enrichment/`.
5. Include a small header or README entry with build date/version (optional).

## Inputs_v2 Enrichment Behavior
- `foodb.py` and `phenol_explorer.py`:
  - Try to load mapping CSV from `inputs_v2/enrichment/`.
  - If present, join on `source_id` (and/or `source_name`) to add:
    - `IdentifierNamespaceCv.FOODON` (identifier)
    - Optional annotation fields like `foodon_name` and `mapping_method`.
  - If missing, skip these annotations silently.

## Open Questions
- Confirm exact `source_id` field for matching (FooDB: `public_id` vs internal `id`).
- Confirm where to store method/name annotations (annotation CV term vs generic).
- Confirm script location preference.
