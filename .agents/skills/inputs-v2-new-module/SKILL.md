---
name: inputs-v2-new-module
description: Create a new pypath inputs_v2 module for an external database. Use for structured research, semantic mapping, CV review, implementation, and poetry-based smoke testing.
compatibility: Requires the pypath repository with Poetry and the inputs_v2 framework.
---

# inputs-v2-new-module

Create a new `pypath.inputs_v2` module in a small number of explicit steps.

## Rules

- Research first; do not start coding immediately.
- Prefer the upstream format with the best semantics and identifiers.
- Design for the silver schema; avoid unsupported deep nesting.
- Add new CV terms only when existing ones are not semantically correct.
- Validate with `poetry run python -c ...` before finishing.

## Workflow

### 1. Research the source
Read:
- `pypath/inputs_v2/README.md`
- related modules in `pypath/inputs_v2/`
- relevant entries in `pypath/resources/urls.py`, `pypath/resources/descriptions.py`, and `pypath/resources/data/resources.json`

Then inspect current upstream data.

Decide:
- what formats exist
- which format is best
- what entity types / identifiers / species / license are available

**Checkpoint:** state the chosen source format and why.

### 2. Design the semantic mapping
Read:
- `pypath/inputs_v2/base.py`
- `pypath/internals/tabular_builder.py`
- `pypath/internals/silver_schema.py`
- one similar module, e.g. `pypath/inputs_v2/signor.py`

Decide:
- dataset names
- top-level entity type per dataset
- identifiers, annotations, and membership
- whether new CV terms are needed

Keep the plan simple and explicit.

**Checkpoint:** state the dataset plan and CV changes before coding.

### 3. Design raw parser output
If base parsers are enough, use them. Otherwise create `pypath/inputs_v2/parsers/<resource>.py`.

Raw records should be:
- flat
- stable
- easy to map with `FieldConfig` / `EntityBuilder`

Prefer helpers for:
- URL resolution
- archive iteration
- xref extraction
- entity type inference

**Checkpoint:** inspect and print at least one raw record.

### 4. Implement
Create or update:
- CV files if needed
- `pypath/inputs_v2/parsers/<resource>.py`
- `pypath/inputs_v2/<resource>.py`

Use:
- `ResourceConfig`
- `Download`
- `Dataset`
- `EntityBuilder`

Start with the highest-value dataset(s). Do not overbuild v1.

### 5. Validate
Run, when possible:

```bash
python -m py_compile pypath/inputs_v2/<resource>.py pypath/inputs_v2/parsers/<resource>.py
```

```bash
poetry run python -c "from pypath.inputs_v2.<resource> import resource; print(next(resource.<dataset>.raw()))"
```

```bash
poetry run python -c "from pypath.inputs_v2.<resource> import resource; print(next(resource.<dataset>()))"
```

Repeat for each dataset.

## What to report at the end

1. Research summary
2. Semantic design
3. Files changed
4. Validation commands run
5. Open questions / next steps

Always include the exact `poetry run python -c ...` commands for the new module.
