---
name: inputs-v2-convert-legacy
description: Convert an existing legacy pypath input module from pypath.inputs into a modern inputs_v2 module. Use for legacy audit, upstream re-checking, semantic redesign, implementation, and parity-oriented smoke testing.
compatibility: Requires the pypath repository with both legacy inputs and inputs_v2 code available.
---

# inputs-v2-convert-legacy

Convert a legacy `pypath.inputs.*` module to `pypath.inputs_v2.*` in a deliberate way.

## Rules

- Do not blindly port old code.
- Treat the legacy module as a reference, not the specification.
- Re-check current upstream data before choosing the source.
- Preserve useful behavior intentionally; document deliberate differences.
- Validate both old and new behavior with small smoke tests.

## Workflow

### 1. Audit the legacy module
Read the old module fully.

Identify:
- public functions
- upstream source(s)
- return type(s)
- identifier mapping, filtering, flattening, and lossy behavior
- helper modules it depends on

If practical, run a small sample from the old module.

Also check resource metadata in:
- `pypath/resources/urls.py`
- `pypath/resources/descriptions.py`
- `pypath/resources/data/resources.json`
- optionally `.agent/omnipath-resources.csv` if useful

**Checkpoint:** state whether v1 should aim for parity, partial parity, or a cleaner redesign.

### 2. Re-check current upstream data
Inspect the current upstream source, not just the legacy URLs.

Decide:
- whether to keep the old source format
- or switch to a better current source
- what the current license, release pattern, and identifier richness are

**Checkpoint:** state the chosen source and why.

### 3. Redesign for `inputs_v2`
Read:
- `pypath/inputs_v2/README.md`
- `pypath/inputs_v2/base.py`
- `pypath/internals/tabular_builder.py`
- `pypath/internals/silver_schema.py`
- one or more similar modern modules such as `reactome.py`, `signor.py`, `intact.py`, `hmdb.py`

Decide:
- dataset names
- top-level entity type per dataset
- identifiers, annotations, membership
- what legacy behavior stays
- what intentionally changes
- whether new CV terms are needed

**Checkpoint:** state the dataset plan, semantic differences, and CV changes before coding.

### 4. Implement
Create or update:
- CV files if needed
- `pypath/inputs_v2/parsers/<resource>.py`
- `pypath/inputs_v2/<resource>.py`

Keep raw records flat and mapping-friendly. Prefer a focused v1 over a full legacy clone.

### 5. Validate and compare
Run, when possible:

```bash
python -m py_compile pypath/inputs_v2/<resource>.py pypath/inputs_v2/parsers/<resource>.py
```

```bash
poetry run python -c "from pypath.inputs import <legacy_module>; x = <legacy_call>; print(type(x)); print(next(iter(x)) if hasattr(x, '__iter__') else x)"
```

```bash
poetry run python -c "from pypath.inputs_v2.<resource> import resource; print(next(resource.<dataset>.raw()))"
```

```bash
poetry run python -c "from pypath.inputs_v2.<resource> import resource; print(next(resource.<dataset>()))"
```

Compare:
- coverage
- identifier richness
- semantics preserved or improved
- intentional differences

If exact parity is not possible or not desirable, say so clearly.

## What to report at the end

1. Legacy audit
2. Current upstream research
3. Semantic redesign
4. Files changed
5. Validation commands run
6. Parity / intentional differences
7. Next steps

Always include the exact `poetry run python -c ...` commands for the new `inputs_v2` module.
