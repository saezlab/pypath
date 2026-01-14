---
description: Translate an old pypath.inputs module to the new inputs_v2 declarative pattern
---

# Translate Input Module to inputs_v2

This workflow converts a legacy `pypath.inputs` module to the new declarative `pypath.inputs_v2` pattern.

## Prerequisites

- Understanding of the module's data source (URL, format, fields)
- Access to the old module at `pypath/pypath/inputs/<module>.py`
- Familiarity with the inputs_v2 pattern (see `pypath/pypath/inputs_v2/README.md`)
- Virtual environment activated: `source /Users/jschaul/Code/omnipath_build/.venv/bin/activate`

---

## Step 1: Analyze the Old Module and Gather Metadata

### 1a. Check existing metadata sources

**IMPORTANT**: Resource metadata (license, URL, pubmed, description) exists in:
- `pypath/pypath/resources/data/resources.json` - Primary source for license, pubmeds, URLs
- `pypath/pypath/resources/descriptions.py` - Detailed descriptions
- `pypath/pypath/resources/urls.py` - Download URLs

// turbo
Search for resource metadata:
```bash
grep -A 50 '"RESOURCE_NAME":' pypath/pypath/resources/data/resources.json | head -60
```

### 1b. Open and analyze `pypath/pypath/inputs/<module>.py`:

1. **Identify data sources**: URLs, POST params, file types
2. **Identify output fields**: What namedtuple/dict fields does it return?
3. **Identify entity types**: Proteins, interactions, complexes, etc.
4. **Identify transformations**: Regex patterns, mappings, filtering

Create a checklist of:
- [ ] Download URL(s)
- [ ] File format (TSV, CSV, JSON, XML, etc.)
- [ ] Entity type(s) (protein, interaction, complex, etc.)
- [ ] Identifier fields (UniProt, SIGNOR ID, etc.)
- [ ] Annotation fields (descriptions, scores, etc.)
- [ ] Membership fields (complex members, interaction participants)

---

## Step 2: Check CV Terms

Ensure necessary CV terms exist in `pypath/pypath/internals/cv_terms/`:

1. **ResourceCv** - `resources.py`: Add resource ID if missing
2. **EntityTypeCv** - `entity_types.py`: Check entity types exist
3. **IdentifierNamespaceCv** - `identifiers.py`: Add namespaces if needed
4. **Annotation CVs** - `annotations.py`: Various annotation types

### CV Term Selection Guidelines

**When to REUSE an existing term:**
- The semantic meaning matches your use case (e.g., `PUBMED` for PubMed IDs)
- The data type is the same (accessions, descriptions, scores)
- You want this data to be queryable/groupable with similar data from other resources

**When to ADD a new term:**
- The data has a specific meaning not captured by existing terms (e.g., `FUNCAT` for FunCat annotations)
- You need to distinguish this data type from similar-looking data
- The source uses a specific ontology/classification scheme

**Examples:**
- ✅ Use `IdentifierNamespaceCv.PUBMED` for PubMed IDs from any source
- ✅ Add `MoleculeAnnotationsCv.FUNCAT` for FunCat-specific annotations
- ❌ Don't use `DESCRIPTION` for data with specific ontology meaning

// turbo
Run to check existing CV terms:
```bash
grep -r "RESOURCE_NAME" pypath/pypath/internals/cv_terms/
```

---

## Step 3: Create ResourceConfig

```python
from pypath.inputs_v2.base import ResourceConfig
from pypath.internals.cv_terms import ResourceCv, LicenseCV, UpdateCategoryCV

config = ResourceConfig(
    id=ResourceCv.<RESOURCE_ID>,
    name='<Human Readable Name>',
    url='<Resource Homepage URL>',
    license=LicenseCV.<LICENSE_TYPE>,  # Usually CC_BY_4_0
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='<PubMed ID>',
    description='<Description of the resource>',
)
```

---

## Step 4: Create FieldConfig

Define extraction patterns for field parsing:

```python
from pypath.internals.tabular_builder import FieldConfig

f = FieldConfig(
    extract={
        'pattern_name': r'<regex_pattern>',  # Named regex patterns
    },
    map={
        'mapping_name': {'value': CV_TERM},  # Named mappings
    },
    transform={
        'transform_name': lambda x: x.upper(),  # Named transforms
    },
    delimiter=',',  # Default delimiter for multi-value fields
)
```

---

## Step 5: Create EntityBuilder Schema

```python
from pypath.internals.tabular_builder import (
    EntityBuilder, IdentifiersBuilder, AnnotationsBuilder, CV
)
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv

schema = EntityBuilder(
    entity_type=EntityTypeCv.<TYPE>,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.<NAMESPACE>, value=f('<column_name>')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('<name_column>')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=<AnnotationCv>.<TYPE>, value=f('<annotation_column>')),
    ),
)
```

For entities with members (complexes, interactions):
```python
from pypath.internals.tabular_builder import (
    MembershipBuilder, Member, MembersFromList
)

schema = EntityBuilder(
    entity_type=EntityTypeCv.COMPLEX,
    identifiers=IdentifiersBuilder(...),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.UNIPROT,
                   value=f('members_column', delimiter=';')),
            ),
        )
    ),
)
```

---

## Step 6: Select or Create Raw Parser

### File Structure Pattern

The inputs_v2 modules follow a separation pattern:
- **Main module** (`inputs_v2/<resource>.py`): Schema, config, Resource definition
- **Parser** (`inputs_v2/parsers/<resource>.py`): Raw file parsing logic

### When to use base parsers (no separate file needed)

For simple tabular formats, use parsers from `pypath.inputs_v2.parsers.base`:
- `iter_tsv` - Tab-separated values
- `iter_csv` - Comma-separated values  
- `iter_json` - JSON files
- `iter_jsonl` - JSON Lines

Example using base parser directly:
```python
from pypath.inputs_v2.parsers.base import iter_tsv
# Use iter_tsv directly as raw_parser in Dataset
```

### When to create a custom parser file

Create `parsers/<resource>.py` when you need:
- Custom file format handling (XML, ZIP extraction, etc.)
- Organism/record filtering during parsing
- Field normalization before schema mapping
- Complex multi-file handling

Custom parser in `parsers/<resource>.py`:
```python
def _raw(opener, organism: int = 9606, **_kwargs):
    if not opener or not opener.result:
        return
    # Custom parsing logic
    for record in parse_data(opener.result):
        if should_include(record, organism):
            yield {'field1': record.x, 'field2': record.y}
```

Then import in main module:
```python
from pypath.inputs_v2.parsers.<resource> import _raw
```

---

## Step 7: Create Resource

```python
from pypath.inputs_v2.base import Dataset, Download, Resource

resource = Resource(
    config,
    <dataset_name>=Dataset(
        download=Download(
            url='<data_url>',
            filename='<local_filename>',
            subfolder='<resource_name>',
            large=True,
            download_kwargs={'post': True, 'query': {...}},  # Optional
        ),
        mapper=schema,
        raw_parser=iter_tsv,  # or custom parser
    ),
)
```

---

## Step 8: Test the Module

First, activate the virtual environment:
```bash
source /Users/jschaul/Code/omnipath_build/.venv/bin/activate
cd /Users/jschaul/Code/omnipath_build/pypath
```

// turbo
Test raw parsing:
```bash
python -c "
from pypath.inputs_v2.<module> import resource
raw = list(resource.<dataset>.raw())
print(f'Raw records: {len(raw)}')
if raw: print(f'Sample fields: {list(raw[0].keys())[:5]}')
"
```

// turbo
Test entity generation:
```bash
python -c "
from pypath.inputs_v2.<module> import resource
entity = next(resource.<dataset>())
print(f'Type: {entity.type}')
print(f'Identifiers: {entity.identifiers[:2]}')
print(f'Members: {len(entity.membership) if entity.membership else 0}')
"
```

Compare with old module output to validate equivalence.

---

## Checklist Template

For translating `<module_name>`:

- [ ] Analyzed old module structure
- [ ] Verified/added CV terms
- [ ] Created ResourceConfig
- [ ] Created FieldConfig (if needed)
- [ ] Created EntityBuilder schema(s)
- [ ] Selected/created raw parser
- [ ] Created Resource with Dataset(s)
- [ ] Tested raw parsing
- [ ] Tested entity generation
- [ ] Compared with old module output
