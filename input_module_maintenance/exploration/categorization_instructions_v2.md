# PyPath Input Module Categorization Instructions (v2.0)
## Updated Based on PyPath Source Code Validation

## Overview
Analyze and categorize pypath input modules according to validated criteria based on PyPath's actual source code definitions. Each module should be examined by reading its source code in `/pypath/inputs/` and cross-referenced with PyPath's official resource lists.

## Ground Truth References
Cross-reference your categorization with these PyPath source files:
- `/pypath/resources/data_formats.py` - Official interaction resource definitions
- `/pypath/core/complex.py` - Official complex resource list (`complex_resources`)
- `/pypath/core/annot.py` - Official annotation resource defaults (`protein_sources_default`)

## Updated Module Classification Framework

### Core PyPath Resource Types
Based on validation against PyPath source code:

#### 1. **Interaction Resources** (33 official resources)
Modules that provide protein-protein, regulatory, or other molecular interactions:
```yaml
# Official PyPath interaction resources include:
# biogrid, intact, signor, hprd, dip, mppi, etc.
interaction_priority: core|extended|utility
```

#### 2. **Complex Resources** (15 official resources)  
Modules that provide protein complex data:
```yaml
# Official PyPath complex resources:
# Signor, Corum, CellPhoneDB, Havugimana, Compleat, ComplexPortal, 
# Pdb, GuideToPharmacology, Humap, Humap2, Icellnet, Kegg, 
# Cellchatdb, Cellinker, Spike
complex_priority: core|extended
```

#### 3. **Annotation Resources** (75 official resources)
Modules that provide protein annotations:
```yaml
# Official PyPath annotation resources include:
# Dgidb, Membranome, Exocarta, Matrisome, Surfaceome, etc.
annotation_priority: core|extended|utility
```

## Enhanced Output Format (YAML)

```yaml
module_name: string
functions:
  - list of function names from modules.json

# New validation fields
pypath_classification:
  official_resource: true|false  # Found in PyPath source definitions
  resource_type: interaction|complex|annotation|utility|mixed
  priority: core|extended|utility  # Based on PyPath usage

maintenance_status:
  category: one_time_paper|discontinued|infrequent|frequent
  last_known_update: YYYY-MM or unknown
  notes: additional context about maintenance status
  validation_notes: any discrepancies with expected status

module_type:
  - primary_data    # Sources that provide biological data
  - utility         # Helper modules for processing/formatting  
  - id_mapping      # Modules for identifier conversion

dataset_size: small|medium|large
species_coverage: human_only|multi_species

technical_aspects:
  access_method: 
    - api            # REST/SOAP API calls
    - file_download  # Direct file downloads
    - web_scraping   # HTML parsing/scraping
  data_format:
    - List actual formats used (e.g., tsv, csv, xml, json, rdf, obo, etc.)

tags:
  entity:
    - protein
    - protein_complex
    - rna
    - metabolite
    - drug
    - small_molecule
    # Add others as needed
  
  annotation:
    - structural               # 3D structures, domains
    - functional              # GO terms, functions
    - tissue_location         # Tissue expression
    - cellular_location       # Cell type annotations
    - subcellular_location    # Compartments, organelles
    - disease_association     # Disease links
    - pathway                 # Pathway membership
    - post_translational_modifications
    - isoforms
    - variants                # Mutations, SNPs
    # Add others as needed
  
  interaction:
    - protein_protein_undirected    # Physical interactions
    - protein_protein_regulation    # Directed regulatory
    - ligand_receptor
    - enzyme_substrate             # Including PTMs
    - protein_small_molecule       # Drug targets
    - transcriptional_regulation   # TF-gene
    - mirna                        # miRNA-target
    - enzyme_metabolite            # Metabolic reactions
    - drug_drug_interaction        # Drug-drug interactions
    - protein_rna                  # Protein-RNA interactions
    # Add others as needed

# New validation tracking
validation_info:
  matches_pypath_definition: true|false
  expected_but_missing: []  # Resources expected but not found
  unexpected_but_present: []  # Resources found but not in official lists
```

## Updated Categorization Guidelines

### PyPath Resource Priority Classification

#### **Core Resources** (High Priority)
- Listed in official PyPath resource definitions
- Essential for OmniPath database construction
- Well-maintained, frequently used
- Examples: BioGRID, IntAct, SIGNOR, CORUM, ComplexPortal

#### **Extended Resources** (Medium Priority)  
- May be used by PyPath but not in core lists
- Specialized datasets or alternative sources
- Examples: HIPPIE, HURI, pathway-specific resources

#### **Utility Resources** (Support Priority)
- Helper modules for data processing
- ID mapping and format conversion
- API access utilities
- Examples: uniprot_db, ebi, pubmed utilities

### Validation-Based Guidelines

#### 1. **Cross-Reference Validation**
Before finalizing categorization:
```bash
# Check if module appears in PyPath source:
grep -r "module_name" /pypath/resources/data_formats.py
grep -r "ModuleName" /pypath/core/complex.py  
grep -r "ModuleName" /pypath/core/annot.py
```

#### 2. **Resource Type Determination**
- **interaction**: Module provides molecular interactions (check data_formats.py)
- **complex**: Module provides protein complex data (check complex.py)
- **annotation**: Module provides protein annotations (check annot.py)
- **utility**: Helper functions, not primary data
- **mixed**: Modules serving multiple purposes

#### 3. **Priority Assessment**
- **core**: Official PyPath resource in source definitions
- **extended**: Used by PyPath but not in core lists
- **utility**: Support functions, format converters, APIs

## Quality Assurance Checklist

### Before Submitting Results:
1. ✅ **Cross-referenced with PyPath source code**
2. ✅ **Verified module exists in expected location**
3. ✅ **Checked for maintenance indicators in code**
4. ✅ **Validated interaction/complex/annotation classifications**
5. ✅ **Confirmed technical access methods**
6. ✅ **Documented any discrepancies or unexpected findings**

## Missing Resource Investigation

### Known Gaps to Investigate:
Based on validation, pay special attention to:

#### Missing Interaction Resources:
- `ccmap` (CancerCellMap), `hi2`, `hiunion`, `lit13`, `yang2016`, `yu2011`

#### Missing Complex Resources:  
- `humap2` (Human Protein Map v2)

#### Missing Annotation Resources:
- UniProt-related: `uniprotfamilies`, `uniprotkeywords`, `uniprotlocations`
- Protein Atlas variants: `humanproteinatlas*`, `cellsurfaceproteinatlas*`
- Pathway modules: `keggpathways*`, `signorpathways`, `netpathpathways`

## Analysis Process (Updated)

1. **Read the module source code**
2. **Cross-reference with PyPath official resource lists**
3. **Identify functions from modules.json**
4. **Determine PyPath classification and priority**
5. **Validate expected vs actual resource type**
6. **Look for maintenance/version clues**
7. **Document validation findings**
8. **Output enhanced YAML with validation info**

## Important Notes

- **Priority**: Focus on core PyPath resources first
- **Validation**: Always note discrepancies between expected and actual findings
- **Coverage**: Aim for 100% coverage of official PyPath resources
- **Documentation**: Record validation process and decision rationale
- **Consistency**: Use standardized naming conventions matching PyPath source