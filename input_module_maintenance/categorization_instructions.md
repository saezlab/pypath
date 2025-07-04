# PyPath Input Module Categorization Instructions

## Overview
Analyze and categorize pypath input modules according to the criteria below. Each module should be examined by reading its source code in `/pypath/inputs/` and understanding its functionality.

## Output Format (YAML)
Each module should be categorized using the following YAML structure:

```yaml
module_name: string
functions:
  - list of function names from modules.json

maintenance_status:
  category: one_time_paper|discontinued|infrequent|frequent
  last_known_update: YYYY-MM or unknown
  notes: additional context about maintenance status

module_type:
  - primary_data    # Sources that provide biological data
  - utility         # Helper modules for processing/formatting
  - id_mapping      # Modules for identifier conversion

dataset_size: small|medium|large
# small: <1MB or <1000 records
# medium: 1MB-100MB or 1000-100k records  
# large: >100MB or >100k records

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
    # Add others as needed
```

## Categorization Guidelines

### Maintenance Status
- **one_time_paper**: Resources from paper supplements, unlikely to be updated
- **discontinued**: Databases that are no longer maintained
- **infrequent**: Updated yearly or less frequently
- **frequent**: Updated monthly or more frequently

Look for clues in:
- Comments mentioning years or versions
- URL patterns (e.g., archived sites, specific paper DOIs)
- Database names (well-known active vs defunct databases)

### Module Type
- **primary_data**: Modules that fetch/parse biological data
- **utility**: Helper modules (formatting, common functions)
- **id_mapping**: Modules specifically for ID conversion

### Dataset Size
Estimate based on:
- Number of interactions/annotations typically provided
- File size references in code
- Known database sizes

### Technical Aspects
- **access_method**: How data is obtained
- **data_format**: Actual file formats processed (don't use predefined list)

### Tags
Apply all relevant tags. A module can have multiple tags in each category.

## Analysis Process
1. Read the module source code
2. Identify the main functions listed in modules.json
3. Understand what data the module provides
4. Look for maintenance/version clues
5. Categorize according to the schema
6. Output clean YAML for each module

## Important Notes
- If a tag doesn't fit the predefined categories, add a new one
- Some modules may have submodules (e.g., hmdb/metabolites) - treat path as module name
- Be concise in the notes field
- When uncertain about maintenance status, look for:
  - Last update dates in URLs or comments
  - Whether the source website is still active
  - Version numbers or publication years