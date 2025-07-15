#!/usr/bin/env python3
"""
Fix tags in YAML files to match PyPath's definitions exactly.
"""

import yaml
import glob
from typing import Dict, List, Set

# PyPath's definitions for module categories
INTERACTION_MODULES = {
    'acsn', 'baccin2019', 'biogrid', 'cancerdrugsdb', 'ccmap', 'cellcall', 
    'cellchatdb-cofactors', 'cellinker', 'cellphonedb', 'cpdb', 'dip', 
    'embrace', 'hi2', 'hi_union', 'hippie', 'hpmr', 'hprd', 'huri', 
    'innatedb', 'intact', 'italk', 'kirouac2010', 'lit13', 'lrdb', 
    'matrixdb', 'mppi', 'nci_pid', 'ramilowski2015', 'reactome', 
    'signor', 'wojtowicz2020', 'yang2016', 'yu2011'
}

COMPLEX_MODULES = {
    'cellphonedb', 'cellchatdb', 'cellinker', 'compleat', 'complexportal',
    'corum', 'guidetopharmacology', 'havugimana', 'humap', 'humap2',
    'icellnet', 'kegg', 'pdb', 'signor', 'spike'
}

ANNOTATION_MODULES = {
    'Adhesome', 'Almen2009', 'Annot', 'Biomart', 'Cancerrxgene', 
    'Cancersea', 'Celltypist', 'Cellphonedb', 'Clinvar', 'ComPPI', 
    'Compartments', 'Cpad', 'CellChatDB', 'DGIdb', 'DisGeNet', 
    'DiseaseOntology', 'DoRothEA', 'Drugbank', 'E-GEOD-5258', 
    'ExoCarta', 'GOInterpro', 'Goslim', 'HGNC', 'HIPPIE', 
    'HPA', 'HumanCellMap', 'IntOGen', 'Italk', 'Kirouac2010', 
    'Lambert2018', 'LRdb', 'Matrisome', 'Matrixdb', 'Membranome', 
    'Msigdb', 'NetPath', 'OmniPath', 'OPM', 'PROGENy', 
    'PanglaoDB', 'Phosphatome', 'PreFoldedPPI', 'ProteinAtlas', 
    'SIGNOR', 'Signaling_Pathway_Impact_Analysis_SPIA', 'Slk2003', 
    'Surfaceome', 'TOPDB', 'TFcensus', 'Tfcheckpoint', 'UniProt', 
    'Vesiclepedia', 'Wang', 'Zhong2015'
}

def normalize_module_name(name: str) -> str:
    """Normalize module names for comparison."""
    # Convert to lowercase for comparison
    return name.lower().replace('_', '').replace('-', '')

def should_have_interaction_tag(module_name: str) -> bool:
    """Check if module should have interaction tags."""
    normalized = normalize_module_name(module_name)
    return any(normalize_module_name(m) == normalized for m in INTERACTION_MODULES)

def should_have_complex_tag(module_name: str) -> bool:
    """Check if module should have protein_complex in entity tags."""
    normalized = normalize_module_name(module_name)
    return any(normalize_module_name(m) == normalized for m in COMPLEX_MODULES)

def should_have_annotation_tag(module_name: str) -> bool:
    """Check if module should have annotation tags."""
    # For annotations, also check exact match with different case
    for ann_module in ANNOTATION_MODULES:
        if module_name.lower() == ann_module.lower():
            return True
    return False

def fix_module_tags(module: Dict) -> Dict:
    """Fix tags for a single module based on PyPath definitions."""
    module_name = module.get('module_name', '')
    
    # Initialize tags if not present
    if 'tags' not in module:
        module['tags'] = {}
    
    # Fix interaction tags
    if 'interaction' in module['tags']:
        if not should_have_interaction_tag(module_name):
            module['tags']['interaction'] = []
    
    # Fix entity tags (specifically protein_complex)
    if 'entity' in module['tags'] and isinstance(module['tags']['entity'], list):
        if 'protein_complex' in module['tags']['entity']:
            if not should_have_complex_tag(module_name):
                # Remove protein_complex from entity list
                module['tags']['entity'] = [e for e in module['tags']['entity'] if e != 'protein_complex']
    
    # Fix annotation tags
    if 'annotation' in module['tags']:
        if not should_have_annotation_tag(module_name):
            module['tags']['annotation'] = []
    
    return module

def process_yaml_file(filepath: str) -> None:
    """Process a single YAML file to fix tags."""
    print(f"Processing {filepath}...")
    
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)
    
    if 'modules' not in data:
        print(f"  Warning: No modules found in {filepath}")
        return
    
    # Track changes
    changes = []
    
    for module in data['modules']:
        module_name = module.get('module_name', 'unnamed')
        original_tags = module.get('tags', {}).copy()
        
        # Fix the module
        module = fix_module_tags(module)
        
        # Check for changes
        new_tags = module.get('tags', {})
        if original_tags != new_tags:
            changes.append({
                'module': module_name,
                'original': original_tags,
                'fixed': new_tags
            })
    
    # Save the file
    with open(filepath, 'w') as f:
        yaml.dump(data, f, sort_keys=False, default_flow_style=False, width=1000)
    
    # Report changes
    if changes:
        print(f"  Fixed {len(changes)} modules:")
        for change in changes:
            print(f"    - {change['module']}:")
            if 'interaction' in change['original'] and change['original']['interaction'] != change['fixed'].get('interaction', []):
                print(f"      interaction: {change['original']['interaction']} -> {change['fixed'].get('interaction', [])}")
            if 'entity' in change['original'] and change['original']['entity'] != change['fixed'].get('entity', []):
                print(f"      entity: {change['original']['entity']} -> {change['fixed'].get('entity', [])}")
            if 'annotation' in change['original'] and change['original']['annotation'] != change['fixed'].get('annotation', []):
                print(f"      annotation: {change['original']['annotation']} -> {change['fixed'].get('annotation', [])}")
    else:
        print(f"  No changes needed")

def main():
    """Main function to process all YAML files."""
    yaml_files = sorted(glob.glob('result-*.yaml'))
    
    print(f"Found {len(yaml_files)} YAML files to process\n")
    
    for filepath in yaml_files:
        process_yaml_file(filepath)
    
    print("\nDone!")

if __name__ == "__main__":
    main()