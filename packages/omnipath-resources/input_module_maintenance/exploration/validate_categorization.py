#!/usr/bin/env python3
"""
Validate PyPath module categorization against actual PyPath source code definitions.
Cross-references our categorization results with the official PyPath resource lists.
"""

import yaml
import re
from collections import defaultdict
from pathlib import Path

def extract_interaction_resources():
    """Extract interaction resources from data_formats.py"""
    resources = set()
    
    # Read the data_formats.py file
    with open('/Users/jschaul/Code/pypath-1/pypath/resources/data_formats.py', 'r') as f:
        content = f.read()
    
    # Find interaction dictionary
    interaction_match = re.search(r'interaction\s*=\s*{(.*?)\n}', content, re.DOTALL)
    if interaction_match:
        interaction_content = interaction_match.group(1)
        # Extract resource names (keys in the dictionary)
        resource_matches = re.findall(r"'([^']+)':\s*input_formats", interaction_content)
        resources.update(resource_matches)
    
    # Also check interaction_misc, interaction_htp, etc.
    misc_patterns = [
        r'interaction_misc\s*=\s*{(.*?)\n}',
        r'interaction_htp\s*=\s*{(.*?)\n}',
        r'ligand_receptor\s*=\s*{(.*?)\n}',
        r'small_molecule_protein\s*=\s*{(.*?)\n}',
    ]
    
    for pattern in misc_patterns:
        match = re.search(pattern, content, re.DOTALL)
        if match:
            misc_content = match.group(1)
            resource_matches = re.findall(r"'([^']+)':\s*input_formats", misc_content)
            resources.update(resource_matches)
    
    return resources

def extract_complex_resources():
    """Extract complex resources from complex.py"""
    resources = set()
    
    with open('/Users/jschaul/Code/pypath-1/pypath/core/complex.py', 'r') as f:
        content = f.read()
    
    # Find complex_resources tuple
    match = re.search(r'complex_resources\s*=\s*\((.*?)\)', content, re.DOTALL)
    if match:
        resources_content = match.group(1)
        # Extract quoted resource names
        resource_matches = re.findall(r"'([^']+)'", resources_content)
        resources.update(resource_matches)
    
    return resources

def extract_annotation_resources():
    """Extract annotation resources from annot.py"""
    resources = set()
    
    with open('/Users/jschaul/Code/pypath-1/pypath/core/annot.py', 'r') as f:
        content = f.read()
    
    # Find protein_sources_default set
    match = re.search(r'protein_sources_default\s*=\s*{(.*?)\}', content, re.DOTALL)
    if match:
        sources_content = match.group(1)
        # Extract quoted resource names
        resource_matches = re.findall(r"'([^']+)'", sources_content)
        resources.update(resource_matches)
    
    return resources

def load_categorization_results():
    """Load all our categorization results"""
    all_modules = []
    
    for i in range(1, 11):
        try:
            with open(f'/Users/jschaul/Code/pypath-1/input_module_maintenance/exploration/result-{i}.yaml', 'r') as f:
                data = yaml.safe_load(f)
                if 'modules' in data:
                    all_modules.extend(data['modules'])
        except FileNotFoundError:
            print(f"Warning: result-{i}.yaml not found")
    
    return all_modules

def normalize_name(name):
    """Normalize resource names for comparison"""
    # Convert to lowercase and remove common variations
    normalized = name.lower().replace('_', '').replace('-', '').replace(' ', '')
    
    # Handle common name variations
    variations = {
        'biogrid': 'biogrid',
        'cellphonedb': 'cellphonedb',
        'corum': 'corum', 
        'intact': 'intact',
        'string': 'string',
        'signor': 'signor',
        'dgidb': 'dgidb',
        'kegg': 'kegg',
        'reactome': 'reactome',
        'drugbank': 'drugbank',
        'guide2pharma': 'guidetopharmacology',
        'guidetopharmacology': 'guidetopharmacology',
        'pathwaycommons': 'pathwaycommons',
        'complexportal': 'complexportal',
        'humap': 'humap',
        'icellnet': 'icellnet',
        'cellchatdb': 'cellchatdb',
        'cellinker': 'cellinker',
        'spike': 'spike',
        'compleat': 'compleat',
        'pdb': 'pdb',
        'havugimana': 'havugimana',
    }
    
    return variations.get(normalized, normalized)

def main():
    print("PyPath Module Categorization Validation")
    print("=" * 50)
    
    # Load PyPath source definitions
    print("\n1. Loading PyPath source resource definitions...")
    interaction_resources = extract_interaction_resources()
    complex_resources = extract_complex_resources()
    annotation_resources = extract_annotation_resources()
    
    print(f"   Interaction resources found: {len(interaction_resources)}")
    print(f"   Complex resources found: {len(complex_resources)}")
    print(f"   Annotation resources found: {len(annotation_resources)}")
    
    # Load our categorization results
    print("\n2. Loading our categorization results...")
    our_modules = load_categorization_results()
    print(f"   Modules categorized: {len(our_modules)}")
    
    # Create mapping of our modules
    our_module_names = {normalize_name(m['module_name']) for m in our_modules}
    our_interaction_modules = set()
    our_complex_modules = set()
    our_annotation_modules = set()
    
    for module in our_modules:
        name = normalize_name(module['module_name'])
        tags = module.get('tags', {})
        
        # Check interaction types
        interactions = tags.get('interaction', [])
        if interactions:
            our_interaction_modules.add(name)
        
        # Check entity types for complexes
        entities = tags.get('entity', [])
        if 'protein_complex' in entities:
            our_complex_modules.add(name)
        
        # Check annotation types
        annotations = tags.get('annotation', [])
        if annotations:
            our_annotation_modules.add(name)
    
    # Normalize PyPath resource names
    pypath_interaction_normalized = {normalize_name(r) for r in interaction_resources}
    pypath_complex_normalized = {normalize_name(r) for r in complex_resources}
    pypath_annotation_normalized = {normalize_name(r) for r in annotation_resources}
    
    print("\n3. Validation Results:")
    print("-" * 30)
    
    # Check interaction resources
    print("\nINTERACTION RESOURCES:")
    missing_interactions = pypath_interaction_normalized - our_module_names
    extra_interactions = our_interaction_modules - pypath_interaction_normalized
    matched_interactions = pypath_interaction_normalized & our_interaction_modules
    
    print(f"  PyPath defines: {len(pypath_interaction_normalized)} interaction resources")
    print(f"  We categorized: {len(our_interaction_modules)} modules as having interactions")
    print(f"  Matched: {len(matched_interactions)}")
    print(f"  Missing from our analysis: {len(missing_interactions)}")
    print(f"  Extra in our analysis: {len(extra_interactions)}")
    
    if missing_interactions:
        print(f"  Missing: {sorted(missing_interactions)}")
    if extra_interactions:
        print(f"  Extra: {sorted(extra_interactions)}")
    
    # Check complex resources
    print("\nCOMPLEX RESOURCES:")
    missing_complexes = pypath_complex_normalized - our_module_names
    extra_complexes = our_complex_modules - pypath_complex_normalized
    matched_complexes = pypath_complex_normalized & our_complex_modules
    
    print(f"  PyPath defines: {len(pypath_complex_normalized)} complex resources")
    print(f"  We categorized: {len(our_complex_modules)} modules as having complexes")
    print(f"  Matched: {len(matched_complexes)}")
    print(f"  Missing from our analysis: {len(missing_complexes)}")
    print(f"  Extra in our analysis: {len(extra_complexes)}")
    
    if missing_complexes:
        print(f"  Missing: {sorted(missing_complexes)}")
    if extra_complexes:
        print(f"  Extra: {sorted(extra_complexes)}")
    
    # Check annotation resources
    print("\nANNOTATION RESOURCES:")
    missing_annotations = pypath_annotation_normalized - our_module_names
    extra_annotations = our_annotation_modules - pypath_annotation_normalized
    matched_annotations = pypath_annotation_normalized & our_annotation_modules
    
    print(f"  PyPath defines: {len(pypath_annotation_normalized)} annotation resources")
    print(f"  We categorized: {len(our_annotation_modules)} modules as having annotations")
    print(f"  Matched: {len(matched_annotations)}")
    print(f"  Missing from our analysis: {len(missing_annotations)}")
    print(f"  Extra in our analysis: {len(extra_annotations)}")
    
    if missing_annotations:
        print(f"  Missing: {sorted(missing_annotations)}")
    if extra_annotations:
        print(f"  Extra: {sorted(extra_annotations)}")
    
    # Overall coverage
    print("\n4. Overall Coverage:")
    print("-" * 20)
    all_pypath_resources = pypath_interaction_normalized | pypath_complex_normalized | pypath_annotation_normalized
    covered_resources = our_module_names & all_pypath_resources
    coverage_percentage = (len(covered_resources) / len(all_pypath_resources)) * 100 if all_pypath_resources else 0
    
    print(f"  Total unique PyPath resources: {len(all_pypath_resources)}")
    print(f"  Resources we analyzed: {len(covered_resources)}")
    print(f"  Coverage: {coverage_percentage:.1f}%")
    
    # Save validation results
    validation_results = {
        'pypath_resources': {
            'interaction': sorted(interaction_resources),
            'complex': sorted(complex_resources),
            'annotation': sorted(annotation_resources),
        },
        'our_categorization': {
            'total_modules': len(our_modules),
            'interaction_modules': len(our_interaction_modules),
            'complex_modules': len(our_complex_modules),
            'annotation_modules': len(our_annotation_modules),
        },
        'validation': {
            'coverage_percentage': round(coverage_percentage, 1),
            'missing_interactions': sorted(missing_interactions),
            'missing_complexes': sorted(missing_complexes),
            'missing_annotations': sorted(missing_annotations),
            'extra_interactions': sorted(extra_interactions),
            'extra_complexes': sorted(extra_complexes),
            'extra_annotations': sorted(extra_annotations),
        }
    }
    
    with open('/Users/jschaul/Code/pypath-1/categorization_validation.yaml', 'w') as f:
        yaml.dump(validation_results, f, default_flow_style=False)
    
    print(f"\nValidation results saved to 'categorization_validation.yaml'")

if __name__ == '__main__':
    main()