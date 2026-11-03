#!/usr/bin/env python
"""
Experiment: Mapping FooDB and Phenol-Explorer foods to FoodOn ontology.

This script explores different approaches to map food terms from these databases
to FoodOn (the Food Ontology) terms.

Approaches:
1. NCBI Taxonomy ID matching - both FooDB and FoodOn use NCBI taxonomy
2. Exact name matching (case-insensitive)
3. Scientific name matching (case-insensitive)

Run this script from the pypath directory:
    python experiment.py
"""

from __future__ import annotations

import csv
import io
import os
import zipfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import requests

# Data paths
PYPATH_DATA = Path(__file__).parent / 'pypath-data'
FOODB_DIR = PYPATH_DATA / 'foodb' / 'foodb_2020_04_07_csv'
PHENOL_EXPLORER_DIR = PYPATH_DATA / 'phenol_explorer'


# =============================================================================
# Data Loading
# =============================================================================

def load_foodb_foods() -> pd.DataFrame:
    """Load FooDB Food.csv data."""
    food_path = FOODB_DIR / 'Food.csv'
    if not food_path.exists():
        print(f"FooDB Food.csv not found at {food_path}")
        return pd.DataFrame()
    
    df = pd.read_csv(food_path, low_memory=False)
    print(f"\nLoaded {len(df)} foods from FooDB")
    
    # Show sample columns
    print(f"Columns: {list(df.columns)}")
    
    return df


def load_phenol_explorer_foods() -> pd.DataFrame:
    """Load Phenol-Explorer foods.csv data."""
    foods_zip = PHENOL_EXPLORER_DIR / 'foods.csv.zip'
    if not foods_zip.exists():
        print(f"Phenol-Explorer foods.csv.zip not found at {foods_zip}")
        return pd.DataFrame()
    
    with zipfile.ZipFile(foods_zip, 'r') as zf:
        # Find the CSV file inside
        csv_names = [n for n in zf.namelist() if n.endswith('.csv')]
        if not csv_names:
            print("No CSV file found in foods.csv.zip")
            return pd.DataFrame()
        
        with zf.open(csv_names[0]) as f:
            df = pd.read_csv(io.TextIOWrapper(f, encoding='utf-8'))
    
    print(f"\nLoaded {len(df)} foods from Phenol-Explorer")
    print(f"Columns: {list(df.columns)}")
    
    return df


# =============================================================================
# FoodOn Loading
# =============================================================================

def download_foodon_terms() -> pd.DataFrame:
    """
    Download FoodOn terms from GitHub OWL file.
    
    Returns a DataFrame with FoodOn term IDs, labels, and synonyms.
    """
    cache_path = PYPATH_DATA / 'foodon_terms.csv'
    owl_cache_path = PYPATH_DATA / 'foodon-full.owl'
    
    if cache_path.exists():
        print(f"\nLoading cached FoodOn terms from {cache_path}")
        return pd.read_csv(cache_path)
    
    print("\nDownloading FoodOn ontology...")
    
    # Download OWL file if not cached
    if not owl_cache_path.exists():
        url = "https://github.com/FoodOntology/foodon/raw/refs/heads/master/foodon-full.owl"
        print(f"  Downloading from: {url}")
        
        try:
            response = requests.get(url, timeout=300, allow_redirects=True, stream=True)
            response.raise_for_status()
            
            # Save to cache
            with open(owl_cache_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            print(f"  Downloaded to: {owl_cache_path}")
            print(f"  File size: {owl_cache_path.stat().st_size / 1024 / 1024:.1f} MB")
            
        except Exception as e:
            print(f"  Download failed: {e}")
            return _get_fallback_foodon_terms()
    else:
        print(f"  Using cached OWL file: {owl_cache_path}")
    
    # Parse OWL with rdflib
    print("\nParsing FoodOn OWL file (this may take a minute)...")
    
    try:
        from rdflib import Graph, Namespace, URIRef
        from rdflib.namespace import RDF, RDFS, OWL
        
        g = Graph()
        g.parse(str(owl_cache_path), format='xml')
        
        print(f"  Loaded {len(g)} triples")
        
        # Namespaces
        FOODON = Namespace("http://purl.obolibrary.org/obo/FOODON_")
        OBO = Namespace("http://purl.obolibrary.org/obo/")
        OBOINOWL = Namespace("http://www.geneontology.org/formats/oboInOwl#")
        NCBITAXON = Namespace("http://purl.obolibrary.org/obo/NCBITaxon_")
        
        # Extract FoodOn terms
        rows = []
        taxon_count = 0
        
        # Find all classes with FOODON IDs
        for s in g.subjects(RDF.type, OWL.Class):
            uri = str(s)
            if 'FOODON_' not in uri:
                continue
            
            # Extract FOODON ID
            foodon_id = uri.split('/')[-1].replace('_', ':')
            
            # Get label
            label = None
            for o in g.objects(s, RDFS.label):
                label = str(o)
                break
            
            # Get synonyms (hasExactSynonym, hasRelatedSynonym, hasBroadSynonym)
            synonyms = []
            for syn_pred in [OBOINOWL.hasExactSynonym, OBOINOWL.hasRelatedSynonym, OBOINOWL.hasBroadSynonym]:
                for o in g.objects(s, syn_pred):
                    synonyms.append(str(o))
            
            # Get NCBI taxon - multiple approaches:
            ncbi_taxon = ''
            
            # Approach 1: hasDbXref with NCBITaxon: prefix
            for o in g.objects(s, OBOINOWL.hasDbXref):
                xref = str(o)
                if xref.startswith('NCBITaxon:'):
                    ncbi_taxon = xref.split(':')[1]
                    break
            
            # Approach 2: Check rdfs:subClassOf pointing to NCBITaxon
            if not ncbi_taxon:
                for o in g.objects(s, RDFS.subClassOf):
                    parent_uri = str(o)
                    if 'NCBITaxon_' in parent_uri:
                        # Extract taxon ID from URI like http://purl.obolibrary.org/obo/NCBITaxon_12345
                        ncbi_taxon = parent_uri.split('NCBITaxon_')[-1]
                        break
            
            # Approach 3: Check for restriction on NCBITaxon
            if not ncbi_taxon:
                for restriction in g.objects(s, RDFS.subClassOf):
                    # Check if it's a restriction
                    if (restriction, RDF.type, OWL.Restriction) in g:
                        for on_class in g.objects(restriction, OWL.someValuesFrom):
                            on_uri = str(on_class)
                            if 'NCBITaxon_' in on_uri:
                                ncbi_taxon = on_uri.split('NCBITaxon_')[-1]
                                break
            
            if ncbi_taxon:
                taxon_count += 1
            
            if label:
                rows.append({
                    'foodon_id': foodon_id,
                    'name': label,
                    'synonyms': '|'.join(synonyms),
                    'ncbi_taxon': ncbi_taxon,
                })
        
        df = pd.DataFrame(rows)
        
        # Cache for next time
        if not df.empty:
            df.to_csv(cache_path, index=False)
            print(f"\nExtracted and cached {len(df)} FoodOn terms")
            print(f"  Terms with NCBI taxon: {taxon_count}")
        
        return df
        
    except Exception as e:
        print(f"  Error parsing OWL: {e}")
        import traceback
        traceback.print_exc()
        return _get_fallback_foodon_terms()


def _get_fallback_foodon_terms() -> pd.DataFrame:
    """Return a minimal set of FoodOn terms for testing."""
    print("Using minimal built-in FoodOn sample for testing...")
    return pd.DataFrame([
        {'foodon_id': 'FOODON:00001002', 'name': 'apple', 'synonyms': 'Malus domestica', 'ncbi_taxon': '3750'},
        {'foodon_id': 'FOODON:00001015', 'name': 'orange', 'synonyms': 'Citrus sinensis', 'ncbi_taxon': '2711'},
        {'foodon_id': 'FOODON:03301710', 'name': 'banana', 'synonyms': 'Musa', 'ncbi_taxon': '4641'},
        {'foodon_id': 'FOODON:03310788', 'name': 'tomato', 'synonyms': 'Solanum lycopersicum', 'ncbi_taxon': '4081'},
        {'foodon_id': 'FOODON:03301103', 'name': 'potato', 'synonyms': 'Solanum tuberosum', 'ncbi_taxon': '4113'},
        {'foodon_id': 'FOODON:03301274', 'name': 'carrot', 'synonyms': 'Daucus carota', 'ncbi_taxon': '79200'},
        {'foodon_id': 'FOODON:00001256', 'name': 'beef', 'synonyms': 'Bos taurus|cattle meat', 'ncbi_taxon': '9913'},
        {'foodon_id': 'FOODON:00001283', 'name': 'chicken', 'synonyms': 'Gallus gallus|poultry', 'ncbi_taxon': '9031'},
        {'foodon_id': 'FOODON:03411457', 'name': 'milk', 'synonyms': 'cow milk|dairy milk', 'ncbi_taxon': ''},
        {'foodon_id': 'FOODON:00002473', 'name': 'wheat', 'synonyms': 'Triticum aestivum', 'ncbi_taxon': '4565'},
    ])


# =============================================================================
# Name Normalization
# =============================================================================

def normalize_name(name: str) -> str:
    """Normalize a food name for matching."""
    if pd.isna(name):
        return ''
    return str(name).lower().strip()


# =============================================================================
# Mapping Approaches
# =============================================================================


def map_by_ncbi_taxonomy(
    foods_df: pd.DataFrame,
    foodon_df: pd.DataFrame,
    taxon_col: str = 'ncbi_taxonomy_id',
) -> dict[str, str]:
    """
    Map foods to FoodOn using NCBI taxonomy IDs.
    
    Returns: dict mapping food name -> FoodOn ID
    """
    if taxon_col not in foods_df.columns:
        print(f"  Column '{taxon_col}' not found")
        return {}
    
    # Build FoodOn taxon lookup
    foodon_by_taxon = {}
    for _, row in foodon_df.iterrows():
        if pd.notna(row.get('ncbi_taxon')) and row['ncbi_taxon']:
            taxon = str(int(float(row['ncbi_taxon']))) if row['ncbi_taxon'] else ''
            if taxon:
                foodon_by_taxon[taxon] = (row['foodon_id'], row['name'])
    
    print(f"  FoodOn terms with NCBI taxon: {len(foodon_by_taxon)}")
    
    # Match
    mappings = {}
    for _, row in foods_df.iterrows():
        food_name = row.get('name', '')
        taxon = row.get(taxon_col)
        
        if pd.notna(taxon):
            taxon_str = str(int(float(taxon))) if taxon else ''
            if taxon_str in foodon_by_taxon:
                foodon_id, foodon_name = foodon_by_taxon[taxon_str]
                mappings[food_name] = {
                    'foodon_id': foodon_id,
                    'foodon_name': foodon_name,
                    'method': 'ncbi_taxonomy',
                }
    
    return mappings


def map_by_exact_name(
    foods_df: pd.DataFrame,
    foodon_df: pd.DataFrame,
    name_col: str = 'name',
) -> dict[str, dict]:
    """
    Map foods to FoodOn using exact name matching (case-insensitive).
    """
    # Build FoodOn name lookup (including synonyms)
    foodon_by_name = {}
    for _, row in foodon_df.iterrows():
        names = [normalize_name(row.get('name', ''))]
        if pd.notna(row.get('synonyms')):
            names.extend([normalize_name(s) for s in str(row['synonyms']).split('|')])
        
        for name in names:
            if name:
                foodon_by_name[name] = (row['foodon_id'], row.get('name', ''))
    
    print(f"  FoodOn terms + synonyms: {len(foodon_by_name)} names")
    
    # Match
    mappings = {}
    for _, row in foods_df.iterrows():
        food_name = row.get(name_col, '')
        normalized = normalize_name(food_name)
        
        if normalized in foodon_by_name:
            foodon_id, foodon_name = foodon_by_name[normalized]
            mappings[food_name] = {
                'foodon_id': foodon_id,
                'foodon_name': foodon_name,
                'method': 'exact_name',
            }
    
    return mappings


def map_by_scientific_name(
    foods_df: pd.DataFrame,
    foodon_df: pd.DataFrame,
    sci_name_col: str = 'name_scientific',
) -> dict[str, dict]:
    """
    Map foods to FoodOn using scientific name matching.
    """
    if sci_name_col not in foods_df.columns:
        print(f"  Column '{sci_name_col}' not found")
        return {}
    
    # Build FoodOn name lookup
    foodon_by_name = {}
    for _, row in foodon_df.iterrows():
        names = [normalize_name(row.get('name', ''))]
        if pd.notna(row.get('synonyms')):
            names.extend([normalize_name(s) for s in str(row['synonyms']).split('|')])
        
        for name in names:
            if name:
                foodon_by_name[name] = (row['foodon_id'], row.get('name', ''))
    
    # Match
    mappings = {}
    for _, row in foods_df.iterrows():
        food_name = row.get('name', '')
        sci_name = row.get(sci_name_col, '')
        normalized = normalize_name(sci_name)
        
        if normalized and normalized in foodon_by_name:
            foodon_id, foodon_name = foodon_by_name[normalized]
            mappings[food_name] = {
                'foodon_id': foodon_id,
                'foodon_name': foodon_name,
                'method': 'scientific_name',
            }
    
    return mappings


# =============================================================================
# Analysis
# =============================================================================

def analyze_mappings(
    source_name: str,
    total_foods: int,
    *mapping_results: tuple[str, dict],
):
    """Analyze and report mapping results."""
    print(f"\n{'='*60}")
    print(f"Mapping Results for {source_name}")
    print(f"{'='*60}")
    print(f"Total foods: {total_foods}")
    
    # Combine all mappings (later methods don't override earlier ones)
    combined = {}
    for method_name, mappings in mapping_results:
        new_mappings = 0
        for food, data in mappings.items():
            if food not in combined:
                combined[food] = data
                new_mappings += 1
        print(f"\n{method_name}:")
        print(f"  Matched: {len(mappings)} ({100*len(mappings)/total_foods:.1f}%)")
        print(f"  New (not matched by earlier methods): {new_mappings}")
    
    print(f"\n{'='*60}")
    print(f"COMBINED COVERAGE: {len(combined)}/{total_foods} ({100*len(combined)/total_foods:.1f}%)")
    print(f"{'='*60}")
    
    # Show some examples
    print("\nSample mappings:")
    for i, (food, data) in enumerate(list(combined.items())[:5]):
        print(f"  {food} -> {data['foodon_name']} ({data['method']})")
    
    # Show unmapped examples
    unmapped = [row['name'] for _, row in foods_df.iterrows() if row['name'] not in combined]
    if unmapped:
        print(f"\nSample UNMAPPED foods ({len(unmapped)} total):")
        for food in unmapped[:5]:
            print(f"  - {food}")
    
    return combined


# =============================================================================
# Main
# =============================================================================

if __name__ == '__main__':
    print("="*60)
    print("FoodOn Mapping Experiment")
    print("="*60)
    
    # Load data
    foodb_df = load_foodb_foods()
    phenol_df = load_phenol_explorer_foods()
    foodon_df = download_foodon_terms()
    
    if foodon_df.empty:
        print("\nFailed to load FoodOn terms. Exiting.")
        exit(1)
    
    print(f"\nFoodOn terms loaded: {len(foodon_df)}")
    print(f"Sample FoodOn terms:")
    print(foodon_df.head())
    
    # =========================================================================
    # FooDB Mapping
    # =========================================================================
    if not foodb_df.empty:
        print("\n" + "="*60)
        print("Mapping FooDB to FoodOn")
        print("="*60)
        
        # Check available columns
        print(f"\nFooDB has taxonomy column: {'ncbi_taxonomy_id' in foodb_df.columns}")
        if 'ncbi_taxonomy_id' in foodb_df.columns:
            has_taxon = foodb_df['ncbi_taxonomy_id'].notna().sum()
            print(f"Foods with NCBI taxonomy: {has_taxon}/{len(foodb_df)}")
        
        # Run mapping approaches
        print("\n--- Approach 1: NCBI Taxonomy ID ---")
        taxon_maps = map_by_ncbi_taxonomy(foodb_df, foodon_df)
        
        print("\n--- Approach 2: Exact Name Match ---")
        exact_maps = map_by_exact_name(foodb_df, foodon_df)
        
        print("\n--- Approach 3: Scientific Name Match ---")
        sci_maps = map_by_scientific_name(foodb_df, foodon_df)
        
        # Analyze
        foods_df = foodb_df
        analyze_mappings(
            "FooDB",
            len(foodb_df),
            ("NCBI Taxonomy", taxon_maps),
            ("Exact Name", exact_maps),
            ("Scientific Name", sci_maps),
        )
    
    # =========================================================================
    # Phenol-Explorer Mapping (with preprocessing)
    # =========================================================================
    if not phenol_df.empty:
        print("\n" + "="*60)
        print("Mapping Phenol-Explorer to FoodOn")
        print("="*60)
        
        # Check columns
        print(f"\nPhenol-Explorer columns: {list(phenol_df.columns)}")
        print(f"\nPhenol-Explorer has taxonomy column: {'ncbi_taxonomy_id' in phenol_df.columns}")
        if 'ncbi_taxonomy_id' in phenol_df.columns:
            has_taxon = phenol_df['ncbi_taxonomy_id'].notna().sum()
            print(f"Foods with NCBI taxonomy: {has_taxon}/{len(phenol_df)}")
        
        # Phenol-Explorer may have different column names
        sci_name_col = 'food_source_scientific_name' if 'food_source_scientific_name' in phenol_df.columns else None

        taxon_maps = {}
        if 'ncbi_taxonomy_id' in phenol_df.columns:
            print("\n--- Approach 1: NCBI Taxonomy ID ---")
            taxon_maps = map_by_ncbi_taxonomy(phenol_df, foodon_df)
        
        print("\n--- Approach 2: Exact Name Match ---")
        exact_maps = map_by_exact_name(phenol_df, foodon_df)
        
        if sci_name_col:
            print(f"\n--- Approach 3: Scientific Name Match ({sci_name_col}) ---")
            sci_maps = map_by_scientific_name(phenol_df, foodon_df, sci_name_col=sci_name_col)
        else:
            sci_maps = {}
        
        # Analyze
        foods_df = phenol_df
        results = []
        if taxon_maps:
            results.append(("NCBI Taxonomy", taxon_maps))
        results.append(("Exact Name", exact_maps))
        if sci_maps:
            results.append(("Scientific Name", sci_maps))
        
        analyze_mappings(
            "Phenol-Explorer",
            len(phenol_df),
            *results,
        )
    
    print("\n" + "="*60)
    print("Experiment complete!")
    print("="*60)
