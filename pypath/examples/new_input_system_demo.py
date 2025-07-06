#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Demo of the new input system with YAML configuration and bronze layer.
"""

import pypath.inputs.input_module as input_module
import pypath.inputs.biogrid_new as biogrid_new


def demo_basic_usage():
    """Basic usage of the new input system."""
    print("=== Basic Input Module Usage ===")
    
    # Method 1: Using the InputModule directly
    module = input_module.InputModule('biogrid', 'biogrid_physical')
    
    # Check if data needs downloading
    if not module.bronze_storage.exists('biogrid_physical', module.source_config):
        print("Data not in bronze layer, will download...")
    else:
        print("Data exists in bronze layer")
    
    # Load processed data (automatically handles download/caching/bronze)
    interactions = module.load_processed()
    print(f"Loaded {len(interactions)} interactions")
    
    # Show first few interactions
    for i, interaction in enumerate(interactions[:3]):
        print(f"  {i+1}. {interaction.partner_a} <-> {interaction.partner_b} (PMID: {interaction.pmid})")


def demo_compatibility():
    """Demo backward compatibility with existing code."""
    print("\n=== Backward Compatibility ===")
    
    # Use the migrated function with the same API as before
    interactions = biogrid_new.biogrid_interactions(
        organism=9606,
        htp_limit=1,
        ltp=True
    )
    
    print(f"Loaded {len(interactions)} human low-throughput interactions")
    
    # The returned data has the same structure as before
    if interactions:
        first = interactions[0]
        print(f"First interaction: {first.partner_a} <-> {first.partner_b}")


def demo_bronze_layer():
    """Demo bronze layer functionality."""
    print("\n=== Bronze Layer Demo ===")
    
    # Create module
    module = input_module.InputModule('hgnc', 'hgnc_genegroups')
    
    # Load raw data as Arrow table
    arrow_table = module.load_raw()
    print(f"Raw data shape: {arrow_table.num_rows} rows, {arrow_table.num_columns} columns")
    print(f"Columns: {arrow_table.column_names}")
    
    # Convert to pandas for inspection
    df = arrow_table.to_pandas()
    print(f"\nFirst few rows:")
    print(df.head())
    
    # List all bronze sources
    bronze_storage = module.bronze_storage
    sources = bronze_storage.list_sources()
    print(f"\nBronze layer contains {len(sources)} sources:")
    for src in sources:
        print(f"  - {src['source_name']}: {src['num_rows']} rows, saved at {src['saved_at']}")


def demo_change_detection():
    """Demo change detection functionality."""
    print("\n=== Change Detection Demo ===")
    
    # This would work with the extended download manager
    from download_manager._change_detection import ChangeDetector
    
    # Example URL to check
    url = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
    
    # Get headers without downloading
    detector = ChangeDetector()
    headers = detector.get_headers(url)
    
    print(f"Resource headers:")
    for key, value in headers.items():
        if value:
            print(f"  {key}: {value}")


def demo_custom_config():
    """Demo creating custom YAML configuration."""
    print("\n=== Custom Configuration Demo ===")
    
    # Create a custom config dict (normally this would be in a YAML file)
    custom_config = {
        'my_custom_source': {
            'url': 'https://example.com/data.tsv',
            'format': 'tsv',
            'field_mapping': {
                'id': 0,
                'name': 1,
                'value': 2,
            },
            'filters': [
                {'field': 'value', 'operator': 'gt', 'value': 0}
            ],
            'description': 'My custom data source',
            'organism': 9606,
        }
    }
    
    # Create module with custom config
    module = input_module.InputModule(custom_config, 'my_custom_source')
    print(f"Created module for: {module.source_config['description']}")
    print(f"URL: {module.source_config['url']}")
    print(f"Fields: {list(module.source_config['field_mapping'].keys())}")


def demo_factory_function():
    """Demo using factory function for simple cases."""
    print("\n=== Factory Function Demo ===")
    
    # Create a simple input function
    load_biogrid = input_module.create_input_function('biogrid', 'biogrid_physical')
    
    # Use it like any other input function
    interactions = load_biogrid(organism=9606)
    print(f"Factory function loaded {len(interactions)} interactions")


if __name__ == '__main__':
    print("PyPath New Input System Demo")
    print("=" * 50)
    
    # Run demos
    demo_basic_usage()
    demo_compatibility()
    demo_bronze_layer()
    demo_change_detection()
    demo_custom_config()
    demo_factory_function()
    
    print("\n" + "=" * 50)
    print("Demo completed!")