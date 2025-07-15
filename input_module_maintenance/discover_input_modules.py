#!/usr/bin/env python3
"""
PyPath Input Module Discovery Script

This script discovers all input modules and their functions, creating a comprehensive
mapping that can be used for testing and documentation purposes.

Usage:
    python discover_input_modules.py --output modules_config.json
    python discover_input_modules.py --print-summary
    python discover_input_modules.py --filter-by-suffix interactions
"""

import os
import sys
import json
import argparse
import importlib
import inspect
from typing import Dict, List, Set
from pathlib import Path

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class InputModuleDiscovery:
    """Discover PyPath input modules and their functions."""
    
    def __init__(self, inputs_dir: str = None):
        if inputs_dir is None:
            self.inputs_dir = Path(__file__).parent.parent / 'pypath' / 'inputs'
        else:
            self.inputs_dir = Path(inputs_dir)
        
        self.modules_info = {}
        self.common_suffixes = set()
        
    def discover_modules(self) -> Dict[str, List[str]]:
        """Discover all input modules and their functions."""
        print(f"Discovering modules in: {self.inputs_dir}")
        
        # Get all Python files in inputs directory
        python_files = []
        for file_path in self.inputs_dir.rglob('*.py'):
            if file_path.name.startswith('_'):
                continue
            if file_path.is_file():
                python_files.append(file_path)
        
        print(f"Found {len(python_files)} Python files")
        
        for file_path in python_files:
            # Calculate module name from file path
            relative_path = file_path.relative_to(self.inputs_dir)
            module_parts = list(relative_path.parts[:-1]) + [relative_path.stem]
            
            if len(module_parts) == 1:
                # Simple module like abs.py
                module_name = module_parts[0]
                module_path = f'pypath.inputs.{module_name}'
            else:
                # Complex module like bindingdb/__init__.py
                module_name = '/'.join(module_parts)
                module_path = f'pypath.inputs.{".".join(module_parts)}'
            
            try:
                functions = self._get_module_functions(module_path, module_name)
                if functions:
                    self.modules_info[module_name] = functions
                    
                    # Track common suffixes
                    for func in functions:
                        if '_' in func:
                            suffix = func.split('_')[-1]
                            self.common_suffixes.add(suffix)
                            
            except Exception as e:
                print(f"Warning: Could not analyze {module_name}: {e}")
        
        return self.modules_info
    
    def _get_module_functions(self, module_path: str, module_name: str) -> List[str]:
        """Get all relevant functions from a module."""
        try:
            module = importlib.import_module(module_path)
            
            # Get all functions
            all_functions = [
                name for name, obj in inspect.getmembers(module)
                if inspect.isfunction(obj) and not name.startswith('_')
            ]
            
            # Filter for functions that belong to this module
            module_functions = []
            for func_name in all_functions:
                func_obj = getattr(module, func_name)
                if hasattr(func_obj, '__module__') and func_obj.__module__ == module_path:
                    module_functions.append(func_name)
            
            # Additional filtering: prefer functions that start with module base name
            base_name = module_name.split('/')[0]  # Get base name for complex modules
            relevant_functions = []
            
            for func_name in module_functions:
                # Include if function starts with module name or contains common patterns
                if (func_name.startswith(base_name) or 
                    any(pattern in func_name for pattern in ['interactions', 'annotations', 'raw', 'mapping', 'enz_sub', 'complexes'])):
                    relevant_functions.append(func_name)
            
            # If no relevant functions found, include all module functions
            if not relevant_functions:
                relevant_functions = module_functions
            
            return sorted(relevant_functions)
            
        except Exception as e:
            raise Exception(f"Failed to import {module_path}: {e}")
    
    def get_modules_by_suffix(self, suffix: str) -> Dict[str, List[str]]:
        """Get modules that have functions with a specific suffix."""
        filtered = {}
        for module_name, functions in self.modules_info.items():
            matching_functions = [f for f in functions if f.endswith(f'_{suffix}')]
            if matching_functions:
                filtered[module_name] = matching_functions
        return filtered
    
    def get_summary(self) -> Dict[str, any]:
        """Get summary statistics."""
        total_modules = len(self.modules_info)
        total_functions = sum(len(funcs) for funcs in self.modules_info.values())
        
        suffix_counts = {}
        for suffix in self.common_suffixes:
            suffix_counts[suffix] = len(self.get_modules_by_suffix(suffix))
        
        return {
            'total_modules': total_modules,
            'total_functions': total_functions,
            'common_suffixes': sorted(self.common_suffixes),
            'suffix_counts': suffix_counts,
            'modules_by_function_count': {
                module: len(functions) 
                for module, functions in sorted(self.modules_info.items(), key=lambda x: len(x[1]), reverse=True)
            }
        }
    
    def save_to_json(self, output_file: str):
        """Save module information to JSON file."""
        data = {
            'modules': self.modules_info,
            'summary': self.get_summary(),
            'generated_at': str(Path(__file__).stat().st_mtime)
        }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"Module information saved to: {output_file}")
    
    def print_summary(self):
        """Print discovery summary."""
        summary = self.get_summary()
        
        print(f"\n{'='*60}")
        print(f"PYPATH INPUT MODULES DISCOVERY SUMMARY")
        print(f"{'='*60}")
        print(f"Total modules discovered: {summary['total_modules']}")
        print(f"Total functions discovered: {summary['total_functions']}")
        print(f"Average functions per module: {summary['total_functions']/summary['total_modules']:.1f}")
        
        print(f"\nCommon function suffixes:")
        for suffix in summary['common_suffixes']:
            count = summary['suffix_counts'].get(suffix, 0)
            print(f"  • {suffix}: {count} modules")
        
        print(f"\nTop 10 modules by function count:")
        for i, (module, count) in enumerate(list(summary['modules_by_function_count'].items())[:10], 1):
            print(f"  {i:2d}. {module}: {count} functions")
        
        print(f"\nSample modules and their functions:")
        for i, (module, functions) in enumerate(list(self.modules_info.items())[:5], 1):
            print(f"  {i}. {module}: {functions}")

def main():
    parser = argparse.ArgumentParser(description='Discover PyPath input modules and functions')
    parser.add_argument('--output', '-o', help='Output JSON file path', default='modules_config.json')
    parser.add_argument('--print-summary', '-s', action='store_true', help='Print discovery summary')
    parser.add_argument('--filter-by-suffix', '-f', help='Filter modules by function suffix')
    parser.add_argument('--inputs-dir', help='Custom inputs directory path')
    
    args = parser.parse_args()
    
    try:
        discoverer = InputModuleDiscovery(args.inputs_dir)
        discoverer.discover_modules()
        
        if args.print_summary:
            discoverer.print_summary()
        
        if args.filter_by_suffix:
            filtered = discoverer.get_modules_by_suffix(args.filter_by_suffix)
            print(f"\nModules with '{args.filter_by_suffix}' functions:")
            for module, functions in filtered.items():
                print(f"  • {module}: {functions}")
        
        if args.output:
            discoverer.save_to_json(args.output)
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()