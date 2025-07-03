#!/usr/bin/env python3
"""
PyPath Input Module Testing Script

This script provides a systematic way to test input modules from pypath.inputs.
It can test individual modules, specific functions, or run batch tests.

Usage:
    python test_input_modules.py --module abs
    python test_input_modules.py --modules abs,acsn,adhesome
    python test_input_modules.py --function abs.abs_interactions
    python test_input_modules.py --batch-test first5
    python test_input_modules.py --list-modules
"""

import sys
import os
import argparse
import importlib
import inspect
import traceback
import json
from typing import Dict, List, Any, Optional
from datetime import datetime

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Dynamic module loading
def load_modules_config(config_file: str = 'input_module_maintenance/modules.json') -> Dict[str, List[str]]:
    """Load module configuration from JSON file."""
    try:
        with open(config_file, 'r') as f:
            data = json.load(f)
        return data.get('modules', {})
    except FileNotFoundError:
        print(f"Warning: Configuration file {config_file} not found. Run discover_input_modules.py first.")
        return {}
    except Exception as e:
        print(f"Error loading configuration: {e}")
        return {}

class InputModuleTester:
    """Test runner for PyPath input modules."""
    
    def __init__(self, verbose: bool = True, timeout: int = 120):
        self.verbose = verbose
        self.timeout = timeout
        self.results = []
        
    def log(self, message: str, level: str = 'INFO'):
        """Log message with timestamp."""
        if self.verbose:
            timestamp = datetime.now().strftime('%H:%M:%S')
            print(f"[{timestamp}] {level}: {message}")
    
    def test_function(self, module_name: str, function_name: str) -> Dict[str, Any]:
        """Test a single function from a module."""
        result = {
            'module': module_name,
            'function': function_name,
            'status': 'UNKNOWN',
            'error': None,
            'record_count': 0,
            'sample_record': None,
            'execution_time': 0
        }
        
        try:
            self.log(f"Testing {module_name}.{function_name}()")
            
            # Import module
            module_path = f'pypath.inputs.{module_name}'
            module = importlib.import_module(module_path)
            
            # Get function
            if not hasattr(module, function_name):
                raise AttributeError(f"Function '{function_name}' not found in module '{module_name}'")
            
            func = getattr(module, function_name)
            
            # Execute function with timing
            start_time = datetime.now()
            data = func()
            end_time = datetime.now()
            
            result['execution_time'] = (end_time - start_time).total_seconds()
            
            # Analyze results
            if data is None:
                result['status'] = 'SUCCESS'
                result['record_count'] = 0
                self.log(f"  → Function returned None")
                
            elif isinstance(data, dict):
                result['status'] = 'SUCCESS'
                result['record_count'] = len(data)
                
                if result['record_count'] > 0:
                    # For dictionaries, show first key-value pair
                    first_key = next(iter(data.keys()))
                    sample_value = str(data[first_key])[:100] + "..." if len(str(data[first_key])) > 100 else str(data[first_key])
                    result['sample_record'] = f"{first_key}: {sample_value}"
                    self.log(f"  → SUCCESS: {result['record_count']} dictionary entries returned")
                    if self.verbose:
                        self.log(f"  → Sample: {result['sample_record']}")
                else:
                    self.log(f"  → SUCCESS: Empty dictionary returned")
                    
            elif isinstance(data, set):
                result['status'] = 'SUCCESS'
                result['record_count'] = len(data)
                
                if result['record_count'] > 0:
                    # For sets, show first element
                    first_item = next(iter(data))
                    result['sample_record'] = str(first_item)[:200] + "..." if len(str(first_item)) > 200 else str(first_item)
                    self.log(f"  → SUCCESS: {result['record_count']} set items returned")
                    if self.verbose:
                        self.log(f"  → Sample: {result['sample_record']}")
                else:
                    self.log(f"  → SUCCESS: Empty set returned")
                    
            elif hasattr(data, '__len__') and not isinstance(data, str):
                result['status'] = 'SUCCESS'
                result['record_count'] = len(data)
                
                if result['record_count'] > 0:
                    result['sample_record'] = str(data[0])[:200] + "..." if len(str(data[0])) > 200 else str(data[0])
                    self.log(f"  → SUCCESS: {result['record_count']} records returned")
                    if self.verbose:
                        self.log(f"  → Sample: {result['sample_record']}")
                else:
                    self.log(f"  → SUCCESS: Empty result returned")
                    
            else:
                result['status'] = 'SUCCESS'
                result['record_count'] = 1
                result['sample_record'] = str(data)[:200] + "..." if len(str(data)) > 200 else str(data)
                self.log(f"  → SUCCESS: Single object returned")
                if self.verbose:
                    self.log(f"  → Data: {result['sample_record']}")
                
        except Exception as e:
            result['status'] = 'ERROR'
            result['error'] = str(e)
            self.log(f"  → ERROR: {e}", 'ERROR')
            if self.verbose:
                self.log(f"  → Traceback: {traceback.format_exc()}", 'ERROR')
        
        self.results.append(result)
        return result
    
    def test_module(self, module_name: str, functions: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """Test all functions in a module or specific functions."""
        self.log(f"Testing module: {module_name}")
        
        try:
            # Import module to get available functions
            module_path = f'pypath.inputs.{module_name}'
            module = importlib.import_module(module_path)
            
            if functions is None:
                # Get all public functions
                functions = [name for name, obj in inspect.getmembers(module) 
                           if inspect.isfunction(obj) and not name.startswith('_')]
                # Filter out non-module functions (like imported ones)
                functions = [f for f in functions if hasattr(getattr(module, f), '__module__') 
                           and getattr(module, f).__module__ == module_path]
            
            if not functions:
                self.log(f"  → No testable functions found in {module_name}")
                return []
            
            self.log(f"  → Found functions: {functions}")
            
            results = []
            for func_name in functions:
                result = self.test_function(module_name, func_name)
                results.append(result)
                
            return results
            
        except Exception as e:
            self.log(f"  → ERROR importing module {module_name}: {e}", 'ERROR')
            return []
    
    def run_batch_test(self, batch_name: str) -> List[Dict[str, Any]]:
        """Run a predefined batch test."""
        modules_config = load_modules_config()
        if not modules_config:
            raise ValueError("No modules configuration found. Run discover_input_modules.py first.")
        
        if batch_name == 'first5':
            # First 5 modules alphabetically
            config = dict(list(sorted(modules_config.items()))[:5])
        elif batch_name == 'interactions':
            # Modules with interaction functions
            config = {}
            for module, functions in modules_config.items():
                interaction_funcs = [f for f in functions if 'interactions' in f]
                if interaction_funcs:
                    config[module] = interaction_funcs
        elif batch_name == 'annotations':
            # Modules with annotation functions
            config = {}
            for module, functions in modules_config.items():
                annotation_funcs = [f for f in functions if 'annotations' in f]
                if annotation_funcs:
                    config[module] = annotation_funcs
        elif batch_name == 'all':
            # All modules
            config = modules_config
        else:
            raise ValueError(f"Unknown batch test '{batch_name}'. Available: first5, interactions, annotations, all")
        
        self.log(f"Running batch test: {batch_name}")
        self.log(f"Testing modules: {list(config.keys())}")
        
        all_results = []
        for module_name, functions in config.items():
            results = self.test_module(module_name, functions)
            all_results.extend(results)
        
        return all_results
    
    def print_summary(self):
        """Print test summary."""
        if not self.results:
            print("No tests run.")
            return
        
        successful = [r for r in self.results if r['status'] == 'SUCCESS']
        failed = [r for r in self.results if r['status'] == 'ERROR']
        
        print(f"\n{'='*60}")
        print(f"TEST SUMMARY")
        print(f"{'='*60}")
        print(f"Total tests: {len(self.results)}")
        print(f"Successful: {len(successful)}")
        print(f"Failed: {len(failed)}")
        print(f"Success rate: {len(successful)/len(self.results)*100:.1f}%")
        
        if successful:
            print(f"\n✅ SUCCESSFUL TESTS:")
            for result in successful:
                print(f"  • {result['module']}.{result['function']}() - {result['record_count']} records ({result['execution_time']:.1f}s)")
        
        if failed:
            print(f"\n❌ FAILED TESTS:")
            for result in failed:
                print(f"  • {result['module']}.{result['function']}() - {result['error']}")
    
    def list_available_modules(self):
        """List all available input modules."""
        inputs_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'pypath', 'inputs')
        
        modules = []
        for file in os.listdir(inputs_dir):
            if file.endswith('.py') and not file.startswith('_'):
                module_name = file[:-3]
                modules.append(module_name)
        
        print(f"Available input modules ({len(modules)}):")
        for i, module in enumerate(sorted(modules), 1):
            print(f"{i:3d}. {module}")
        
        return modules

def main():
    parser = argparse.ArgumentParser(description='Test PyPath input modules')
    parser.add_argument('--module', '-m', help='Test specific module')
    parser.add_argument('--modules', help='Test multiple modules (comma-separated)')
    parser.add_argument('--function', '-f', help='Test specific function (format: module.function)')
    parser.add_argument('--batch-test', '-b', help='Run batch test', choices=['first5', 'interactions', 'annotations', 'all'])
    parser.add_argument('--list-modules', '-l', action='store_true', help='List available modules')
    parser.add_argument('--quiet', '-q', action='store_true', help='Quiet mode')
    parser.add_argument('--timeout', '-t', type=int, default=120, help='Timeout for each test (seconds)')
    
    args = parser.parse_args()
    
    tester = InputModuleTester(verbose=not args.quiet, timeout=args.timeout)
    
    if args.list_modules:
        tester.list_available_modules()
        return
    
    try:
        if args.function:
            # Test specific function
            if '.' not in args.function:
                print("Error: Function must be in format 'module.function'")
                return
            module_name, function_name = args.function.split('.', 1)
            tester.test_function(module_name, function_name)
            
        elif args.module:
            # Test specific module
            tester.test_module(args.module)
            
        elif args.modules:
            # Test multiple modules
            modules = [m.strip() for m in args.modules.split(',')]
            for module in modules:
                tester.test_module(module)
                
        elif args.batch_test:
            # Run batch test
            tester.run_batch_test(args.batch_test)
            
        else:
            # Default: run first5 batch test
            print("No specific test specified. Running 'first5' batch test...")
            tester.run_batch_test('first5')
        
        tester.print_summary()
        
    except KeyboardInterrupt:
        print("\nTest interrupted by user")
    except Exception as e:
        print(f"Error: {e}")
        if not args.quiet:
            traceback.print_exc()

if __name__ == '__main__':
    main()