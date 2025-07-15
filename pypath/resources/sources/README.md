# PyPath Input Module Configuration System

This directory contains YAML configuration files for PyPath data sources. The new system simplifies adding and maintaining data sources by using declarative configuration instead of code.

## Overview

The new input system provides:

1. **YAML-based configuration** - Define data sources declaratively
2. **Bronze layer storage** - Automatic caching in Parquet format for performance
3. **Change detection** - Only re-download when sources change (ETag, Last-Modified, checksums)
4. **Unified data loading** - Consistent interface for all file formats
5. **Backward compatibility** - Existing code continues to work

## Configuration Structure

Each YAML file can contain one or more data source configurations:

```yaml
source_name:
  # Required fields
  url: "https://example.com/data.tsv"
  format: "tsv"  # tsv, csv, excel, xml, json, rda
  field_mapping:
    field_name: column_index_or_path
  
  # Optional fields
  description: "Human-readable description"
  organism: 9606  # NCBI taxonomy ID
  skip_header: 1
  filters: [...]  # Row filters
  # ... see schema.yaml for all options
```

## File Formats Supported

- **TSV/CSV** - High-performance loading with DuckDB
- **Excel** - XLS/XLSX files via pandas
- **JSON** - Records or nested structures
- **XML** - Automatic flattening with xmltodict
- **RDA** - R data files via pyreadr
- **Parquet** - Direct loading of bronze layer files

## Usage Examples

### Basic Usage

```python
from pypath.inputs.input_module import InputModule

# Load a configured source
module = InputModule('biogrid', 'biogrid_physical')
interactions = module.load_processed()
```

### Creating a New Source

1. Create a YAML file in this directory:

```yaml
# mysource.yaml
my_data:
  url: "https://example.com/data.tsv.gz"
  format: "tsv"
  compression: "gzip"
  field_mapping:
    protein_a: 0
    protein_b: 1
    score: 2
  filters:
    - field: score
      operator: gt
      value: 0.5
  organism: 9606
```

2. Use it in code:

```python
module = InputModule('mysource', 'my_data')
data = module.load_processed()
```

### Migrating Existing Modules

To migrate an existing input module:

1. Create a YAML config with the same parameters
2. Use `create_input_function()` for a drop-in replacement:

```python
from pypath.inputs.input_module import create_input_function

# Creates a function with the old API
my_input_func = create_input_function('mysource', 'my_data')
data = my_input_func(organism=9606)
```

## Change Detection

The system automatically tracks changes using:

- **ETag headers** - HTTP entity tags
- **Last-Modified headers** - Timestamp of last modification
- **Content-Length** - File size changes
- **Checksums** - MD5/SHA256 verification

Configure in YAML:

```yaml
my_source:
  url: "https://example.com/data.tsv"
  check_etag: true  # default
  check_last_modified: true  # default
  checksum_url: "https://example.com/data.tsv.md5"
  checksum_type: "md5"
```

## Bronze Layer

Downloaded data is automatically stored in Parquet format for fast subsequent access:

```python
# First run: downloads and converts to parquet
module = InputModule('large_dataset')
data = module.load_processed()  # May take time

# Subsequent runs: loads from parquet
data = module.load_processed()  # Much faster!
```

## Advanced Features

### Custom Transformations

```yaml
my_source:
  # ... basic config ...
  transform: "mymodule.process_data"
  transform_args:
    normalize: true
```

### Partitioned Datasets

```yaml
large_source:
  # ... basic config ...
  partition_by: ["year", "month"]
  bronze_path: "/custom/path/to/bronze"
```

### Complex Field Mappings

```yaml
json_source:
  format: "json"
  field_mapping:
    id: "data.items[0].id"
    name: "data.items[0].properties.name"
```

## Benefits

1. **Less Code** - Most sources now just need YAML config
2. **Faster Development** - Add new sources without writing Python
3. **Better Performance** - Automatic Parquet caching
4. **Smarter Downloads** - Only fetch when data changes
5. **Consistent Interface** - Same pattern for all sources

## See Also

- `schema.yaml` - Full configuration schema
- `examples/` - Example configurations
- `pypath/inputs/input_module.py` - Implementation details