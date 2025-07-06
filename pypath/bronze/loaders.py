#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
Data loaders for various file formats to Apache Arrow tables.
"""

from __future__ import annotations

import os
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, Dict, Any, Union, List

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb
import xmltodict

try:
    import pyreadr
except ImportError:
    pyreadr = None

import pypath.share.session as session

_logger = session.Logger(name='bronze_loaders')
_log = _logger._log

# Initialize DuckDB connection
DUCK = duckdb.connect(':memory:')


def load_tsv(
    path: Union[str, Path],
    separator: str = '\t',
    skip_header: int = 0,
    encoding: str = 'utf-8',
    column_names: Optional[List[str]] = None,
    **kwargs
) -> pa.Table:
    """
    Load TSV file using DuckDB for performance.
    
    Args:
        path: Path to TSV file
        separator: Field separator (default: tab)
        skip_header: Number of header rows to skip
        encoding: File encoding
        column_names: Optional column names to use
        **kwargs: Additional arguments passed to DuckDB's read_csv
        
    Returns:
        PyArrow Table
    """
    _log(f'Loading TSV file: {path}')
    
    # Build DuckDB query
    read_csv_args = {
        'delim': separator,
        'skip': skip_header,
        'header': column_names is None,
        'auto_detect': True,
    }
    
    if column_names:
        read_csv_args['columns'] = {f'column{i}': dtype 
                                    for i, dtype in enumerate(['VARCHAR'] * len(column_names))}
        read_csv_args['names'] = column_names
    
    # Add any additional kwargs
    read_csv_args.update(kwargs)
    
    # Format arguments for SQL
    args_str = ', '.join(f"{k}={repr(v)}" for k, v in read_csv_args.items())
    
    try:
        query = f"SELECT * FROM read_csv_auto('{path}', {args_str})"
        result = DUCK.execute(query).arrow()
        _log(f'Successfully loaded {result.num_rows} rows from TSV')
        return result
    except Exception as e:
        _logger._log_traceback()
        _log(f'Failed to load TSV with DuckDB, falling back to pandas: {e}')
        # Fallback to pandas
        df = pd.read_csv(
            path,
            sep=separator,
            skiprows=skip_header,
            encoding=encoding,
            names=column_names,
            header=0 if column_names is None else None,
        )
        return pa.Table.from_pandas(df, preserve_index=False)


def load_csv(
    path: Union[str, Path],
    separator: str = ',',
    skip_header: int = 0,
    encoding: str = 'utf-8',
    column_names: Optional[List[str]] = None,
    **kwargs
) -> pa.Table:
    """
    Load CSV file using DuckDB for performance.
    
    Args:
        path: Path to CSV file
        separator: Field separator (default: comma)
        skip_header: Number of header rows to skip
        encoding: File encoding
        column_names: Optional column names to use
        **kwargs: Additional arguments passed to DuckDB's read_csv
        
    Returns:
        PyArrow Table
    """
    return load_tsv(
        path=path,
        separator=separator,
        skip_header=skip_header,
        encoding=encoding,
        column_names=column_names,
        **kwargs
    )


def load_excel(
    path: Union[str, Path],
    sheet: Union[str, int] = 0,
    skip_header: int = 0,
    column_names: Optional[List[str]] = None,
    **kwargs
) -> pa.Table:
    """
    Load Excel file using pandas and convert to Arrow.
    
    Args:
        path: Path to Excel file
        sheet: Sheet name or index
        skip_header: Number of header rows to skip
        column_names: Optional column names to use
        **kwargs: Additional arguments passed to pandas read_excel
        
    Returns:
        PyArrow Table
    """
    _log(f'Loading Excel file: {path}, sheet: {sheet}')
    
    try:
        df = pd.read_excel(
            path,
            sheet_name=sheet,
            skiprows=skip_header,
            names=column_names,
            header=0 if column_names is None else None,
            engine='openpyxl',
            **kwargs
        )
        
        result = pa.Table.from_pandas(df, preserve_index=False)
        _log(f'Successfully loaded {result.num_rows} rows from Excel')
        return result
    except Exception as e:
        _logger._log_traceback()
        raise ValueError(f'Failed to load Excel file: {e}')


def load_json(
    path: Union[str, Path],
    orient: str = 'records',
    encoding: str = 'utf-8',
    **kwargs
) -> pa.Table:
    """
    Load JSON file and convert to Arrow.
    
    Args:
        path: Path to JSON file
        orient: JSON orientation (records, index, columns, values, table)
        encoding: File encoding
        **kwargs: Additional arguments passed to pandas read_json
        
    Returns:
        PyArrow Table
    """
    _log(f'Loading JSON file: {path}')
    
    try:
        # For simple records orientation, use direct JSON loading
        if orient == 'records':
            with open(path, 'r', encoding=encoding) as f:
                data = json.load(f)
            if isinstance(data, list) and data and isinstance(data[0], dict):
                df = pd.DataFrame(data)
            else:
                df = pd.read_json(path, orient=orient, encoding=encoding, **kwargs)
        else:
            df = pd.read_json(path, orient=orient, encoding=encoding, **kwargs)
        
        result = pa.Table.from_pandas(df, preserve_index=False)
        _log(f'Successfully loaded {result.num_rows} rows from JSON')
        return result
    except Exception as e:
        _logger._log_traceback()
        raise ValueError(f'Failed to load JSON file: {e}')


def load_xml(
    path: Union[str, Path],
    record_path: Optional[str] = None,
    encoding: str = 'utf-8',
    flatten: bool = True,
    **kwargs
) -> pa.Table:
    """
    Load XML file and convert to Arrow.
    
    Args:
        path: Path to XML file
        record_path: Path to records in XML structure
        encoding: File encoding
        flatten: Whether to flatten nested structures
        **kwargs: Additional arguments
        
    Returns:
        PyArrow Table
    """
    _log(f'Loading XML file: {path}')
    
    try:
        with open(path, 'r', encoding=encoding) as f:
            content = f.read()
        
        if flatten:
            # Use xmltodict for simple flattening
            data = xmltodict.parse(content)
            
            # If record_path is specified, navigate to it
            if record_path:
                for key in record_path.split('.'):
                    data = data.get(key, data)
            
            # Convert to DataFrame
            if isinstance(data, list):
                df = pd.DataFrame(data)
            elif isinstance(data, dict):
                # Try to find the main data list
                for key, value in data.items():
                    if isinstance(value, list):
                        df = pd.DataFrame(value)
                        break
                else:
                    # Single record
                    df = pd.DataFrame([data])
            else:
                df = pd.json_normalize(data)
        else:
            # Use pandas XML reader for more complex cases
            df = pd.read_xml(path, encoding=encoding, **kwargs)
        
        result = pa.Table.from_pandas(df, preserve_index=False)
        _log(f'Successfully loaded {result.num_rows} rows from XML')
        return result
    except Exception as e:
        _logger._log_traceback()
        raise ValueError(f'Failed to load XML file: {e}')


def load_rda(
    path: Union[str, Path],
    which: Union[int, str] = 0,
    **kwargs
) -> pa.Table:
    """
    Load R data file (.rda, .RData) and convert to Arrow.
    
    Args:
        path: Path to RDA file
        which: Which object to load (index or name)
        **kwargs: Additional arguments
        
    Returns:
        PyArrow Table
    """
    if pyreadr is None:
        raise ImportError(
            'pyreadr is required for loading R data files. '
            'Install with: pip install pyreadr'
        )
    
    _log(f'Loading RDA file: {path}')
    
    try:
        # Read R data file
        result = pyreadr.read_r(str(path))
        
        # Get the specified dataframe
        if isinstance(which, int):
            df_name = list(result.keys())[which]
            df = result[df_name]
        else:
            df = result[which]
        
        # Convert to Arrow
        table = pa.Table.from_pandas(df, preserve_index=False)
        _log(f'Successfully loaded {table.num_rows} rows from RDA')
        return table
    except Exception as e:
        _logger._log_traceback()
        raise ValueError(f'Failed to load RDA file: {e}')


def load_parquet(path: Union[str, Path], **kwargs) -> pa.Table:
    """
    Load existing parquet file.
    
    Args:
        path: Path to parquet file
        **kwargs: Additional arguments passed to pq.read_table
        
    Returns:
        PyArrow Table
    """
    _log(f'Loading parquet file: {path}')
    
    try:
        table = pq.read_table(path, **kwargs)
        _log(f'Successfully loaded {table.num_rows} rows from parquet')
        return table
    except Exception as e:
        _logger._log_traceback()
        raise ValueError(f'Failed to load parquet file: {e}')


def detect_format(path: Union[str, Path]) -> str:
    """
    Detect file format from extension.
    
    Args:
        path: File path
        
    Returns:
        Format string (tsv, csv, excel, json, xml, rda, parquet)
    """
    path = Path(path)
    ext = path.suffix.lower()
    
    format_map = {
        '.tsv': 'tsv',
        '.tab': 'tsv',
        '.txt': 'tsv',  # Assume tab-separated for .txt
        '.csv': 'csv',
        '.xls': 'excel',
        '.xlsx': 'excel',
        '.xlsm': 'excel',
        '.json': 'json',
        '.xml': 'xml',
        '.rda': 'rda',
        '.rdata': 'rda',
        '.parquet': 'parquet',
        '.pq': 'parquet',
    }
    
    return format_map.get(ext, 'tsv')  # Default to TSV


def load_file(
    path: Union[str, Path],
    format: Optional[str] = None,
    **kwargs
) -> pa.Table:
    """
    Load any supported file format to Arrow table.
    
    Args:
        path: File path
        format: File format (auto-detected if None)
        **kwargs: Format-specific arguments
        
    Returns:
        PyArrow Table
    """
    path = Path(path)
    
    if format is None:
        format = detect_format(path)
    
    loaders = {
        'tsv': load_tsv,
        'csv': load_csv,
        'excel': load_excel,
        'json': load_json,
        'xml': load_xml,
        'rda': load_rda,
        'parquet': load_parquet,
    }
    
    loader = loaders.get(format)
    if loader is None:
        raise ValueError(f'Unsupported format: {format}')
    
    return loader(path, **kwargs)