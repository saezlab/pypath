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
Generic input module base class for simplified data source management.
"""

from __future__ import annotations

import os
import yaml
import importlib
import collections
from pathlib import Path
from typing import Optional, Dict, Any, List, Union, Callable, Iterator

import pandas as pd
import pyarrow as pa

from download_manager import DownloadManager
from download_manager._manager_extended import DownloadManagerExtended

import pypath.share.session as session
import pypath.share.settings as settings
import pypath.bronze.loaders as loaders
import pypath.bronze.storage as storage
import pypath.utils.mapping as mapping

_logger = session.Logger(name='input_module')
_log = _logger._log


class InputModule:
    """
    Generic input module that reads YAML configuration and handles
    downloading, caching, bronze layer storage, and data processing.
    """
    
    def __init__(
        self,
        config: Union[str, Dict[str, Any]],
        source_name: Optional[str] = None,
        use_bronze: bool = True,
        use_cache: bool = True,
        bronze_storage: Optional[storage.BronzeStorage] = None,
        download_manager: Optional[DownloadManager] = None,
    ):
        """
        Initialize input module.
        
        Args:
            config: Path to YAML config file or config dict
            source_name: Name of source in config (if config has multiple)
            use_bronze: Whether to use bronze layer storage
            use_cache: Whether to use download cache
            bronze_storage: Bronze storage instance (creates new if None)
            download_manager: Download manager instance (creates new if None)
        """
        self.config = self._load_config(config, source_name)
        self.source_name = source_name or list(self.config.keys())[0]
        self.source_config = self.config[self.source_name]
        
        self.use_bronze = use_bronze
        self.use_cache = use_cache
        
        # Initialize storage backends
        self.bronze_storage = bronze_storage or storage.BronzeStorage()
        self.download_manager = download_manager or DownloadManagerExtended(
            pkg='pypath',
            config={'backend': 'requests'}
        )
        
        # Cache for loaded data
        self._data_cache = None
        self._arrow_cache = None
    
    def _load_config(
        self,
        config: Union[str, Dict[str, Any]],
        source_name: Optional[str] = None
    ) -> Dict[str, Any]:
        """Load configuration from YAML file or dict."""
        if isinstance(config, str):
            config_path = Path(config)
            if not config_path.is_absolute():
                # Look in resources/sources directory
                config_path = (
                    Path(__file__).parent.parent / 'resources' / 'sources' / config
                )
                if not config_path.suffix:
                    config_path = config_path.with_suffix('.yaml')
            
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
        
        if source_name and source_name in config:
            return {source_name: config[source_name]}
        
        return config
    
    def download(self, force: bool = False) -> Optional[str]:
        """
        Download source data if needed.
        
        Args:
            force: Force download even if unchanged
            
        Returns:
            Path to downloaded file or None
        """
        config = self.source_config
        url = config['url']
        
        # Check for changes if not forcing
        check_config = {
            'check_etag': config.get('check_etag', True),
            'check_last_modified': config.get('check_last_modified', True),
            'checksum_url': config.get('checksum_url'),
            'checksum_type': config.get('checksum_type', 'md5'),
        }
        
        # Download with change detection
        if hasattr(self.download_manager, 'download_if_changed'):
            path, was_downloaded, change_info = self.download_manager.download_if_changed(
                url=url,
                dest=True,  # Use cache
                force=force,
                **check_config,
                **self._get_download_params()
            )
            
            if was_downloaded:
                _log(f'Downloaded {self.source_name} from {url}')
            else:
                _log(f'Using cached {self.source_name}, no changes detected')
            
            return path
        else:
            # Fallback to regular download
            return self.download_manager.download(
                url=url,
                dest=True,
                **self._get_download_params()
            )
    
    def _get_download_params(self) -> Dict[str, Any]:
        """Extract download parameters from config."""
        params = {}
        
        # HTTP method and headers
        if 'method' in self.source_config:
            params['method'] = self.source_config['method']
        if 'headers' in self.source_config:
            params['headers'] = self.source_config['headers']
        if 'params' in self.source_config:
            params['params'] = self.source_config['params']
        if 'data' in self.source_config:
            params['post'] = self.source_config['data']
        if 'json_body' in self.source_config:
            params['json'] = self.source_config['json_body']
        
        return params
    
    def load_raw(self, force_download: bool = False) -> pa.Table:
        """
        Load raw data into Arrow table, using bronze layer if available.
        
        Args:
            force_download: Force new download
            
        Returns:
            PyArrow table with raw data
        """
        if self._arrow_cache is not None and not force_download:
            return self._arrow_cache
        
        # Check bronze layer first
        if self.use_bronze and not force_download:
            bronze_data = self.bronze_storage.load(
                self.source_name,
                self.source_config
            )
            if bronze_data is not None:
                _log(f'Loaded {self.source_name} from bronze layer')
                self._arrow_cache = bronze_data
                return bronze_data
        
        # Download if needed
        file_path = self.download(force=force_download)
        if not file_path:
            raise ValueError(f'Failed to download {self.source_name}')
        
        # Load file based on format
        format_config = self._get_format_config()
        arrow_table = loaders.load_file(file_path, **format_config)
        
        # Save to bronze layer
        if self.use_bronze:
            self.bronze_storage.save(
                arrow_table,
                self.source_name,
                self.source_config,
                partition_by=self.source_config.get('partition_by')
            )
        
        self._arrow_cache = arrow_table
        return arrow_table
    
    def _get_format_config(self) -> Dict[str, Any]:
        """Extract file format configuration."""
        config = self.source_config
        format_config = {
            'format': config.get('format', 'tsv'),
        }
        
        # Add format-specific options
        if 'separator' in config:
            format_config['separator'] = config['separator']
        if 'sheet' in config:
            format_config['sheet'] = config['sheet']
        if 'encoding' in config:
            format_config['encoding'] = config['encoding']
        if 'skip_header' in config:
            format_config['skip_header'] = config['skip_header']
        
        return format_config
    
    def load_processed(
        self,
        force_download: bool = False,
        raw: bool = False
    ) -> Union[List[Any], pd.DataFrame]:
        """
        Load and process data according to configuration.
        
        Args:
            force_download: Force new download
            raw: Return raw DataFrame instead of processed records
            
        Returns:
            Processed data (list of named tuples or DataFrame)
        """
        if self._data_cache is not None and not force_download:
            return self._data_cache
        
        # Load raw data
        arrow_table = self.load_raw(force_download)
        df = arrow_table.to_pandas()
        
        # Apply filters
        df = self._apply_filters(df)
        
        # Apply field mapping
        df = self._apply_field_mapping(df)
        
        # Parse subfields
        df = self._parse_subfields(df)
        
        # Apply custom transform if specified
        if 'transform' in self.source_config:
            df = self._apply_transform(df)
        
        if raw:
            self._data_cache = df
            return df
        
        # Convert to list of named tuples
        records = self._to_records(df)
        self._data_cache = records
        return records
    
    def _apply_filters(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply configured filters to DataFrame."""
        filters = self.source_config.get('filters', [])
        
        for filter_config in filters:
            field = filter_config['field']
            operator = filter_config['operator']
            value = filter_config['value']
            
            if field not in df.columns:
                _log(f'Warning: filter field {field} not in data')
                continue
            
            if operator == 'eq':
                df = df[df[field] == value]
            elif operator == 'ne':
                df = df[df[field] != value]
            elif operator == 'gt':
                df = df[df[field] > value]
            elif operator == 'lt':
                df = df[df[field] < value]
            elif operator == 'gte':
                df = df[df[field] >= value]
            elif operator == 'lte':
                df = df[df[field] <= value]
            elif operator == 'in':
                df = df[df[field].isin(value)]
            elif operator == 'not_in':
                df = df[~df[field].isin(value)]
            elif operator == 'regex':
                df = df[df[field].str.match(value, na=False)]
        
        return df
    
    def _apply_field_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply field mapping configuration."""
        mapping_config = self.source_config.get('field_mapping', {})
        
        # For column index mapping
        if all(isinstance(v, int) for v in mapping_config.values()):
            # Select and rename columns by index
            selected_cols = []
            rename_map = {}
            
            for field_name, col_idx in mapping_config.items():
                if col_idx < len(df.columns):
                    old_name = df.columns[col_idx]
                    selected_cols.append(old_name)
                    rename_map[old_name] = field_name
            
            df = df[selected_cols].rename(columns=rename_map)
        else:
            # For column name mapping or complex mappings
            new_df = pd.DataFrame()
            
            for field_name, field_spec in mapping_config.items():
                if isinstance(field_spec, (int, str)):
                    # Simple column reference
                    if isinstance(field_spec, int):
                        if field_spec < len(df.columns):
                            new_df[field_name] = df.iloc[:, field_spec]
                    else:
                        if field_spec in df.columns:
                            new_df[field_name] = df[field_spec]
                elif isinstance(field_spec, dict):
                    # Complex mapping (e.g., for XML/JSON paths)
                    # This would need more sophisticated handling
                    pass
            
            df = new_df
        
        return df
    
    def _parse_subfields(self, df: pd.DataFrame) -> pd.DataFrame:
        """Parse subfields based on configuration."""
        subfield_config = self.source_config.get('subfield_separator', {})
        
        for field_name, separator in subfield_config.items():
            if field_name in df.columns:
                # Convert to string and split
                df[field_name] = df[field_name].astype(str).str.split(separator)
        
        return df
    
    def _apply_transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply custom transform function."""
        transform_spec = self.source_config['transform']
        transform_args = self.source_config.get('transform_args', {})
        
        # Import the transform function
        if '.' in transform_spec:
            module_name, func_name = transform_spec.rsplit('.', 1)
            module = importlib.import_module(module_name)
            transform_func = getattr(module, func_name)
        else:
            raise ValueError(f'Invalid transform spec: {transform_spec}')
        
        # Apply transform
        return transform_func(df, **transform_args)
    
    def _to_records(self, df: pd.DataFrame) -> List[Any]:
        """Convert DataFrame to list of named tuples."""
        # Create named tuple class
        fields = list(df.columns)
        record_name = f"{self.source_name.title()}Record"
        Record = collections.namedtuple(record_name, fields)
        
        # Convert to records
        records = []
        for _, row in df.iterrows():
            records.append(Record(**row.to_dict()))
        
        return records
    
    def __iter__(self) -> Iterator[Any]:
        """Iterate over processed records."""
        return iter(self.load_processed())
    
    def __len__(self) -> int:
        """Get number of records."""
        return len(self.load_processed())


def create_input_function(
    config_name: str,
    source_name: Optional[str] = None,
    record_processor: Optional[Callable] = None
) -> Callable:
    """
    Factory function to create input module functions with the old API.
    
    Args:
        config_name: Name of YAML config file
        source_name: Source name in config
        record_processor: Optional function to process records
        
    Returns:
        Function compatible with old input module API
    """
    def input_function(
        organism: Optional[int] = None,
        **kwargs
    ) -> List[Any]:
        """Generated input function."""
        # Create input module
        module = InputModule(config_name, source_name)
        
        # Load data
        records = module.load_processed(**kwargs)
        
        # Filter by organism if specified
        if organism is not None:
            organism = str(organism)
            records = [
                r for r in records
                if hasattr(r, 'organism') and str(r.organism) == organism
            ]
        
        # Apply additional processing if provided
        if record_processor:
            records = record_processor(records)
        
        return records
    
    # Set function metadata
    input_function.__name__ = f"{config_name}_{source_name or 'data'}"
    input_function.__doc__ = f"Load data from {config_name}"
    
    return input_function