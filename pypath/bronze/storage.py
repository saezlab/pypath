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
Bronze layer storage management for parquet files.
"""

from __future__ import annotations

import os
import json
import hashlib
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, Union, List

import pyarrow as pa
import pyarrow.parquet as pq
import platformdirs

import pypath.share.session as session
import pypath.share.settings as settings

_logger = session.Logger(name='bronze_storage')
_log = _logger._log


class BronzeStorage:
    """
    Manages bronze layer parquet file storage with metadata.
    """
    
    def __init__(
        self,
        base_path: Optional[Union[str, Path]] = None,
        pkg: str = 'pypath',
        use_pypath_data: bool = True,
    ):
        """
        Initialize bronze storage.
        
        Args:
            base_path: Base directory for bronze storage
            pkg: Package name for default directory
            use_pypath_data: If True and base_path is None, use pypath/data directory
        """
        if base_path is None:
            if use_pypath_data:
                # Use pypath/data directory
                import pypath
                pypath_dir = Path(pypath.__file__).parent
                base_path = pypath_dir / 'data'
            else:
                # Use platform-specific data directory
                base_path = platformdirs.user_data_dir(pkg, 'omnipathdb')
        
        self.base_path = Path(base_path) / 'bronze'
        self.base_path.mkdir(parents=True, exist_ok=True)
        
        self.metadata_file = self.base_path / 'metadata.json'
        self.metadata = self._load_metadata()
    
    def _load_metadata(self) -> Dict[str, Any]:
        """Load metadata from JSON file."""
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        return {}
    
    def _save_metadata(self):
        """Save metadata to JSON file."""
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2, default=str)
    
    def _generate_key(self, source_name: str, config: Dict[str, Any]) -> str:
        """
        Generate unique key for bronze file based on source and config.
        
        Args:
            source_name: Name of the data source
            config: Source configuration
            
        Returns:
            Unique key string
        """
        # Create a stable hash of relevant config parameters
        key_parts = [
            source_name,
            config.get('url', ''),
            config.get('version', ''),
            str(config.get('organism', '')),
        ]
        
        # Include filter values in the key to differentiate cached data
        filters = config.get('filters', [])
        for f in sorted(filters, key=lambda x: x.get('field', '')):
            if f.get('field') in ('organism', 'taxid', 'taxid_a', 'taxid_b'):
                key_parts.append(f'{f.get("field")}={f.get("value")}')
        
        key_str = '|'.join(key_parts)
        key_hash = hashlib.md5(key_str.encode()).hexdigest()[:8]
        
        return f"{source_name}_{key_hash}"
    
    def get_path(
        self,
        source_name: str,
        config: Dict[str, Any],
        partition_by: Optional[List[str]] = None,
    ) -> Path:
        """
        Get path for bronze parquet file.
        
        Args:
            source_name: Name of the data source
            config: Source configuration
            partition_by: Fields to partition by
            
        Returns:
            Path to parquet file/directory
        """
        # Check if custom path is specified
        if 'bronze_path' in config:
            return Path(config['bronze_path'])
        
        key = self._generate_key(source_name, config)
        
        if partition_by:
            # Partitioned dataset directory
            return self.base_path / f"{key}_partitioned"
        else:
            # Single parquet file
            return self.base_path / f"{key}.parquet"
    
    def save(
        self,
        table: pa.Table,
        source_name: str,
        config: Dict[str, Any],
        source_metadata: Optional[Dict[str, Any]] = None,
        partition_by: Optional[List[str]] = None,
        **kwargs
    ) -> Path:
        """
        Save Arrow table to bronze layer.
        
        Args:
            table: PyArrow table to save
            source_name: Name of the data source
            config: Source configuration
            source_metadata: Additional metadata to store
            partition_by: Fields to partition by
            **kwargs: Additional arguments for parquet writer
            
        Returns:
            Path where data was saved
        """
        path = self.get_path(source_name, config, partition_by)
        partition_by = partition_by or config.get('partition_by')
        
        _log(f'Saving bronze data to: {path}')
        
        # Prepare metadata
        key = self._generate_key(source_name, config)
        metadata = {
            'source_name': source_name,
            'url': config.get('url'),
            'format': config.get('format'),
            'saved_at': datetime.now().isoformat(),
            'num_rows': table.num_rows,
            'num_columns': table.num_columns,
            'schema': table.schema.to_string(),
            'config': config,
        }
        
        if source_metadata:
            metadata['source_metadata'] = source_metadata
        
        # Save parquet file
        if partition_by:
            pq.write_to_dataset(
                table,
                root_path=path,
                partition_cols=partition_by,
                **kwargs
            )
        else:
            pq.write_table(table, path, **kwargs)
        
        # Update metadata
        self.metadata[key] = metadata
        self._save_metadata()
        
        _log(f'Successfully saved {table.num_rows} rows to bronze layer')
        
        return path
    
    def load(
        self,
        source_name: str,
        config: Dict[str, Any],
        columns: Optional[List[str]] = None,
        filters: Optional[List[tuple]] = None,
        **kwargs
    ) -> Optional[pa.Table]:
        """
        Load data from bronze layer.
        
        Args:
            source_name: Name of the data source
            config: Source configuration
            columns: Columns to load (None for all)
            filters: PyArrow filters for partitioned datasets
            **kwargs: Additional arguments for parquet reader
            
        Returns:
            PyArrow table or None if not found
        """
        path = self.get_path(source_name, config)
        
        if not path.exists():
            _log(f'Bronze file not found: {path}')
            return None
        
        _log(f'Loading bronze data from: {path}')
        
        try:
            if path.is_dir():
                # Partitioned dataset
                table = pq.read_table(
                    path,
                    columns=columns,
                    filters=filters,
                    **kwargs
                )
            else:
                # Single file
                table = pq.read_table(
                    path,
                    columns=columns,
                    **kwargs
                )
            
            _log(f'Successfully loaded {table.num_rows} rows from bronze layer')
            return table
            
        except Exception as e:
            _logger._log_traceback()
            _log(f'Failed to load bronze data: {e}')
            return None
    
    def exists(self, source_name: str, config: Dict[str, Any]) -> bool:
        """
        Check if bronze data exists.
        
        Args:
            source_name: Name of the data source
            config: Source configuration
            
        Returns:
            True if data exists
        """
        path = self.get_path(source_name, config)
        return path.exists()
    
    def get_metadata(
        self,
        source_name: str,
        config: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """
        Get metadata for bronze data.
        
        Args:
            source_name: Name of the data source
            config: Source configuration
            
        Returns:
            Metadata dict or None if not found
        """
        key = self._generate_key(source_name, config)
        return self.metadata.get(key)
    
    def delete(self, source_name: str, config: Dict[str, Any]) -> bool:
        """
        Delete bronze data.
        
        Args:
            source_name: Name of the data source
            config: Source configuration
            
        Returns:
            True if deleted successfully
        """
        path = self.get_path(source_name, config)
        key = self._generate_key(source_name, config)
        
        if path.exists():
            if path.is_dir():
                import shutil
                shutil.rmtree(path)
            else:
                path.unlink()
            
            # Remove metadata
            if key in self.metadata:
                del self.metadata[key]
                self._save_metadata()
            
            _log(f'Deleted bronze data: {path}')
            return True
        
        return False
    
    def list_sources(self) -> List[Dict[str, Any]]:
        """
        List all sources in bronze storage.
        
        Returns:
            List of metadata dicts
        """
        return list(self.metadata.values())
    
    def clean_old(self, days: int = 30) -> int:
        """
        Clean bronze files older than specified days.
        
        Args:
            days: Delete files older than this many days
            
        Returns:
            Number of files deleted
        """
        from datetime import timedelta
        
        cutoff = datetime.now() - timedelta(days=days)
        deleted = 0
        
        for key, meta in list(self.metadata.items()):
            saved_at = datetime.fromisoformat(meta['saved_at'])
            if saved_at < cutoff:
                config = meta.get('config', {})
                source_name = meta.get('source_name')
                if source_name and self.delete(source_name, config):
                    deleted += 1
        
        _log(f'Cleaned {deleted} old bronze files')
        return deleted