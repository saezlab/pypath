"""Mapper -- central ID translation manager."""

from __future__ import annotations

import os
import re
import logging
import threading
from collections.abc import Iterable

from omnipath_utils._constants import DEFAULT_ORGANISM
from omnipath_utils.mapping._table import MappingTable, MappingTableKey
from omnipath_utils.mapping._reader import MapReader
from omnipath_utils.mapping._id_types import IdTypeRegistry

_log = logging.getLogger(__name__)

# UniProt AC format
RE_UNIPROT = re.compile(
    r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$'
    r'|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$',
)

# Source ID types where case fallbacks should NOT be applied
_UNIPROT_SOURCE_TYPES = frozenset(
    {
        'uniprot',
        'swissprot',
        'trembl',
        'uniprot-sec',
        'uniprot-pri',
    }
)


class Mapper:
    """Central ID translation manager.

    Manages loaded MappingTables, dispatches loading to MapReader,
    and provides the public translation API.

    In memory mode (default): downloads data from upstream APIs,
    caches locally in pickle files, holds translation dicts in memory.
    """

    _instance: Mapper | None = None
    _lock = threading.Lock()

    def __init__(
        self,
        ncbi_tax_id: int | None = None,
        lifetime: int = 300,
        cachedir: str | None = None,
    ):
        self.ncbi_tax_id = ncbi_tax_id or DEFAULT_ORGANISM
        self.lifetime = lifetime
        self._cachedir = cachedir or self._default_cachedir()
        self.tables: dict[MappingTableKey, MappingTable] = {}
        self._id_types = IdTypeRegistry.get()

    @classmethod
    def get(cls, **kwargs: object) -> Mapper:
        """Singleton access."""

        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls(**kwargs)

        return cls._instance

    @staticmethod
    def _default_cachedir() -> str:
        """Default cache directory for mapping pickles."""

        import platformdirs

        return os.path.join(
            platformdirs.user_cache_dir('omnipath_utils'),
            'mapping',
        )

    def _direct_lookup(
        self,
        name: str,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
    ) -> set[str]:
        """Raw table lookup without any fallback logic."""

        table = self.which_table(id_type, target_id_type, ncbi_tax_id)

        if table:
            return table[name]

        return set()

    def _reverse_lookup(
        self,
        name: str,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
    ) -> set[str]:
        """Reverse lookup: scan target -> source tables for the name."""

        table = self.which_table(
            target_id_type,
            id_type,
            ncbi_tax_id,
        )

        if table:
            result = set()

            for src, targets in table.data.items():
                if name in targets:
                    result.add(src)

            if result:
                return result

        return set()

    def map_name(
        self,
        name: str,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
        strict: bool = False,
        uniprot_cleanup_flag: bool = True,
    ) -> set[str]:
        """Translate a single identifier.

        Args:
            name: The identifier to translate.
            id_type: Source ID type (e.g. genesymbol).
            target_id_type: Target ID type (e.g. uniprot).
            ncbi_tax_id: Organism. Defaults to human (9606).
            strict: If True, skip fuzzy gene symbol fallbacks.
            uniprot_cleanup_flag: If True and target is UniProt,
                run cleanup (secondary -> primary, SwissProt preference).

        Returns:
            Set of target identifiers. Empty set if no mapping found.
        """

        if not name:
            return set()

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        id_type = self._id_types.resolve(id_type) or id_type
        target_id_type = (
            self._id_types.resolve(target_id_type) or target_id_type
        )

        _target_is_uniprot = target_id_type == 'uniprot'

        # Same type -- still do cleanup if target is uniprot
        if id_type == target_id_type:
            if _target_is_uniprot and uniprot_cleanup_flag:
                from omnipath_utils.mapping._cleanup import uniprot_cleanup

                return uniprot_cleanup({name}, ncbi_tax_id, mapper=self)

            return {name}

        # Direct lookup
        result = self._direct_lookup(
            name,
            id_type,
            target_id_type,
            ncbi_tax_id,
        )

        if result:
            if _target_is_uniprot and uniprot_cleanup_flag:
                from omnipath_utils.mapping._cleanup import uniprot_cleanup

                result = uniprot_cleanup(result, ncbi_tax_id, mapper=self)

            return result

        # --- Fallback strategies (no direct result) ---

        # Gene symbol fallbacks -- skip for UniProt/TrEMBL source types
        if (
            id_type in ('genesymbol', 'genesymbol-syn')
            and id_type not in _UNIPROT_SOURCE_TYPES
        ):
            from omnipath_utils.mapping._special import (
                map_genesymbol_fallbacks,
            )

            result = map_genesymbol_fallbacks(
                name,
                target_id_type,
                ncbi_tax_id,
                self,
                strict,
            )

            if result:
                if _target_is_uniprot and uniprot_cleanup_flag:
                    from omnipath_utils.mapping._cleanup import uniprot_cleanup

                    result = uniprot_cleanup(result, ncbi_tax_id, mapper=self)

                return result

        # RefSeq version handling
        if id_type.startswith('refseq'):
            from omnipath_utils.mapping._special import map_refseq

            result = map_refseq(
                name,
                id_type,
                target_id_type,
                ncbi_tax_id,
                self,
                strict,
            )

            if result:
                if _target_is_uniprot and uniprot_cleanup_flag:
                    from omnipath_utils.mapping._cleanup import uniprot_cleanup

                    result = uniprot_cleanup(result, ncbi_tax_id, mapper=self)

                return result

        # Ensembl version stripping
        if id_type.startswith('ens') and '.' in name:
            from omnipath_utils.mapping._special import (
                map_ensembl_strip_version,
            )

            result = map_ensembl_strip_version(
                name,
                id_type,
                target_id_type,
                ncbi_tax_id,
                self,
            )

            if result:
                if _target_is_uniprot and uniprot_cleanup_flag:
                    from omnipath_utils.mapping._cleanup import uniprot_cleanup

                    result = uniprot_cleanup(result, ncbi_tax_id, mapper=self)

                return result

        # miRNA reciprocal fallback
        if id_type.startswith('mir-'):
            from omnipath_utils.mapping._special import map_mirna_fallback

            result = map_mirna_fallback(
                name,
                id_type,
                target_id_type,
                ncbi_tax_id,
                self,
            )

            if result:
                return result

        # CURIE prefix stripping
        if ':' in name:
            from omnipath_utils.mapping._special import strip_prefix

            result = strip_prefix(
                name,
                id_type,
                target_id_type,
                ncbi_tax_id,
                self,
            )

            if result:
                if _target_is_uniprot and uniprot_cleanup_flag:
                    from omnipath_utils.mapping._cleanup import uniprot_cleanup

                    result = uniprot_cleanup(result, ncbi_tax_id, mapper=self)

                return result

        # Chain translation via uniprot
        if id_type != 'uniprot' and target_id_type != 'uniprot':
            from omnipath_utils.mapping._special import chain_map

            result = chain_map(
                name,
                id_type,
                target_id_type,
                ncbi_tax_id,
                self,
            )

            if result:
                return result

        # Try reverse lookup
        result = self._reverse_lookup(
            name,
            id_type,
            target_id_type,
            ncbi_tax_id,
        )

        if result:
            if _target_is_uniprot and uniprot_cleanup_flag:
                from omnipath_utils.mapping._cleanup import uniprot_cleanup

                result = uniprot_cleanup(result, ncbi_tax_id, mapper=self)

            return result

        return set()

    def map_names(
        self,
        names: Iterable[str],
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
    ) -> set[str]:
        """Translate multiple identifiers, return union of results."""

        result = set()

        for name in names:
            result.update(
                self.map_name(name, id_type, target_id_type, ncbi_tax_id),
            )

        return result

    def map_name0(
        self,
        name: str,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
    ) -> str | None:
        """Translate, returning a single result.

        If the mapping is ambiguous, an arbitrary result is picked.
        """

        result = self.map_name(
            name,
            id_type,
            target_id_type,
            ncbi_tax_id,
        )

        return next(iter(result)) if result else None

    def translate(
        self,
        identifiers: Iterable[str],
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
    ) -> dict[str, set[str]]:
        """Batch translate. Returns dict mapping source -> set of targets."""

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id
        id_type = self._id_types.resolve(id_type) or id_type
        target_id_type = (
            self._id_types.resolve(target_id_type) or target_id_type
        )

        table = self.which_table(id_type, target_id_type, ncbi_tax_id)

        if not table:
            return {name: set() for name in identifiers}

        return {name: table[name] for name in identifiers}

    def which_table(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
        load: bool = True,
        backend: str | None = None,
    ) -> MappingTable | None:
        """Find or load the appropriate mapping table.

        Args:
            id_type: Source ID type.
            target_id_type: Target ID type.
            ncbi_tax_id: Organism.
            load: If True, attempt to load from backends when not cached.
            backend: Force a specific backend. None = auto-select.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id
        key = MappingTableKey(id_type, target_id_type, ncbi_tax_id)

        if not backend and key in self.tables:
            return self.tables[key]

        if not load:
            return None

        table = self._load_table(
            id_type,
            target_id_type,
            ncbi_tax_id,
            backend=backend,
        )

        if table:
            self.tables[key] = table

        return table

    def _load_table(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        backend: str | None = None,
    ) -> MappingTable | None:
        """Try loading a table from available backends.

        Args:
            id_type: Source ID type.
            target_id_type: Target ID type.
            ncbi_tax_id: Organism.
            backend: Force a specific backend. None = try all.
        """

        if backend:
            backends = [backend]
        else:
            backends = self._find_backends(id_type, target_id_type)

        for backend_name in backends:
            try:
                reader = MapReader(
                    id_type=id_type,
                    target_id_type=target_id_type,
                    ncbi_tax_id=ncbi_tax_id,
                    backend=backend_name,
                    cachedir=self._cachedir,
                )
                data = reader.load()

                if data:
                    return MappingTable(
                        data=data,
                        id_type=id_type,
                        target_id_type=target_id_type,
                        ncbi_tax_id=ncbi_tax_id,
                        lifetime=self.lifetime,
                    )
            except Exception as e:
                _log.warning(
                    'Failed to load %s -> %s from %s: %s',
                    id_type,
                    target_id_type,
                    backend_name,
                    e,
                )

        return None

    # Map from backend registry name to id_types.yaml backend key
    _BACKEND_YAML_KEYS: dict[str, str] = {
        'uniprot': 'uniprot',
        'uniprot_ftp': 'uniprot',
        'biomart': 'ensembl',
    }

    # Backends that override read() and do not use yaml_key column lookup
    _CUSTOM_BACKENDS: tuple[str, ...] = (
        'mirbase',
        'unichem',
        'ramp',
        'hmdb',
    )

    def _find_backends(
        self,
        id_type: str,
        target_id_type: str,
    ) -> list[str]:
        """Find backends that can provide this ID type pair."""

        backends = []

        for backend_name, yaml_key in self._BACKEND_YAML_KEYS.items():
            src_col = self._id_types.backend_column(
                id_type,
                yaml_key,
            )
            tgt_col = self._id_types.backend_column(
                target_id_type,
                yaml_key,
            )

            if src_col and tgt_col:
                backends.append(backend_name)

        # Custom backends that handle their own ID type support check
        backends.extend(self._CUSTOM_BACKENDS)

        return backends

    def translation_table(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
    ) -> dict[str, set[str]]:
        """Get the full translation table as a dict."""

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id
        id_type = self._id_types.resolve(id_type) or id_type
        target_id_type = (
            self._id_types.resolve(target_id_type) or target_id_type
        )

        table = self.which_table(id_type, target_id_type, ncbi_tax_id)

        return dict(table.data) if table else {}

    def remove_expired(self) -> None:
        """Remove tables that have not been used recently."""

        expired = [k for k, t in self.tables.items() if t.expired]

        for k in expired:
            _log.debug('Removing expired table: %s', k)
            del self.tables[k]

    def id_types(self) -> list[str]:
        """List all known ID types."""

        return self._id_types.all_names()

    def __repr__(self) -> str:
        return (
            f'<Mapper tables={len(self.tables)}, organism={self.ncbi_tax_id}>'
        )
