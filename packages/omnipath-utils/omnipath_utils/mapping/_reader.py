"""MapReader -- loads mapping data from backends and manages caching."""

from __future__ import annotations

import os
import pickle
import hashlib
import logging

_log = logging.getLogger(__name__)


class MapReader:
    """Reads ID mapping data from a backend and returns raw data.

    Handles pickle caching: each unique combination of parameters
    produces a deterministic cache filename via MD5 hash.
    """

    def __init__(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        backend: str,
        cachedir: str | None = None,
        **backend_params: object,
    ):
        self.id_type = id_type
        self.target_id_type = target_id_type
        self.ncbi_tax_id = ncbi_tax_id
        self.backend = backend
        self.backend_params = backend_params
        self._cachedir = cachedir

    @property
    def cache_path(self) -> str | None:
        """Path for the pickle cache file, or None if no cachedir."""

        if not self._cachedir:
            return None

        os.makedirs(self._cachedir, exist_ok=True)

        key_str = (
            f'{self.id_type}_{self.target_id_type}'
            f'_{self.ncbi_tax_id}_{self.backend}'
        )

        if self.backend_params:
            key_str += '_' + str(sorted(self.backend_params.items()))

        md5 = hashlib.md5(key_str.encode()).hexdigest()[:12]

        return os.path.join(
            self._cachedir,
            f'mapping_{self.id_type}__{self.target_id_type}'
            f'__{self.ncbi_tax_id}__{md5}.pickle',
        )

    def load(self) -> dict[str, set[str]]:
        """Load mapping data, using cache if available."""

        cache = self.cache_path

        if cache and os.path.exists(cache):
            _log.debug('Loading mapping from cache: %s', cache)

            try:
                with open(cache, 'rb') as f:
                    return pickle.load(f)
            except Exception:
                _log.warning('Failed to load cache %s, reloading', cache)

        _log.info(
            'Loading mapping %s -> %s (organism %d) from %s',
            self.id_type,
            self.target_id_type,
            self.ncbi_tax_id,
            self.backend,
        )
        data = self._read_from_backend()

        if cache and data:
            try:
                with open(cache, 'wb') as f:
                    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
                _log.debug('Saved mapping to cache: %s', cache)
            except Exception:
                _log.warning('Failed to save cache %s', cache)

        return data

    def _read_from_backend(self) -> dict[str, set[str]]:
        """Dispatch to the appropriate backend reader."""

        from omnipath_utils.mapping.backends import get_backend

        backend_reader = get_backend(self.backend)

        if not backend_reader:
            _log.warning('Unknown backend: %s', self.backend)
            return {}

        return backend_reader.read(
            id_type=self.id_type,
            target_id_type=self.target_id_type,
            ncbi_tax_id=self.ncbi_tax_id,
            **self.backend_params,
        )
