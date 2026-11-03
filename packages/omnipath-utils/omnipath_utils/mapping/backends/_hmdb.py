"""HMDB mapping backend for metabolite identifiers."""

from __future__ import annotations

import logging

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)


class HmdbBackend(MappingBackend):
    name = 'hmdb'
    yaml_key = 'hmdb'

    def read(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        """HMDB mappings are human-only."""
        try:
            return self._read_via_pypath(id_type, target_id_type, ncbi_tax_id)
        except ImportError:
            _log.debug('pypath not available for hmdb backend')
            return {}

    def _read_via_pypath(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        from pypath.inputs.hmdb.metabolites import mapping as hmdb_mapping

        _log.info('HMDB: %s -> %s', id_type, target_id_type)

        try:
            data = hmdb_mapping(id_type, target_id_type)
        except Exception as e:
            _log.debug('HMDB: %s', e)
            return {}

        if not data:
            return {}

        return {k: (v if isinstance(v, set) else {v}) for k, v in data.items()}

    def _read_direct(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        raise ImportError('hmdb requires pypath')


register('hmdb', HmdbBackend)
