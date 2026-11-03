"""RaMP mapping backend for metabolite/compound identifiers."""

from __future__ import annotations

import logging

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)


class RampBackend(MappingBackend):
    name = 'ramp'
    yaml_key = 'ramp'

    def read(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        """RaMP mappings are not organism-specific (compound IDs)."""
        try:
            return self._read_via_pypath(id_type, target_id_type, ncbi_tax_id)
        except ImportError:
            _log.debug('pypath not available for ramp backend')
            return {}

    def _read_via_pypath(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        from pypath.inputs.ramp import ramp_mapping, ramp_synonym_mapping

        _log.info('RaMP: %s -> %s', id_type, target_id_type)

        data = {}

        # 1. Try the standard source-to-source mapping
        try:
            source_data = ramp_mapping(id_type, target_id_type)
            if source_data:
                data.update(source_data)
        except Exception as e:
            _log.debug('RaMP source mapping failed: %s', e)

        # 2. Add synonym mappings if applicable
        # ramp_synonym_mapping(id_type) returns {id: set(synonyms)}
        # This is useful for name resolution (synonym -> formal ID)
        try:
            if target_id_type == 'synonym':
                syn_data = ramp_synonym_mapping(id_type)
                if syn_data:
                    for k, v in syn_data.items():
                        data.setdefault(k, set()).update(v)
            elif id_type == 'synonym':
                from pkg_infra.utils import swap_dict

                syn_data = ramp_synonym_mapping(target_id_type)
                if syn_data:
                    swapped = swap_dict(syn_data, force_sets=True)
                    for k, v in swapped.items():
                        data.setdefault(k, set()).update(v)
        except Exception as e:
            _log.debug('RaMP synonym mapping: %s', e)

        if not data:
            return {}

        return {k: (v if isinstance(v, set) else {v}) for k, v in data.items()}

    def _read_direct(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        raise ImportError('ramp requires pypath')


register('ramp', RampBackend)
