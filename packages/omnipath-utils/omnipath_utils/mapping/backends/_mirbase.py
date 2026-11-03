"""miRBase mapping backend."""

from __future__ import annotations

import logging
from collections import defaultdict

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)


class MirBaseBackend(MappingBackend):
    name = 'mirbase'
    yaml_key = 'mirbase'

    # miRBase doesn't use the standard yaml_key backend column lookup.
    # It provides its own mapping tables. Override read() directly.

    def read(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        """Load miRBase ID mappings.

        Supports mappings between:
        - mir-pre (precursor name) <-> mirbase (miRBase accession MI*)
        - mir-mat-name (mature name) <-> mirbase
        - mir-pre <-> mir-mat-name
        """
        try:
            return self._read_via_pypath(id_type, target_id_type, ncbi_tax_id)
        except ImportError:
            _log.debug('pypath not available for mirbase backend')
            return {}

    def _read_via_pypath(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        from pypath.inputs import mirbase as mirbase_input
        from pkg_infra.utils import swap_dict

        # Map our canonical type names to mirbase function names
        # mirbase functions yield (key, value) tuples:
        # mirbase_precursor: (accession, name) => mirbase -> mir-pre
        # mirbase_mature: (accession, name) => mirbase -> mir-mat-name
        # mirbase_precursor_to_mature: (pre_name, mat_name) => mir-pre -> mir-mat-name
        _TABLES = {
            ('mirbase', 'mir-pre'): lambda org: mirbase_input.mirbase_precursor(
                org
            ),
            ('mirbase', 'mir-mat-name'): lambda org: (
                mirbase_input.mirbase_mature(org)
            ),
            ('mir-pre', 'mir-mat-name'): lambda org: (
                mirbase_input.mirbase_precursor_to_mature(org)
            ),
        }

        key = (id_type, target_id_type)
        rev_key = (target_id_type, id_type)

        if key in _TABLES:
            raw = _TABLES[key](ncbi_tax_id)
            return self._to_dict(raw)
        elif rev_key in _TABLES:
            raw = _TABLES[rev_key](ncbi_tax_id)
            return swap_dict(self._to_dict(raw), force_sets=True)

        _log.debug('mirbase: no table for %s -> %s', id_type, target_id_type)
        return {}

    def _read_direct(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        raise ImportError('mirbase requires pypath')

    @staticmethod
    def _to_dict(data) -> dict[str, set[str]]:
        """Convert mirbase data (generator of tuples) to our format."""
        result = defaultdict(set)
        if isinstance(data, dict):
            for k, v in data.items():
                if isinstance(v, (set, list, tuple)):
                    result[k].update(str(i) for i in v)
                else:
                    result[k].add(str(v))
        else:
            # Generator of (key, value) tuples
            for k, v in data:
                result[str(k)].add(str(v))
        return dict(result)


register('mirbase', MirBaseBackend)
