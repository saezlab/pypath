"""BiGG mapping backend for metabolite identifiers.

Reads BiGG universal metabolite cross-references from the BiGG Models
bulk TSV file (via ``pypath.inputs.bigg.bigg_metabolite_mapping()``).
Covers all 9,090 universal metabolites across 85+ BiGG models.
"""

from __future__ import annotations

import logging

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)


class BiggBackend(MappingBackend):
    name = "bigg"
    yaml_key = "bigg"

    def read(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        """BiGG mappings are organism-agnostic."""

        try:
            return self._read_via_pypath(
                id_type, target_id_type, ncbi_tax_id,
            )
        except ImportError:
            _log.debug("pypath not available for bigg backend")
            return {}

    def _read_via_pypath(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        from pypath.inputs.bigg import bigg_metabolite_mapping

        _log.info("BiGG: %s -> %s", id_type, target_id_type)

        # Determine direction: bigg->X or X->bigg
        if id_type == "bigg":
            target = target_id_type
        elif target_id_type == "bigg":
            target = id_type
        else:
            _log.debug("BiGG backend requires bigg as source or target")
            return {}

        try:
            data = bigg_metabolite_mapping(target)
        except Exception as e:
            _log.warning("BiGG: %s", e)
            return {}

        if not data:
            return {}

        # Reverse if target is bigg (invert the mapping)
        if target_id_type == "bigg":
            from collections import defaultdict

            reversed_data: dict[str, set[str]] = defaultdict(set)

            for bigg_id, target_ids in data.items():
                for tid in target_ids:
                    reversed_data[tid].add(bigg_id)

            return dict(reversed_data)

        return data

    def _read_direct(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        raise ImportError("bigg requires pypath")


register("bigg", BiggBackend)
