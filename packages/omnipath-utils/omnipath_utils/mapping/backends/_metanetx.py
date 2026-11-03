"""MetaNetX/MNXref mapping backend for metabolite cross-references.

MetaNetX provides a reconciled namespace for metabolic networks.  The
MNXref ``chem_xref.tsv`` file maps external database identifiers (ChEBI,
HMDB, KEGG, BiGG, LipidMaps, etc.) to MetaNetX compound IDs (MNXM*).

This backend delegates to ``pypath.inputs.metanetx.metanetx_mapping()``
which reads the cross-reference file and builds pairwise mappings between
any two supported metabolite ID types using MetaNetX IDs as a bridge.
"""

from __future__ import annotations

import logging

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)


class MetaNetXBackend(MappingBackend):
    name = "metanetx"
    yaml_key = "metanetx"

    def read(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        """MetaNetX mappings are organism-agnostic."""

        try:
            return self._read_via_pypath(
                id_type, target_id_type, ncbi_tax_id,
            )
        except ImportError:
            _log.debug("pypath not available for metanetx backend")
            return {}

    def _read_via_pypath(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        from pypath.inputs.metanetx import metanetx_mapping

        _log.info("MetaNetX: %s -> %s", id_type, target_id_type)

        try:
            data = metanetx_mapping(id_type, target_id_type)
        except Exception as e:
            _log.warning("MetaNetX: %s", e)
            return {}

        if not data:
            return {}

        # metanetx_mapping returns dict[str, set[str]] — correct format
        return data

    def _read_direct(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        raise ImportError("metanetx requires pypath")


register("metanetx", MetaNetXBackend)
