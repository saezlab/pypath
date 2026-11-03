"""KEGG Compound → ChEBI mapping backend for omnipath-utils.

Merges two independent sources:

1. MetaNetX ``chem_xref.tsv`` kegg→chebi column (~45 K pairs), loaded via
   ``pypath.inputs.metanetx.metanetx_mapping("kegg", "chebi")``.

2. KEGG REST ``/conv/chebi/compound`` endpoint (~17 K pairs), loaded via
   ``pypath.inputs.kegg_api._kegg_conv``.

Both sources are merged by union per KEGG ID.  Any KEGG ID that maps to
**more than one distinct ChEBI** after merging is excluded (1→many abort);
the count of excluded pairs is logged at ``INFO`` level.

ChEBI IDs are normalised to the ``CHEBI:<number>`` format throughout.

Wire-up in ``omnipath_utils/data/id_types.yaml``:
  kegg.backends.kegg_compound = kegg
  chebi.backends.kegg_compound = chebi
"""

from __future__ import annotations

import logging
from collections import defaultdict

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)


class KeggCompoundBackend(MappingBackend):
    name = "kegg_compound"
    yaml_key = "kegg_compound"

    def read(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        """KEGG compound mappings are organism-agnostic."""
        try:
            return self._read_via_pypath(id_type, target_id_type, ncbi_tax_id)
        except ImportError:
            _log.debug("pypath not available for kegg_compound backend")
            return {}

    def _read_via_pypath(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        from pypath.inputs.metanetx import metanetx_mapping

        try:
            from pypath.inputs.kegg_api import _kegg_conv as _rest_fn
        except ImportError:
            _log.warning(
                "KeggCompound: pypath.inputs.kegg_api unavailable "                "(pycurl missing?); KEGG REST source disabled"
            )
            _rest_fn = lambda *a, **kw: {}

        _log.info("KeggCompound: %s -> %s", id_type, target_id_type)

        if id_type == "kegg" and target_id_type == "chebi":
            return self._load_kegg_to_chebi(metanetx_mapping, _rest_fn)
        elif id_type == "chebi" and target_id_type == "kegg":
            fwd = self._load_kegg_to_chebi(metanetx_mapping, _rest_fn)
            return self._reverse(fwd)
        else:
            _log.debug(
                "KeggCompound backend only supports kegg<->chebi; got %s->%s",
                id_type,
                target_id_type,
            )
            return {}

    def _load_kegg_to_chebi(
        self,
        metanetx_mapping,
        _kegg_conv,
    ) -> dict[str, set[str]]:
        merged: dict[str, set[str]] = defaultdict(set)

        # Source 1: MetaNetX chem_xref.tsv kegg→chebi (~45 K pairs)
        try:
            for kegg_id, chebi_set in metanetx_mapping("kegg", "chebi").items():
                merged[kegg_id].update(chebi_set)
        except Exception as exc:
            _log.warning("KeggCompound: MetaNetX source failed: %s", exc)

        # Source 2: KEGG REST /conv/chebi/compound (~17 K pairs)
        # After source_split / target_split the pairs are e.g. C00022 → 15361
        try:
            rest = _kegg_conv(
                "compound", "chebi",
                source_split=True,
                target_split=True,
            )
            for kegg_id, raw_set in rest.items():
                for raw_id in raw_set:
                    chebi_id = (
                        raw_id
                        if raw_id.upper().startswith("CHEBI:")
                        else f"CHEBI:{raw_id}"
                    )
                    merged[kegg_id].add(chebi_id)
        except Exception as exc:
            _log.warning("KeggCompound: KEGG REST source failed: %s", exc)

        # 1→many abort
        result: dict[str, set[str]] = {}
        excluded = 0
        for kegg_id, chebi_set in merged.items():
            if len(chebi_set) > 1:
                excluded += 1
            else:
                result[kegg_id] = chebi_set

        _log.info(
            "KeggCompound: %d kegg→chebi pairs loaded; %d excluded (1→many)",
            len(result),
            excluded,
        )
        return result

    @staticmethod
    def _reverse(fwd: dict[str, set[str]]) -> dict[str, set[str]]:
        rev: dict[str, set[str]] = defaultdict(set)
        for kegg_id, chebi_set in fwd.items():
            for chebi_id in chebi_set:
                rev[chebi_id].add(kegg_id)
        return dict(rev)

    def _read_direct(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        raise ImportError("kegg_compound requires pypath")


register("kegg_compound", KeggCompoundBackend)
