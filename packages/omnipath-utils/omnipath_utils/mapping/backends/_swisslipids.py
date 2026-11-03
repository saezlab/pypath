"""SwissLipids structure-bearing mapping backend (Milestone I) via the inputs_v2 adapter.

SwissLipids provides SMILES + a LIPID MAPS cross-reference (no InChIKey of its
own); it broadens SMILES coverage and bridges to lipidmaps' InChIKeys.
"""

from __future__ import annotations

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._inputs_v2_adapter import InputsV2Backend


class SwissLipidsBackend(InputsV2Backend):
    name = 'swisslipids'
    yaml_key = 'swisslipids'
    resource_module = 'swisslipids'
    dataset = 'lipids'

    def _read_via_pypath(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        *,
        src_col: str,
        tgt_col: str,
        limit: int | None = None,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        """SwissLipids' raw ``InChI key (pH7.3)`` column carries the literal
        ``InChIKey=`` format prefix (unlike the id_types.yaml schema's own
        EntityBuilder transform, which strips it -- the inputs_v2 adapter
        reads raw rows directly, bypassing that transform). Stored verbatim,
        this produces a malformed structure key: 37 characters instead of 27,
        never matching another namespace's bare InChIKey on join (T085/T086).
        """
        from omnipath_utils.mapping._id_types import normalize_identifier

        mapping = super()._read_via_pypath(
            id_type, target_id_type, ncbi_tax_id,
            src_col=src_col, tgt_col=tgt_col, limit=limit, **kwargs,
        )
        if target_id_type != 'inchikey':
            return mapping
        cleaned: dict[str, set[str]] = {}
        for source, targets in mapping.items():
            keys = {
                normalized
                for value in targets
                if (normalized := normalize_identifier('inchikey', value))
            }
            if keys:
                cleaned[source] = keys
        return cleaned


register('swisslipids', SwissLipidsBackend)
