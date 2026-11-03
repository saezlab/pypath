"""LIPID MAPS structure-bearing mapping backend (Milestone I) via the inputs_v2 adapter."""

from __future__ import annotations

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._inputs_v2_adapter import InputsV2Backend


class LipidMapsBackend(InputsV2Backend):
    name = 'lipidmaps'
    yaml_key = 'lipidmaps'
    resource_module = 'lipidmaps'
    dataset = 'lipids'


register('lipidmaps', LipidMapsBackend)
