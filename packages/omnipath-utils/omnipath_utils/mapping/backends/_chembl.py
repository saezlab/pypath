"""ChEMBL structure-bearing mapping backend (Milestone I) via the inputs_v2 adapter."""

from __future__ import annotations

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._inputs_v2_adapter import InputsV2Backend


class ChemblBackend(InputsV2Backend):
    name = 'chembl'
    yaml_key = 'chembl'
    resource_module = 'chembl'
    dataset = 'molecules'


register('chembl', ChemblBackend)
