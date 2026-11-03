"""ChEBI structure-bearing mapping backend (Milestone I) via the inputs_v2 adapter.

ChEBI is the highest-leverage chemical namespace (chebi_id + inchikey/inchi/smiles
+ many xrefs) and the one forced onto the adapter (no legacy chebi module).
"""

from __future__ import annotations

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._inputs_v2_adapter import InputsV2Backend


class ChebiBackend(InputsV2Backend):
    name = 'chebi'
    yaml_key = 'chebi'
    resource_module = 'chebi'
    dataset = 'molecules'


register('chebi', ChebiBackend)
