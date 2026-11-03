"""COSMOS prior-knowledge network, projected from the reaction stars.

One build-time table, ``cosmos_edge``, holding the binary metabolite/enzyme
edges the COSMOS R package consumes: a reactant feeds the enzyme that turns it
over, that enzyme produces the products, and a connector ties each bare
identifier to the per-reaction gene node that carries it.

The step reads the build database and nothing else. The interactions phase has
already assembled every reaction into an N-ary header with its participants in
role, so the projection is a read of that layer rather than a second pass over
the resources.

A reaction the resources call reversible is emitted twice, the second half with
the arrow turned round and a `_rev` gene node of its own. The direction comes
from what the resources published about the event and is never defaulted: a
reaction nobody stated a direction for keeps a null, because directionless in
the source is a different claim from known to run one way.

The labels are put into the namespaces the formalism reads — UniProt on the
gene side, ChEBI on the metabolite side — by a translation stage that pushes
the identifiers down to the mapping database and brings the answers back. What
it cannot translate keeps the identifier the build canonicalised it to, and
every endpoint records the namespace its label actually carries, so a consumer
tells a translated label from a fallen-back one from the data rather than by
parsing a string. A build with no reachable mapping database falls every label
back and is complete either way.
"""

from __future__ import annotations

from omnipath_build.cosmos.build import (
    TABLE,
    CosmosBuildStats,
    build_cosmos_projection,
    build_id,
    ensure_cosmos_edge_table,
)
from omnipath_build.cosmos.translate import (
    UTILS_SCHEMA,
    TranslationStats,
    stage_label_translation,
)

__all__ = [
    'TABLE',
    'UTILS_SCHEMA',
    'CosmosBuildStats',
    'TranslationStats',
    'build_cosmos_projection',
    'build_id',
    'ensure_cosmos_edge_table',
    'stage_label_translation',
]
