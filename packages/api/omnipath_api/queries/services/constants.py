from biolink_model.datamodel.model import slots
from omnipath_core.biolink import slot_name

# Public filter fields are the canonical Biolink qualifier slot names.
RELATION_QUALIFIER_FILTERS = tuple(
    map(
        slot_name,
        (
            slots.object_aspect_qualifier,
            slots.object_direction_qualifier,
            slots.causal_mechanism_qualifier,
        ),
    )
)
