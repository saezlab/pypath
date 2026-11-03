"""Small interaction vocabulary expressed as native Biolink qualified statements.

These profiles describe source assertions, not activity thresholds or predictions.
Binding assay magnitudes/comparators remain observations, not a binary activity call.
"""

BINDING_QUALIFIERS = (("causal_mechanism_qualifier", "binding"),)
TRANSPORT_QUALIFIERS = (
    ("object_aspect_qualifier", "transport"),
    ("causal_mechanism_qualifier", "relocalization"),
)


EFFECT_QUALIFIERS = (
    "object_direction_qualifier",
    "object_aspect_qualifier",
    "causal_mechanism_qualifier",
)


def relation_qualifiers(annotations=None) -> dict[str, str]:
    """Effect qualifiers of a relation ("increased", "activity", "binding"), first value each."""
    found: dict[str, str] = {}
    for a in annotations or []:
        term, value = a.get("term"), a.get("value")
        if term in EFFECT_QUALIFIERS and value and a.get("scope", "relation") == "relation":
            found.setdefault(term, str(value))
    return found


def interaction_label(predicate: str, annotations=None) -> str | None:
    """Optional qualified label; never infer a profile from the resource name."""
    values = {
        a["term"]: a.get("value")
        for a in annotations or []
        if a.get("scope", "relation") == "relation"
    }
    if predicate == "interacts_with" and values.get("causal_mechanism_qualifier") == "binding":
        return "Binds"
    if predicate == "affects" and all(values.get(k) == v for k, v in TRANSPORT_QUALIFIERS):
        # A transport inhibitor is still an effect, not a transporter.
        if not values.get("object_direction_qualifier"):
            return "Transports"
    return None
