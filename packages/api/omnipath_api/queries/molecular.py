"""Molecular form filters over relation evidence; no resolution at read time."""


def has_form_filters(filters):
    return any(
        filters.get(k)
        for k in ("protein_entity_keys", "transcript_entity_keys", "isoform_identifiers")
    )


def form_match_sql(filters):
    """All specified identity constraints match one endpoint of one evidence row
    (``relation_evidence``); returns the condition and its parameters."""
    params = {"subject": [], "object": []}
    sides = []
    for side in ("subject", "object"):
        form = f"{side}_molecular_form"
        clauses = []
        for key, value in [
            ("protein_entity_keys", f"{form}.protein_entity_key"),
            ("transcript_entity_keys", f"{form}.transcript_entity_key"),
            (
                "isoform_identifiers",
                f"{form}.isoform_identifier.ns || ':' || {form}.isoform_identifier.id",
            ),
        ]:
            if filters.get(key):
                clauses.append(f"list_contains(?::VARCHAR[], {value})")
                params[side].append(filters[key])
        sides.append("(" + (" AND ".join(clauses) or "TRUE") + ")")
    mode = filters.get("molecular_endpoint_mode", "any")
    if mode == "source":
        return sides[0], params["subject"]
    if mode == "target":
        return sides[1], params["object"]
    joined = f"({sides[0]} {'AND' if mode == 'both' else 'OR'} {sides[1]})"
    return joined, params["subject"] + params["object"]


def matching_evidence(items, filters):
    def endpoint(form):
        form = form or {}
        for key, field in [
            ("protein_entity_keys", "protein_entity_key"),
            ("transcript_entity_keys", "transcript_entity_key"),
        ]:
            if filters.get(key) and form.get(field) not in filters[key]:
                return False
        iso = form.get("isoform_identifier") or {}
        return (
            not filters.get("isoform_identifiers")
            or f"{iso.get('ns')}:{iso.get('id')}" in filters["isoform_identifiers"]
        )

    result = []
    for item in items or []:
        subject, obj = (
            endpoint(item.get("subject_molecular_form")),
            endpoint(item.get("object_molecular_form")),
        )
        mode = filters.get("molecular_endpoint_mode", "any")
        match = (
            subject
            if mode == "source"
            else obj
            if mode == "target"
            else subject and obj
            if mode == "both"
            else subject or obj
        )
        if match:
            result.append(item)
    return result
