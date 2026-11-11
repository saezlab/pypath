-- Each relation key's copies as one relation: the copies agree on the statement (ends,
-- predicate, direction, sign, category, qualifiers); their sources and resources are
-- united, their evidence counts summed, and the reference keys taken from the first
-- copy (by resource and id).
SELECT n.relation_n, r.relation_key,
    any_value(r.statement_kind) AS statement_kind,
    s.entity_n AS subject_n,
    o.entity_n AS object_n,
    arg_min(r.subject_reference_entity_key, struct_pack(a := r._resource, b := r._resource_version, c := r.relation_id))
        AS subject_reference_entity_key,
    arg_min(r.object_reference_entity_key, struct_pack(a := r._resource, b := r._resource_version, c := r.relation_id))
        AS object_reference_entity_key,
    any_value(r.subject_type) AS subject_type,
    any_value(r.predicate) AS predicate,
    any_value(r.object_type) AS object_type,
    any_value(r.taxon) AS taxon,
    any_value(r.is_directed) AS is_directed,
    any_value(r.sign) AS sign,
    any_value(r.category) AS category,
    any_value(r.interaction_class) AS interaction_class,
    any_value(r.object_aspect_qualifier) AS object_aspect_qualifier,
    any_value(r.object_direction_qualifier) AS object_direction_qualifier,
    any_value(r.causal_mechanism_qualifier) AS causal_mechanism_qualifier,
    list_sort(list_distinct(flatten(list(coalesce(r.sources, []))))) AS sources,
    list(DISTINCT r._resource ORDER BY r._resource) AS resources,
    sum(coalesce(r.evidence_count, 0))::BIGINT AS evidence_count
FROM work.relation r
JOIN work.relation_numbers n USING (relation_key)
LEFT JOIN work.entity_numbers s ON s.entity_key = r.subject_entity_key
LEFT JOIN work.entity_numbers o ON o.entity_key = r.object_entity_key
WHERE {{partition r.relation_key}}
GROUP BY n.relation_n, r.relation_key, s.entity_n, o.entity_n, r.subject_entity_key, r.object_entity_key
