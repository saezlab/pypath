-- Each relation key's copies as one relation: the copies agree on the statement (ends,
-- predicate, direction, sign, category, qualifiers); their sources and resources are
-- united, their evidence counts summed, and the reference keys taken from the first
-- copy (by resource and id). Most keys have one copy, which is the relation as it is;
-- only the others are grouped (list aggregates hold their groups in memory).
WITH copies AS (
    SELECT relation_key, count(*) AS n FROM work.relation GROUP BY relation_key
), single AS (
    SELECT r.relation_key, r.subject_entity_key, r.object_entity_key, r.statement_kind,
        r.subject_reference_entity_key, r.object_reference_entity_key, r.subject_type, r.predicate,
        r.object_type, r.taxon, r.is_directed, r.sign, r.category, r.interaction_class,
        r.object_aspect_qualifier, r.object_direction_qualifier, r.causal_mechanism_qualifier,
        list_sort(list_distinct(coalesce(r.sources, []))) AS sources, [r._resource] AS resources,
        coalesce(r.evidence_count, 0)::BIGINT AS evidence_count
    FROM work.relation r JOIN copies c USING (relation_key) WHERE c.n = 1
), several AS (
    SELECT r.relation_key, r.subject_entity_key, r.object_entity_key,
        any_value(r.statement_kind) AS statement_kind,
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
    FROM work.relation r JOIN copies c USING (relation_key) WHERE c.n > 1
    GROUP BY r.relation_key, r.subject_entity_key, r.object_entity_key
), merged AS (
    SELECT * FROM single UNION ALL BY NAME SELECT * FROM several
)
SELECT n.relation_n, m.relation_key, m.statement_kind, s.entity_n AS subject_n, o.entity_n AS object_n,
    m.subject_reference_entity_key, m.object_reference_entity_key, m.subject_type, m.predicate,
    m.object_type, m.taxon, m.is_directed, m.sign, m.category, m.interaction_class,
    m.object_aspect_qualifier, m.object_direction_qualifier, m.causal_mechanism_qualifier,
    m.sources, m.resources, m.evidence_count
FROM merged m
JOIN work.relation_numbers n USING (relation_key)
LEFT JOIN work.entity_numbers s ON s.entity_key = m.subject_entity_key
LEFT JOIN work.entity_numbers o ON o.entity_key = m.object_entity_key
