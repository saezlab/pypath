-- The relation ends that no resource has an entity row for, one row per resource stating
-- them: what the relation itself says of the end (type, label, reference key, taxon,
-- from the resource's first such relation), so that every end is an entity.
WITH ends AS (
    SELECT subject_entity_key AS entity_key, subject_type AS entity_type, subject_label AS label,
        subject_reference_entity_key AS reference_entity_key, taxon, _resource, _resource_version, relation_id
    FROM work.relation
    UNION ALL
    SELECT object_entity_key, object_type, object_label, object_reference_entity_key, taxon,
        _resource, _resource_version, relation_id
    FROM work.relation
)
SELECT entity_key, _resource, any_value(_resource_version) AS _resource_version,
    arg_min(struct_pack(entity_type := entity_type, label := label, reference_entity_key := reference_entity_key,
        taxon := taxon), relation_id) AS stub
FROM ends
WHERE entity_key IS NOT NULL AND entity_key NOT IN (SELECT entity_key FROM work.entity WHERE entity_key IS NOT NULL)
GROUP BY entity_key, _resource
