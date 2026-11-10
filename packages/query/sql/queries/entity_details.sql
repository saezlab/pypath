-- query: entity_details
-- The explorer's entity details: the entity merged from its resource copies, its first
-- 50 structural relationships with both endpoints, its relation count and up to 100
-- association relations, as the API's JSON. NULL when no selected resource has it.
-- param: entity_key VARCHAR
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_key": "missing"}

-- step: copies AS table
SELECT * EXCLUDE (key_rank, _resource_version), _resource || '/' || _resource_version AS resource
FROM entity_by_key
WHERE entity_key = {{entity_key}} AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: copy_ids AS list
SELECT _eid FROM __copies ORDER BY _eid

-- step: identifier_rows AS table
SELECT * FROM entity_identifier WHERE _eid IN ({{copy_ids*}})

-- step: annotation_rows AS table
SELECT * FROM entity_annotation WHERE _eid IN ({{copy_ids*}})

-- step: main AS table
WITH copy_identifiers AS (
    SELECT _eid, list(struct_pack(ns := ns, id := id, is_canonical := is_canonical, source := source)
        ORDER BY ordinal) AS identifiers
    FROM __identifier_rows GROUP BY ALL
), copy_annotations AS (
    SELECT _eid, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset) ORDER BY ordinal) AS annotations
    FROM __annotation_rows GROUP BY ALL
), nested AS (
    SELECT c.*, 0 AS position, coalesce(i.identifiers, []) AS identifiers, coalesce(a.annotations, []) AS annotations
    FROM __copies c LEFT JOIN copy_identifiers i USING (_eid) LEFT JOIN copy_annotations a USING (_eid)
)
SELECT * FROM {{include:merged_entity nested}}

-- step: endpoint_entries AS table
SELECT DISTINCT _rid, _resource, resource, relation_id, predicate, category
FROM relation_endpoint
WHERE key = {{entity_key}} AND key_kind = 'entity'
    AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: structural_page AS table
SELECT _rid, resource, relation_id, count(*) OVER () AS total,
    row_number() OVER (ORDER BY predicate, resource, relation_id) AS position
FROM (SELECT DISTINCT _rid, resource, relation_id, predicate FROM __endpoint_entries
      WHERE predicate IN ('has_input', 'has_output', 'enabled_by', 'has_member', 'has_part'))
QUALIFY position <= 50

-- step: association_page AS table
SELECT _rid, resource, relation_id
FROM (SELECT DISTINCT _rid, resource, relation_id, predicate FROM __endpoint_entries
      WHERE category = 'association')
ORDER BY predicate, resource, relation_id LIMIT 100

-- step: page_ids AS list
SELECT _rid FROM (SELECT _rid FROM __structural_page UNION SELECT _rid FROM __association_page) ORDER BY _rid

-- step: page_relations AS table
SELECT * FROM relation WHERE _rid IN ({{page_ids*}})

-- step: structural_ids AS list
SELECT _rid FROM __structural_page ORDER BY _rid

-- step: structural_annotation_rows AS table
SELECT * FROM relation_annotation WHERE _rid IN ({{structural_ids*}})

-- step: structural AS table
SELECT p.position, r.* FROM __structural_page p JOIN __page_relations r USING (_rid)

-- step: endpoint_keys AS list
SELECT DISTINCT k FROM (
    SELECT subject_entity_key AS k FROM __structural UNION SELECT object_entity_key FROM __structural)
WHERE k IS NOT NULL ORDER BY k

-- step: endpoint_rows AS table
SELECT * FROM entity_by_key WHERE entity_key IN ({{endpoint_keys*}})

-- step: summaries AS table
WITH endpoints AS (
    SELECT side AS item, s.position,
        CASE side WHEN 'subject' THEN s.subject_entity_key ELSE s.object_entity_key END AS key,
        CASE side WHEN 'subject' THEN s.subject_type ELSE s.object_type END AS side_type,
        CASE side WHEN 'subject' THEN s.subject_label ELSE s.object_label END AS side_label,
        CASE side WHEN 'subject' THEN s.subject_reference_entity_key
            ELSE s.object_reference_entity_key END AS side_reference,
        s._resource, s.taxon AS relation_taxon
    FROM __structural s, (VALUES ('subject'), ('object')) AS sides(side)
), joined AS (
    SELECT p.*, e.entity_key AS e_entity_key, e.namespace AS e_namespace, e.identifier AS e_identifier,
        e.label AS e_label, e.entity_type AS e_entity_type, e.taxon AS e_taxon,
        e.has_hierarchy AS e_has_hierarchy, e.parent_count AS e_parent_count,
        e.child_count AS e_child_count, e.relation_count AS e_relation_count,
        e.gene_reference_keys AS e_gene_reference_keys, e.reference_entity_key AS e_reference_entity_key
    FROM endpoints p LEFT JOIN __endpoint_rows e ON e._resource = p._resource AND e.entity_key = p.key
), to_shape AS (
    SELECT * FROM __main
    UNION ALL BY NAME
    SELECT * FROM {{include:endpoint_entity joined}}
)
SELECT * FROM {{include:entity_summary to_shape}}

-- step: structural_annotations AS table
SELECT _rid, list(struct_pack(term := term, value := value,
    quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
        quantity_has_binary_relation, quantity_source_field, quantity_comparator),
    source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
FROM __structural_annotation_rows GROUP BY ALL

-- step: response AS value when (SELECT count(*) FROM __copies) > 0
WITH relationship_items AS (
    SELECT s.position, struct_pack(
        relation := op_relation_json(s, coalesce(a.annotations, [])),
        subjectEntity := subj.summary,
        objectEntity := obj.summary,
        annotations := coalesce(a.annotations, [])) AS item
    FROM __structural s
    LEFT JOIN __structural_annotations a USING (_rid)
    JOIN __summaries subj ON subj.item = 'subject' AND subj.position = s.position
    JOIN __summaries obj ON obj.item = 'object' AND obj.position = s.position
)
SELECT to_json(struct_pack(
    -- the main entity with its evidence count, merged as JSON: DuckDB 2.0 drops the
    -- name of struct_insert's field inside a macro
    entity := json_merge_patch(to_json((SELECT summary FROM __summaries WHERE item = 'main')),
        to_json(struct_pack(molecularEvidenceTotal := (SELECT evidence_total FROM __main)))),
    relationships := coalesce((SELECT list(item ORDER BY position) FROM relationship_items), []),
    relationshipsTotal := coalesce((SELECT any_value(total) FROM __structural_page), 0),
    summary := struct_pack(interactionCount := (SELECT count(*) FROM __endpoint_entries)),
    annotations := coalesce((SELECT list(struct_pack(relationPk := r.relation_key, predicate := r.predicate)
            ORDER BY r.relation_key, r._resource || '/' || r._resource_version)
        FROM __association_page p JOIN __page_relations r USING (_rid)), [])
))
