-- An entity's resource copies merged into one row, as the API merges them
-- (_merge_entity_row_group): the richest copy's columns, relation counts summed,
-- identifiers united by (ns, id) and annotations by value. Reads {{input}}: the copies
-- (entity columns and resource) with nested identifiers and annotations lists; gives
-- the columns the entity_summary fragment reads, and evidence_total.
WITH merged AS (
    -- The richest copy's columns, relation counts summed, identifiers united by (ns, id)
    -- and annotations by value when there are several copies.
    SELECT
        any_value(entity_key) AS entity_key,
        arg_max(struct_pack(namespace := namespace, identifier := identifier, label := label,
            entity_type := entity_type, taxon := taxon, has_hierarchy := has_hierarchy,
            parent_count := parent_count, child_count := child_count), struct_pack(
            a := len(annotations),
            i := len(identifiers),
            g := regexp_full_match(coalesce(label, ''), '[A-Z][A-Z0-9-]{1,14}')::INTEGER,
            n := (coalesce(label, '') <> '' AND coalesce(label, '') <> coalesce(identifier, '')
                AND NOT regexp_full_match(coalesce(label, ''), '[0-9]+'))::INTEGER,
            f := lower(coalesce(label, '')),
            l := coalesce(label, ''),
            k := entity_key,
            r := resource,
            h := coalesce(has_hierarchy, false),
            p := coalesce(parent_count, 0),
            c := coalesce(child_count, 0))) AS best,
        sum(coalesce(relation_count, 0))::BIGINT AS relation_count,
        sum(coalesce(evidence_count, 0))::BIGINT AS evidence_total,
        list_sort(list_distinct(flatten(list(coalesce(gene_reference_keys, []))))) AS gene_reference_keys,
        list(DISTINCT reference_entity_key)
            FILTER (WHERE coalesce(reference_entity_key, '') <> '') AS reference_keys,
        flatten(list(identifiers ORDER BY resource, entity_id)) AS all_identifiers,
        flatten(list(annotations ORDER BY resource, entity_id)) AS all_annotations,
        count(*) AS copies
    FROM {{input}}
    HAVING count(*) > 0
), main_entity AS (
    SELECT
        'main' AS item, 0 AS position,
        entity_key AS pk, best.namespace AS raw_namespace, best.identifier AS raw_identifier,
        best.label AS raw_label, best.entity_type AS raw_type, best.taxon AS taxon,
        best.has_hierarchy AS has_hierarchy, best.parent_count AS parent_count,
        best.child_count AS child_count, relation_count, gene_reference_keys,
        CASE
            WHEN list_contains(coalesce(reference_keys, []),
                    coalesce(best.namespace, 'None') || ':' || coalesce(best.identifier, 'None'))
                OR len(coalesce(reference_keys, [])) > 1
            THEN coalesce(best.namespace, 'None') || ':' || coalesce(best.identifier, 'None')
            WHEN len(coalesce(reference_keys, [])) = 1 THEN reference_keys[1]
        END AS reference_entity_key,
        list_transform(list_filter(
            list_transform(all_identifiers, lambda x, i: struct_pack(x := x, i := i)),
            lambda s: list_position(list_transform(all_identifiers,
                lambda x: struct_pack(n := coalesce(x.ns, chr(1)), d := coalesce(x.id, chr(1)))),
                struct_pack(n := coalesce(s.x.ns, chr(1)), d := coalesce(s.x.id, chr(1)))) = s.i),
            lambda s: s.x) AS identifiers,
        CASE WHEN copies > 1 THEN list_transform(list_filter(
                list_transform(all_annotations, lambda x, i: struct_pack(x := x, i := i)),
                lambda s: list_position(all_annotations, s.x) = s.i), lambda s: s.x)
            ELSE all_annotations END AS annotations,
        [entity_key] AS source_keys,
        evidence_total
    FROM merged
)
SELECT * FROM main_entity
