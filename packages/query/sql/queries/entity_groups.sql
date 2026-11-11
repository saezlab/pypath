-- query: entity_groups
-- The grouped explorer's page: the matching entities grouped, chemicals by structure
-- (the connectivity block their copies agree on), other entities by a shared Entrez
-- gene reference; an entity whose copies disagree, or that has no group, alone.
-- Most relations first, then most members, then the group key. Each group as a card
-- (its representative entity), with its first members. group_key: one group only.
-- param: query VARCHAR = ""
-- param: strategy VARCHAR = "chemical_connectivity"
-- param: sources VARCHAR[] = NULL
-- param: entity_types VARCHAR[] = NULL
-- param: taxonomy_ids VARCHAR[] = NULL
-- param: reference_entity_keys VARCHAR[] = NULL
-- param: entity_ids VARCHAR[] = NULL
-- param: limit INTEGER = 20
-- param: cursor VARCHAR = NULL
-- param: group_key VARCHAR = NULL
-- param: member_cursor VARCHAR = NULL
-- param: member_limit INTEGER = 5
-- param: include_details BOOLEAN = false
-- param: include_member_keys BOOLEAN = false
-- param: group_filters JSON = NULL
-- param: group_resources VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"query": "TP53", "strategy": "auto"}

-- step: checked AS value
SELECT CASE WHEN {{strategy}} IN ('auto', 'chemical_connectivity', 'gene_reference') THEN true
    ELSE error('invalid: Unknown grouping strategy: ' || {{strategy}}) END

-- step: q AS value
SELECT trim(coalesce({{query}}, ''))

-- step: ql AS value
SELECT lower({{q}})

-- step: q_upper AS value
SELECT CASE WHEN {{ql}} = '' THEN NULL
    ELSE left({{ql}}, length({{ql}}) - 1) || chr(unicode(right({{ql}}, 1)) + 1) END

-- step: curie_ns AS value
SELECT CASE WHEN regexp_full_match({{q}}, '([A-Za-z][\w.-]*)[:|](.+)')
    THEN lower(regexp_extract({{q}}, '^([A-Za-z][\w.-]*)[:|](.+)$', 1)) END

-- step: curie_id AS value
SELECT CASE WHEN regexp_full_match({{q}}, '([A-Za-z][\w.-]*)[:|](.+)')
    THEN regexp_extract({{q}}, '^([A-Za-z][\w.-]*)[:|](.+)$', 2) END

-- step: phases AS value
SELECT CASE WHEN {{q}} = '' THEN ['all'] WHEN regexp_full_match({{q}}, '[0-9a-fA-F]{64}') THEN ['key']
    ELSE ['prefix', 'contains'] END

-- step: type_names AS value
SELECT (SELECT list(name) FROM op_biolink('class', {{entity_types}}))

-- step: scope AS value
-- the selected resources, else the source filter
SELECT CASE
    WHEN op_unset({{resources}}) AND op_unset({{sources}}) THEN NULL
    WHEN op_unset({{sources}}) THEN {{resources}}
    WHEN op_unset({{resources}}) THEN {{sources}}
    ELSE coalesce(nullif(list_intersect({{resources}}, {{sources}}), []), [''])
END

-- step: cursor_ok AS value
-- a page cursor: [relations, member count, group key] of the last group shown
SELECT CASE WHEN {{cursor}} IS NULL OR {{group_key}} IS NOT NULL THEN false
    WHEN coalesce(try(json_array_length({{cursor}}::JSON)) = 3
        AND json_type({{cursor}}::JSON, '$[0]') IN ('UBIGINT', 'BIGINT') AND ({{cursor}}::JSON->>'$[0]')::BIGINT >= 0
        AND json_type({{cursor}}::JSON, '$[1]') IN ('UBIGINT', 'BIGINT') AND ({{cursor}}::JSON->>'$[1]')::BIGINT >= 1
        AND json_type({{cursor}}::JSON, '$[2]') = 'VARCHAR', false)
    THEN true
    ELSE error('invalid: Invalid group cursor') END

-- step: cursor_relations AS value
SELECT CASE WHEN {{cursor_ok}} THEN ({{cursor}}::JSON->>'$[0]')::BIGINT END

-- step: cursor_members AS value
SELECT CASE WHEN {{cursor_ok}} THEN ({{cursor}}::JSON->>'$[1]')::BIGINT END

-- step: cursor_key AS value
SELECT CASE WHEN {{cursor_ok}} THEN {{cursor}}::JSON->>'$[2]' END

-- step: key_tokens AS list
SELECT unnest(coalesce({{entity_ids}}, []))

-- step: key_scope AS value
SELECT {{resources}}

-- steps: entity_keys

-- step: group_members_rows AS table when {{group_key}} IS NOT NULL
-- one group: the entities listed under it
SELECT DISTINCT _eid FROM entity_group
WHERE {{group_key}} IS NOT NULL
    AND group_key = substr({{group_key}}, strpos({{group_key}}, ':') + 1)
    AND kind = CASE split_part({{group_key}}, ':', 1) WHEN 'connectivity' THEN 'connectivity' WHEN 'gene' THEN 'reference' END
    AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))

-- step: group_member_ids AS list
SELECT _eid FROM __group_members_rows ORDER BY _eid

-- step: candidate_keys AS value
SELECT CASE
    WHEN {{group_key}} IS NULL THEN NULL
    WHEN starts_with({{group_key}}, 'entity:') THEN [substr({{group_key}}, 8)]
    WHEN split_part({{group_key}}, ':', 1) NOT IN ('connectivity', 'gene') THEN []::VARCHAR[]
    ELSE coalesce((SELECT list(DISTINCT entity_key ORDER BY entity_key) FROM entity
        WHERE _eid IN ({{group_member_ids*}})), []::VARCHAR[])
END

-- step: match_keys AS value
SELECT CASE
    WHEN op_unset({{entity_ids}}) THEN {{candidate_keys}}
    WHEN {{candidate_keys}} IS NULL THEN coalesce({{entity_keys}}, [])
    ELSE list_intersect(coalesce({{entity_keys}}, []), {{candidate_keys}})
END

-- step: match_scope AS value
SELECT {{scope}}

-- steps: entity_matches

-- step: keys AS table
-- one row per entity key: the group its copies agree on, or none
WITH copies AS (
    SELECT entity_key, list_contains((SELECT list(entity_type) FROM op_chemical_types), entity_type) AS chemical,
        group_connectivity, reference_entity_key, relation_count
    FROM {{include:entity_match_rows m}}
)
SELECT entity_key,
    CASE WHEN {{strategy}} IN ('auto', 'chemical_connectivity') AND bool_or(coalesce(chemical, false))
        AND min(hash(group_connectivity)) FILTER (WHERE group_connectivity IS NOT NULL)
            = max(hash(group_connectivity)) FILTER (WHERE group_connectivity IS NOT NULL)
        THEN any_value(group_connectivity) FILTER (WHERE group_connectivity IS NOT NULL) END AS connectivity,
    CASE WHEN {{strategy}} IN ('auto', 'gene_reference') AND NOT bool_or(coalesce(chemical, false))
        AND bool_and(coalesce(starts_with(reference_entity_key, 'entrez:'), false))
        AND min(hash(reference_entity_key)) FILTER (WHERE reference_entity_key IS NOT NULL)
            = max(hash(reference_entity_key)) FILTER (WHERE reference_entity_key IS NOT NULL)
        THEN any_value(reference_entity_key) FILTER (WHERE reference_entity_key IS NOT NULL) END AS reference_entity_key,
    sum(coalesce(relation_count, 0))::BIGINT AS relations
FROM copies GROUP BY entity_key

-- step: page AS table
WITH counts AS (
    SELECT 'connectivity:' || connectivity AS group_key, connectivity, NULL::VARCHAR AS reference_entity_key,
        count(*) AS member_count, sum(relations) AS relations
    FROM __keys WHERE connectivity IS NOT NULL GROUP BY connectivity
    UNION ALL
    SELECT 'gene:' || reference_entity_key, NULL, reference_entity_key, count(*), sum(relations)
    FROM __keys WHERE reference_entity_key IS NOT NULL GROUP BY reference_entity_key
    UNION ALL
    SELECT 'entity:' || entity_key, NULL, NULL, 1, relations
    FROM __keys WHERE connectivity IS NULL AND reference_entity_key IS NULL
)
SELECT *, row_number() OVER (ORDER BY relations DESC, member_count DESC, group_key) AS position
FROM counts
WHERE CASE
    WHEN {{group_key}} IS NOT NULL THEN group_key = {{group_key}}
    WHEN NOT {{cursor_ok}} THEN true
    ELSE relations < {{cursor_relations}} OR (relations = {{cursor_relations}}
        AND (member_count < {{cursor_members}} OR (member_count = {{cursor_members}}
            AND group_key > {{cursor_key}})))
END
QUALIFY position <= CASE WHEN {{group_key}} IS NOT NULL THEN 1 ELSE {{limit}} + 1 END

-- step: page_groups AS table
SELECT * FROM __page WHERE position <= CASE WHEN {{group_key}} IS NOT NULL THEN 1 ELSE {{limit}} END

-- step: members AS table
SELECT p.group_key, k.entity_key FROM __page_groups p JOIN __keys k
    ON (p.connectivity IS NOT NULL AND k.connectivity = p.connectivity)
    OR (p.reference_entity_key IS NOT NULL AND k.reference_entity_key = p.reference_entity_key)
    OR (p.group_key = 'entity:' || k.entity_key)

-- step: member_keys AS list
SELECT DISTINCT entity_key FROM __members ORDER BY entity_key

-- step: listed_ids AS list
-- the copies a group's entity_group rows list (in scope), for the groups on the page
SELECT DISTINCT g._eid FROM entity_group g
JOIN __page_groups p ON (p.connectivity IS NOT NULL AND g.kind = 'connectivity' AND g.group_key = p.connectivity)
    OR (p.reference_entity_key IS NOT NULL AND g.kind = 'reference' AND g.group_key = p.reference_entity_key)
WHERE g.group_key IN (SELECT coalesce(connectivity, reference_entity_key) FROM __page_groups)
    AND (op_unset({{scope}}) OR list_contains({{scope}}, g._resource))
ORDER BY 1

-- step: singleton_keys AS list
SELECT DISTINCT m.entity_key FROM __members m JOIN __page_groups p USING (group_key)
WHERE p.connectivity IS NULL AND p.reference_entity_key IS NULL ORDER BY 1

-- step: member_copies AS table when {{include_details}}
-- the members' copies in scope, for structure group cards with details (one entity
-- merged from all their members' copies): those their groups list, and every copy of a
-- single entity (not only those the search matched)
SELECT * EXCLUDE (key_rank, _resource_version), _resource || '/' || _resource_version AS resource
FROM (
    SELECT e.* FROM entity e WHERE _eid IN ({{listed_ids*}})
    UNION ALL BY NAME
    SELECT * FROM entity_by_key
    WHERE entity_key IN ({{singleton_keys*}}) AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource)))
WHERE list_contains({{member_keys}}, entity_key)

-- step: copy_ids AS list when {{include_details}}
SELECT _eid FROM __member_copies ORDER BY _eid

-- step: identifier_rows AS table
SELECT * FROM entity_identifier WHERE _eid IN ({{copy_ids*}})

-- step: annotation_rows AS table
SELECT * FROM entity_annotation WHERE _eid IN ({{copy_ids*}})

-- step: member_merged AS table
-- the members as merged at build (when a copy is in scope)
SELECT * FROM entity_merged
WHERE entity_key IN ({{member_keys*}}) AND (op_unset({{scope}}) OR list_has_any(resources, {{scope}}))

-- step: member_summaries AS table
-- the stored summaries, or with details built with the identifier and annotation lists
WITH placed AS (SELECT *, 0 AS position FROM __member_merged WHERE {{include_details}}),
listed AS (SELECT * FROM {{include:merged_rows placed}}),
summaries AS (
    SELECT op_summary(summary) AS summary FROM __member_merged WHERE NOT {{include_details}}
    UNION ALL
    SELECT summary FROM {{include:entity_summary listed}}
)
SELECT s.summary, m.entity_key, m.relation_count, m.entity_type, m.taxon, m.label
FROM summaries s JOIN __member_merged m ON m.entity_key = s.summary.entityPk

-- step: chemical_bases AS table when {{include_details}}
-- with details, a structure group's card merges all its members' copies (as one entity
-- of the first member's key); else it is its first member
WITH chosen AS (
    SELECT p.group_key, c.* REPLACE ((SELECT min(entity_key) FROM __members WHERE group_key = p.group_key) AS entity_key)
    FROM __page_groups p JOIN __members mb USING (group_key) JOIN __member_copies c USING (entity_key)
    WHERE p.connectivity IS NOT NULL
), copy_identifiers AS (
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
    -- each group's chosen copies as one entity (the group key stands for the entity key)
    SELECT c.* EXCLUDE (group_key), c.group_key AS card_group, 0 AS position,
        CASE WHEN {{include_details}} THEN coalesce(i.identifiers, []) END AS identifiers,
        CASE WHEN {{include_details}} THEN coalesce(a.annotations, []) END AS annotations
    FROM chosen c LEFT JOIN copy_identifiers i USING (_eid) LEFT JOIN copy_annotations a USING (_eid)
), merged AS (
    SELECT * FROM {{include:merged_entity nested}}
), shaped AS (
    SELECT * FROM {{include:entity_summary merged}}
)
SELECT g.card_group AS group_key, s.summary FROM shaped s
JOIN (SELECT DISTINCT entity_key, card_group FROM nested) g ON g.entity_key = s.summary.entityPk

-- step: groups AS table
WITH member_list AS (
    SELECT mb.group_key, ms.*, row_number() OVER (PARTITION BY mb.group_key ORDER BY ms.entity_key) AS n
    FROM __members mb JOIN __member_summaries ms USING (entity_key)
), visible AS (
    SELECT *, row_number() OVER (PARTITION BY group_key ORDER BY entity_key) AS v
    FROM member_list WHERE entity_key > coalesce({{member_cursor}}, '')
), stats AS (
    SELECT group_key,
        list(summary ORDER BY entity_key) AS summaries,
        sum(coalesce(relation_count, 0))::BIGINT AS relations,
        count(*) AS member_count,
        count(DISTINCT coalesce(taxon, chr(1))) AS taxa,
        -- the label of the most connected well-named member, else the preferred label
        min(struct_pack(r := -coalesce(relation_count, 0), k := name_rank, l := summary.label))
            FILTER (WHERE coalesce(summary.label, '') <> '' AND name_rank.score < 8) AS best_named,
        min(name_rank) FILTER (WHERE trim(coalesce(summary.label, '')) <> '') AS best_label,
        list(DISTINCT summary.entityType ORDER BY summary.entityType) AS member_types,
        min(entity_key) FILTER (WHERE summary.entityType = 'gene') AS first_gene,
        min(entity_key) AS first_key
    -- the rank per row: DuckDB 2.0 cannot copy the macro's subquery inside an aggregate
    FROM (SELECT *, op_name_rank(summary.label) AS name_rank FROM member_list)
    GROUP BY group_key
), group_sources AS (
    SELECT group_key, list(DISTINCT s ORDER BY s) AS sources
    FROM (SELECT group_key, unnest(summary.sources) AS s FROM member_list)
    GROUP BY group_key
), unions AS (
    -- the members' identifiers by (type, id) and attributes by value, first seen first
    -- each identity at its first place, with its last value (as a dict keeps them)
    SELECT group_key, list(item ORDER BY first) AS identifiers
    FROM (
        SELECT group_key, min(struct_pack(n := n, p := ip)) AS first, arg_max(i, struct_pack(n := n, p := ip)) AS item
        FROM (SELECT group_key, n, unnest(summary.identifiers) AS i, generate_subscripts(summary.identifiers, 1) AS ip FROM member_list)
        GROUP BY group_key, i.identifierType, i.identifier)
    GROUP BY group_key
), attributes AS (
    SELECT group_key, list(item ORDER BY first) AS attributes
    FROM (
        SELECT group_key, min(struct_pack(n := n, p := ap)) AS first, arg_max(a, struct_pack(n := n, p := ap)) AS item
        FROM (SELECT group_key, n, unnest(coalesce(summary.entityAttributes, [])) AS a,
            generate_subscripts(coalesce(summary.entityAttributes, []), 1) AS ap FROM member_list)
        GROUP BY group_key, a::VARCHAR)
    GROUP BY group_key
)
SELECT p.*, st.summaries, st.relations AS member_relations, st.member_count AS shown_members, st.taxa,
    coalesce(gs.sources, []) AS sources, st.best_named, st.best_label, st.member_types, st.first_gene, st.first_key,
    coalesce(u.identifiers, []) AS identifiers, a.attributes,
    (SELECT list(summary ORDER BY v) FROM visible WHERE visible.group_key = p.group_key AND v <= {{member_limit}}) AS shown,
    (SELECT max(entity_key) FROM visible WHERE visible.group_key = p.group_key AND v = {{member_limit}}) AS member_cursor_at,
    (SELECT count(*) FROM visible WHERE visible.group_key = p.group_key) AS visible_count,
    (SELECT list(entity_key ORDER BY v) FROM visible WHERE visible.group_key = p.group_key AND v <= {{member_limit}}) AS shown_keys,
    (SELECT list(entity_key ORDER BY entity_key) FROM member_list WHERE member_list.group_key = p.group_key) AS all_keys,
    coalesce(cb.summary, (SELECT summary FROM __member_summaries WHERE entity_key = st.first_key)) AS chemical_base,
    (SELECT summary FROM __member_summaries WHERE entity_key = coalesce(st.first_gene, st.first_key)) AS gene_base
FROM __page_groups p
LEFT JOIN stats st USING (group_key) LEFT JOIN group_sources gs USING (group_key)
LEFT JOIN unions u USING (group_key) LEFT JOIN attributes a USING (group_key)
LEFT JOIN __chemical_bases cb USING (group_key)

-- step: response AS value
WITH cards AS (
    SELECT position, group_key, connectivity, reference_entity_key, member_count, relations, shown,
        connectivity IS NOT NULL OR reference_entity_key IS NOT NULL AS is_group,
        CASE WHEN visible_count > {{member_limit}} THEN member_cursor_at END AS next_member_cursor,
        CASE
        WHEN connectivity IS NOT NULL THEN op_json_patch(op_json_patch(op_json_patch(to_json(chemical_base), to_json(struct_pack(
            identifiers := identifiers, identifiersTotal := len(identifiers), sources := sources,
            displayName := coalesce(best_named.l, best_label.name, chemical_base.displayName),
            canonicalIdentifier := connectivity, canonicalIdentifierType := 'connectivity',
            primaryNamespace := 'connectivity', primaryIdentifier := connectivity,
            groupStrategy := 'chemical_connectivity', ontologyHierarchy := NULL::JSON))),
            CASE WHEN {{include_details}} THEN to_json(struct_pack(entityAttributes := CASE WHEN len(coalesce(attributes, [])) > 0 THEN attributes END)) ELSE '{}'::JSON END),
            CASE WHEN taxa > 1 THEN to_json(struct_pack(taxonomyId := NULL::VARCHAR, taxonomyName := NULL::VARCHAR, taxon := NULL::VARCHAR)) ELSE '{}'::JSON END)
        WHEN reference_entity_key IS NOT NULL THEN op_json_patch(to_json(gene_base), to_json(struct_pack(
            displayName := CASE WHEN first_gene IS NOT NULL THEN gene_base.displayName
                ELSE coalesce(best_named.l, best_label.name, gene_base.displayName) END,
            entityAttributes := CASE WHEN len(coalesce(attributes, [])) > 0 THEN attributes END,
            identifiers := identifiers, identifiersTotal := len(identifiers), sources := sources,
            referenceEntityKey := reference_entity_key, entityType := 'gene', memberEntityTypes := member_types,
            canonicalIdentifier := split_part(reference_entity_key, ':', 2), canonicalIdentifierType := 'entrez',
            groupStrategy := 'gene_reference')))
        ELSE to_json(summaries[1])
        END AS card,
        shown_keys, all_keys, shown_members, member_relations
    FROM __groups
), entities AS (
    SELECT *, CASE WHEN is_group THEN op_json_patch(card, to_json(struct_pack(
            relationCount := member_relations, entityPk := group_key,
            groupMemberKeys := CASE WHEN {{include_member_keys}} THEN all_keys ELSE coalesce(shown_keys, []) END,
            groupMemberCount := shown_members, groupDetailsLoaded := {{include_details}},
            groupQuery := {{query}}, groupFilters := {{group_filters}},
            groupResources := coalesce({{group_resources}}, []::VARCHAR[]))))
        ELSE card END AS entity_json
    FROM cards
)
SELECT to_json(struct_pack(
    groups := coalesce((SELECT list(to_json(struct_pack(
            group_key := group_key, connectivity := connectivity, reference_entity_key := reference_entity_key,
            member_count := member_count, relations := relations,
            entity := op_json_patch(entity_json, to_json(struct_pack(groupMemberCursor := next_member_cursor))),
            members := CASE WHEN is_group THEN coalesce(shown, [])
                ELSE [op_json_patch(entity_json, to_json(struct_pack(groupMemberCursor := next_member_cursor)))] END,
            is_group := is_group,
            group_label := coalesce(connectivity, reference_entity_key),
            nextMemberCursor := next_member_cursor)) ORDER BY position) FROM entities), []),
    nextCursor := CASE WHEN {{group_key}} IS NULL AND (SELECT count(*) FROM __page) > {{limit}}
        THEN (SELECT '[' || relations || ', ' || member_count || ', ' || to_json(group_key) || ']'
              FROM __page WHERE position = {{limit}}) END))
