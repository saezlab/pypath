-- query: entity_search
-- Entities matching a name or identifier, merged across resources and ranked as in the
-- explorer: an empty query lists every entity, an entity key finds that entity, other
-- text finds prefix matches of labels and identifiers (and a ns:id CURIE that
-- identifier, or the entity it is the reference of), else labels containing it. Prefix
-- matches rank exact labels, exact identifiers, label prefixes, identifier prefixes,
-- then shorter terms, then entities with more relations (over every selected copy).
-- param: query VARCHAR = ""
-- param: sources VARCHAR[] = NULL
-- param: entity_types VARCHAR[] = NULL
-- param: taxonomy_ids VARCHAR[] = NULL
-- param: reference_entity_keys VARCHAR[] = NULL
-- param: entity_ids VARCHAR[] = NULL
-- param: limit INTEGER = 20
-- param: cursor_phase VARCHAR = NULL
-- param: cursor_offset INTEGER = NULL
-- param: cursor_key VARCHAR = NULL
-- param: resources VARCHAR[] = NULL
-- output: table rows
-- output: json response
-- example: {"query": "TP53"}

-- step: q AS value
SELECT trim(coalesce({{query}}, ''))

-- step: ql AS value
SELECT lower({{q}})

-- step: scope AS value
-- the resources read: the selected ones and the source filter's; [''] selects none
SELECT CASE
    WHEN op_unset({{resources}}) AND op_unset({{sources}}) THEN NULL
    WHEN op_unset({{sources}}) THEN {{resources}}
    WHEN op_unset({{resources}}) THEN {{sources}}
    ELSE coalesce(nullif(list_intersect({{resources}}, {{sources}}), []), [''])
END

-- step: query_key AS value
SELECT sha256(to_json(struct_pack(query := {{q}}, sources := {{sources}}, entity_types := {{entity_types}},
    taxonomy_ids := {{taxonomy_ids}}, reference_entity_keys := {{reference_entity_keys}},
    entity_ids := {{entity_ids}}, resources := {{resources}}))::VARCHAR)

-- step: checked AS value
SELECT CASE WHEN {{cursor_key}} IS NOT NULL AND {{cursor_key}} <> {{query_key}}
    THEN error('invalid: Entity cursor belongs to another query') ELSE true END

-- step: offset_ AS value
SELECT coalesce({{cursor_offset}}, 0)

-- step: count_ AS value
SELECT {{offset_}} + {{limit}} + 1

-- step: phases AS value
SELECT CASE
    WHEN {{cursor_phase}} IS NOT NULL THEN [{{cursor_phase}}]
    WHEN {{q}} = '' THEN ['all']
    WHEN regexp_full_match({{q}}, '[0-9a-fA-F]{64}') THEN ['key']
    ELSE ['prefix', 'contains']
END

-- step: type_names AS value
SELECT (SELECT list(name) FROM op_biolink('class', {{entity_types}}))

-- step: key_tokens AS list
SELECT unnest(coalesce({{entity_ids}}, []))

-- step: key_scope AS value
SELECT coalesce({{resources}}, {{sources}})

-- steps: entity_keys

-- step: keys_filter AS value
SELECT CASE WHEN op_unset({{entity_ids}}) THEN NULL ELSE coalesce({{entity_keys}}, []) END

-- step: q_upper AS value
SELECT CASE WHEN {{ql}} = '' THEN NULL
    ELSE left({{ql}}, length({{ql}}) - 1) || chr(unicode(right({{ql}}, 1)) + 1) END

-- step: curie_ns AS value
SELECT CASE WHEN regexp_full_match({{q}}, '([A-Za-z][\w.-]*)[:|](.+)')
    THEN lower(regexp_extract({{q}}, '^([A-Za-z][\w.-]*)[:|](.+)$', 1)) END

-- step: curie_id AS value
SELECT CASE WHEN regexp_full_match({{q}}, '([A-Za-z][\w.-]*)[:|](.+)')
    THEN regexp_extract({{q}}, '^([A-Za-z][\w.-]*)[:|](.+)$', 2) END

-- repeat: hits_limit = greatest(4 * {{count_}}, 64) next {{hits_limit}} * 4 until search_done max 12

-- step: hits AS table when list_contains({{phases}}, 'prefix')
-- term hits by entity row: the best rank, shortest and least term, and the namespaces
-- an entity must have (a CURIE's), or any
SELECT _eid, min(match_rank) AS match_rank, min(length(term)) AS term_length, min(term) AS term,
    list(DISTINCT namespace) FILTER (WHERE namespace IS NOT NULL) AS namespaces,
    bool_or(namespace IS NULL) AS any_namespace
FROM (
    SELECT _eid, CASE WHEN term = {{ql}} THEN (CASE WHEN kind = 'label' THEN 1 ELSE 2 END)
            ELSE (CASE WHEN kind = 'label' THEN 3 ELSE 4 END) END AS match_rank,
        term, NULL::VARCHAR AS namespace
    FROM entity_term
    WHERE term >= {{ql}} AND term < {{q_upper}} AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))
    UNION ALL
    SELECT _eid, CASE WHEN term = lower({{curie_id}}) THEN (CASE WHEN kind = 'label' THEN 1 ELSE 2 END)
            ELSE (CASE WHEN kind = 'label' THEN 3 ELSE 4 END) END,
        term, {{curie_ns}}
    FROM entity_term
    WHERE {{curie_ns}} IS NOT NULL AND term = lower({{curie_id}})
        AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))
    UNION ALL
    SELECT _eid, 2, {{ql}}, NULL FROM entity_group
    WHERE {{curie_ns}} IS NOT NULL AND kind = 'reference' AND group_key = {{curie_ns}} || ':' || {{curie_id}}
        AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))
)
GROUP BY _eid
ORDER BY match_rank, term_length, term, _eid
LIMIT {{hits_limit}}

-- step: hit_ids AS list
SELECT _eid FROM __hits ORDER BY _eid

-- step: hit_rows AS table
-- only what ranking and the filters read
SELECT _eid, entity_key, namespace, _resource, entity_type, taxon, reference_entity_key
FROM entity WHERE _eid IN ({{hit_ids*}})

-- step: ranked AS table
-- each matching entity's best (rank, term length, term) among the hits read
SELECT e.entity_key, min(struct_pack(r := h.match_rank, l := h.term_length, t := h.term)) AS best
FROM __hits h JOIN __hit_rows e USING (_eid)
WHERE op_entity_ok(e, {{scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, {{keys_filter}})
    AND (h.any_namespace OR list_contains(h.namespaces, lower(coalesce(e.namespace, ''))))
GROUP BY e.entity_key

-- step: search_done AS value
-- every hit read, or the hits read cover every tie of the last entity needed
SELECT NOT list_contains({{phases}}, 'prefix')
    OR (SELECT count(*) FROM __hits) < {{hits_limit}}
    OR ((SELECT count(*) FROM __ranked) >= {{count_}}
        AND (SELECT struct_pack(r := match_rank, l := term_length) FROM __hits
             ORDER BY match_rank DESC, term_length DESC LIMIT 1)
          > (SELECT struct_pack(r := best.r, l := best.l) FROM __ranked
             ORDER BY best.r, best.l LIMIT 1 OFFSET {{count_}} - 1))

-- end

-- step: candidates AS list
SELECT entity_key FROM __ranked ORDER BY entity_key

-- step: candidate_copies AS table
-- the candidates' merged rows that pass the filters, for their relation counts
SELECT entity_key, relation_count FROM entity_merged
WHERE entity_key IN ({{candidates*}})
    AND op_merged_ok(entity_merged, {{scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, {{keys_filter}})

-- step: prefix_page AS table
SELECT r.entity_key,
    row_number() OVER (ORDER BY r.best.r, r.best.l, c.relations DESC, r.best.t, r.entity_key) AS position
FROM __ranked r
JOIN (SELECT entity_key, sum(coalesce(relation_count, 0)) AS relations FROM __candidate_copies GROUP BY ALL) c
    USING (entity_key)
QUALIFY position <= {{count_}}

-- step: phase AS value
SELECT CASE
    WHEN list_contains({{phases}}, 'prefix') AND (SELECT count(*) FROM __prefix_page) > 0 THEN 'prefix'
    WHEN list_contains({{phases}}, 'contains') THEN 'contains'
    ELSE {{phases}}[1]
END

-- step: scan_page AS table when {{phase}} <> 'prefix'
WITH matched AS (
    SELECT entity_key, label, relation_count FROM entity_by_key
    WHERE op_entity_ok(entity_by_key, {{scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, {{keys_filter}})
        AND CASE {{phase}}
            WHEN 'key' THEN entity_key = {{q}}
            WHEN 'contains' THEN contains(lower(label), {{ql}})
            ELSE true END
), page AS (
    SELECT entity_key, row_number() OVER (ORDER BY
        CASE WHEN {{phase}} = 'contains' THEN length(sort_label) END,
        CASE WHEN {{phase}} = 'contains' THEN relations END DESC,
        CASE WHEN {{phase}} = 'contains' THEN sort_label END,
        entity_key) AS position
    FROM (SELECT entity_key, coalesce(min(label), '') AS sort_label, sum(relation_count) AS relations
          FROM matched GROUP BY entity_key)
    QUALIFY position <= {{count_}}
)
SELECT entity_key, position FROM page

-- step: page AS table
SELECT entity_key, position FROM __prefix_page WHERE {{phase}} = 'prefix'
UNION ALL
SELECT entity_key, position FROM __scan_page

-- step: page_keys AS list
SELECT entity_key FROM __page ORDER BY entity_key

-- step: summaries AS table
-- the page's entities as merged at build
SELECT p.position, op_summary(m.summary) AS summary
FROM __page p JOIN (SELECT entity_key, summary FROM entity_merged WHERE entity_key IN ({{page_keys*}})) m
    USING (entity_key)

-- step: rows AS table
SELECT summary.entityPk AS entity_key, summary.label AS label, summary.entityType AS entity_type,
    summary.canonicalIdentifierType AS namespace, summary.canonicalIdentifier AS identifier,
    summary.taxonomyId AS taxon, summary.relationCount::BIGINT AS relation_count,
    summary.referenceEntityKey AS reference_entity_key, summary.geneReferenceKeys AS gene_reference_keys
FROM __summaries
WHERE position > {{offset_}} AND position <= {{offset_}} + {{limit}}
ORDER BY position

-- step: response AS value
SELECT to_json(struct_pack(
    entities := coalesce((SELECT list(summary ORDER BY position) FROM __summaries
        WHERE position > {{offset_}} AND position <= {{offset_}} + {{limit}}), []),
    nextCursor := CASE WHEN (SELECT count(*) FROM __summaries) > {{offset_}} + {{limit}}
        THEN struct_pack(phase := {{phase}}, "offset" := {{offset_}} + {{limit}}, queryKey := {{query_key}}) END,
    total := {{offset_}} + (SELECT count(*) FROM __summaries WHERE position > {{offset_}} AND position <= {{offset_}} + {{limit}})
        + ((SELECT count(*) FROM __summaries) > {{offset_}} + {{limit}})::INTEGER))
