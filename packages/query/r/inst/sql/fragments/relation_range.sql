-- Steps: when fast, the selected keys' range of the endpoint table, one row per relation
-- (__range_selected, with the relation's page-order, filter and facet columns): the
-- entity keys at the ends the gene mode names, or the reference keys at either end, in
-- the resources read. Needs the relation_candidates steps.

-- step: range_keys AS list when {{fast}}
SELECT unnest(CASE WHEN op_unset({{reference_entity_keys}}) THEN coalesce({{entity_key_list}}, [])
    ELSE {{reference_entity_keys}} END) AS k ORDER BY k

-- step: range_rows AS table when {{fast}}
SELECT * FROM relation_endpoint
WHERE key IN ({{range_keys*}})
    AND key_kind = CASE WHEN op_unset({{reference_entity_keys}}) THEN 'entity' ELSE 'reference' END

-- step: range_selected AS table when {{fast}}
-- a relation with both ends among the keys is one relation
WITH mode AS (SELECT CASE WHEN op_unset({{reference_entity_keys}}) THEN {{gene_mode}} ELSE 'any' END AS m),
sides AS (
    SELECT _rid FROM __range_rows, mode GROUP BY _rid
    HAVING any_value(m) <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
)
SELECT DISTINCT ON (r._rid) r.* FROM __range_rows r, mode
WHERE (op_unset({{scope}}) OR list_contains({{scope}}, r._resource))
    AND (m NOT IN ('source', 'source_only') OR r.side = 'subject')
    AND (m NOT IN ('target', 'target_only') OR r.side = 'object')
    AND r._rid IN (SELECT _rid FROM sides)
