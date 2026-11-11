-- Steps: when fast, the selected keys' range of the endpoint table, one row per relation
-- (__range_selected, with the relation's page-order, filter and facet columns): the
-- entities at the ends the gene mode names, or the reference keys at either end, in
-- the resources read. Needs the relation_candidates steps.

-- step: range_entity_rows AS table when {{fast}} AND op_unset({{reference_entity_keys}})
-- one entity at a time
SELECT * EXCLUDE (entity_n) FROM relation_endpoint WHERE entity_n IN ({{entity_n_list+}})

-- step: range_reference_rows AS table when {{fast}} AND NOT op_unset({{reference_entity_keys}})
SELECT * EXCLUDE (reference_key) FROM relation_reference WHERE reference_key IN ({{reference_entity_keys+}})

-- step: range_selected AS table when {{fast}}
-- a relation with both ends among the keys is one relation
WITH rows_ AS (SELECT * FROM __range_entity_rows UNION ALL BY NAME SELECT * FROM __range_reference_rows),
mode AS (SELECT CASE WHEN op_unset({{reference_entity_keys}}) THEN {{gene_mode}} ELSE 'any' END AS m),
sides AS (
    SELECT relation_n FROM rows_, mode GROUP BY relation_n
    HAVING any_value(m) <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
)
SELECT DISTINCT ON (r.relation_n) r.* FROM rows_ r, mode
WHERE (op_unset({{scope}}) OR list_has_any(r.resources, {{scope}}))
    AND (m NOT IN ('source', 'source_only') OR r.side = 'subject')
    AND (m NOT IN ('target', 'target_only') OR r.side = 'object')
    AND r.relation_n IN (SELECT relation_n FROM sides)
