-- Under a resource selection and an evidence count filter, the relations with fewer
-- evidence records of the selected resources than the minimum (relation_n), which the
-- relation_conditions fragment drops whatever their merged count: of the relations of a list
-- the includer names (as in {{include:relation_fold_low fold_ids}}), or, without endpoint
-- filters, of every relation. No rows without a selection or a minimum.
SELECT relation_n FROM relation_evidence
WHERE NOT op_unset({{resources}}) AND coalesce({{min_evidence_count}}, 0) > 1 AND {{keyed}}
    AND relation_n IN ({{{{input}}*}}) AND list_contains({{resources}}, resource)
GROUP BY relation_n HAVING count(*) < {{min_evidence_count}}
UNION ALL
SELECT relation_n FROM relation_evidence
WHERE NOT op_unset({{resources}}) AND coalesce({{min_evidence_count}}, 0) > 1 AND NOT {{keyed}}
    AND list_contains({{resources}}, resource)
GROUP BY relation_n HAVING count(*) < {{min_evidence_count}}
