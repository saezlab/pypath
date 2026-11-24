-- Under a resource selection, the number of each relation's evidence records of the
-- selected resources (relation_n, n) for the relations of a list the includer names (as in
-- {{include:relation_scope_counts page_ids}}: a list step of relation numbers, read by the
-- relation number); no rows without a selection.
SELECT relation_n, count(*) AS n FROM relation_evidence
WHERE NOT op_unset({{resources}}) AND relation_n IN ({{{{input}}*}}) AND list_contains({{resources}}, resource)
GROUP BY relation_n
