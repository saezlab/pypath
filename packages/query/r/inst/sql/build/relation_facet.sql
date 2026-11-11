-- Each combination of the relation filter columns (and resources), with the largest
-- evidence count and the number of its relations: which filter options a selection
-- leaves, without reading the relations.
SELECT {{facet_columns}}, resources, max(evidence_count) AS evidence_count, count(*) AS relation_count
FROM omnipath.relation GROUP BY ALL
