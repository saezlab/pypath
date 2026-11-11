-- The ontologies' hierarchies as child -> parent edges, with the ends' labels.
SELECT CASE WHEN {{flip}} THEN o.entity_key ELSE s.entity_key END AS child_key,
    CASE WHEN {{flip}} THEN o.label ELSE s.label END AS child_label,
    r.predicate,
    CASE WHEN {{flip}} THEN s.entity_key ELSE o.entity_key END AS parent_key,
    CASE WHEN {{flip}} THEN s.label ELSE o.label END AS parent_label,
    r.resources
FROM omnipath.relation r
JOIN omnipath.entity s ON s.entity_n = r.subject_n
JOIN omnipath.entity o ON o.entity_n = r.object_n
WHERE r.statement_kind = 'ontology' AND r.predicate IN ({{hierarchy_predicates}})
