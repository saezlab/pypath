-- Steps: the relations the evidence count filter reads the selected resources' count of
-- (the relation_fold_low fragment, which relation_conditions includes): the candidates
-- (or the range's relations, when fast) that other resources state too. The merged count
-- is theirs and the other resources', so for the others it is the selected resources' count
-- as it is. Needs the relation_candidates and relation_range steps.

-- step: fold_ids AS list when NOT op_unset({{resources}}) AND coalesce({{min_evidence_count}}, 0) > 1 AND {{keyed}}
SELECT DISTINCT relation_n FROM (
    SELECT relation_n, resources FROM __range_selected WHERE {{fast}}
    UNION ALL
    SELECT relation_n, resources FROM __candidate_rows WHERE NOT {{fast}})
WHERE NOT list_has_all({{resources}}, resources) ORDER BY relation_n
