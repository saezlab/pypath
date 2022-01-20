SELECT
    term.acc AS go_id,
    ancestor.acc AS ancestor_id
FROM graph_path
INNER JOIN term
    ON graph_path.term2_id = term.id
INNER JOIN term AS ancestor
    ON graph_path.term1_id = ancestor.id
%s
GROUP BY graph_path.id;
