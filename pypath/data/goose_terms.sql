SELECT
    term.name AS go_name,
    term.term_type AS go_aspect,
    term.acc AS go_id
FROM term
%s
GROUP BY term.acc;
