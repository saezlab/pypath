-- Each relation touching an ontology term, from the term's side: what the term
-- annotates.
WITH terms AS (SELECT entity_n FROM omnipath.entity WHERE entity_type = 'ontology_class')
SELECT object_n AS term_n, subject_n AS entity_n, relation_n, resources FROM omnipath.relation
WHERE object_n IN (SELECT entity_n FROM terms)
UNION ALL
SELECT subject_n, object_n, relation_n, resources FROM omnipath.relation
WHERE subject_n IN (SELECT entity_n FROM terms)
