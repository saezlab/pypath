-- query: exclusion_metabolites
-- One page of the exclusion metabolite table: a metabolite that mass spectrometry does not
-- measure, in one class, with the assignment that put it there. A metabolite assigned to a
-- class by two sources has two rows. The release's exclusion subset holds the table; a
-- release without it has no answer.
-- param: classes VARCHAR[] = NULL
-- param: sources VARCHAR[] = NULL
-- param: origins VARCHAR[] = NULL
-- param: match_kinds VARCHAR[] = NULL
-- param: entity_given BOOLEAN = false
-- param: entity_keys VARCHAR[] = NULL
-- param: inchikeys VARCHAR[] = NULL
-- param: hmdb VARCHAR[] = NULL
-- param: chebi VARCHAR[] = NULL
-- param: kegg VARCHAR[] = NULL
-- param: pubchem VARCHAR[] = NULL
-- param: fields VARCHAR[] = ["entity", "metabolite_label", "inchikey", "hmdb", "pubchem", "chebi", "kegg", "exclusion_class", "assignment_source", "origin", "match_kind"]
-- param: limit INTEGER = 1000
-- param: offset INTEGER = 0
-- param: total BOOLEAN = false
-- output: json response
-- example: {"classes": ["cofactors"], "limit": 5}

-- step: present AS value
SELECT CASE WHEN {{include:exclusion_present x}} THEN true
    ELSE error('invalid: this release has no exclusion metabolite subset') END

-- step: ident AS table
-- The excluded metabolites, named as the MetSigDB table names its metabolites.
SELECT e.entity_key, any_value(e.label) AS label, any_value(e.entity_type) AS entity_type,
    max(i.id) FILTER (WHERE i.ns = 'inchikey') AS inchikey,
    max(CASE WHEN regexp_matches(i.id, '^\s*(?i:HMDB)[\s:_-]*[0-9]+\s*$')
             THEN 'HMDB' || lpad(ltrim(regexp_extract(i.id, '([0-9]+)\s*$', 1), '0'), 7, '0')
             ELSE i.id END) FILTER (WHERE i.ns = 'hmdb') AS hmdb,
    max(i.id) FILTER (WHERE i.ns = 'pubchem') AS pubchem,
    max(CASE WHEN regexp_matches(i.id, '^(?i:CHEBI):[0-9]+$') THEN substr(i.id, 7) ELSE i.id END)
        FILTER (WHERE i.ns = 'chebi') AS chebi,
    max(i.id) FILTER (WHERE i.ns = 'kegg') AS kegg
FROM (SELECT DISTINCT entity_key FROM exclusion_metabolite) k
JOIN entity e ON e.entity_key = k.entity_key
LEFT JOIN entity_identifier i ON i.entity_n = e.entity_n
GROUP BY e.entity_key

-- step: page AS table
SELECT x.*, i.label, i.entity_type, i.inchikey, i.hmdb, i.pubchem, i.chebi, i.kegg,
    row_number() OVER (ORDER BY x.exclusion_class, x.entity_key, x.assignment_source,
        x.match_kind, coalesce(x.seed_identifier::VARCHAR, ''), coalesce(x.via_entity_key::VARCHAR, '')) AS n
FROM (
    SELECT x.* FROM exclusion_metabolite x
    WHERE {{include:exclusion_conditions x}}
    ORDER BY x.exclusion_class, x.entity_key, x.assignment_source, x.match_kind,
        coalesce(x.seed_identifier::VARCHAR, ''), coalesce(x.via_entity_key::VARCHAR, '')
    LIMIT {{limit}} + 1 OFFSET {{offset}}
) x
LEFT JOIN __ident i ON i.entity_key = x.entity_key

-- step: matched AS value when {{total}}
SELECT count(*) FROM exclusion_metabolite x WHERE {{include:exclusion_conditions x}}

-- step: response AS value
SELECT to_json(struct_pack(
    count := count(*) FILTER (WHERE n <= {{limit}}),
    has_more := count(*) > {{limit}},
    rows := coalesce(list(op_exclusion_fields(p, {{fields}}) ORDER BY n) FILTER (WHERE n <= {{limit}}), []),
    total := {{matched}}))
FROM __page p
