-- query: metsigdb_rows
-- One page of MetSigDB memberships, one row per metabolite in a set, in the order of the
-- table's key. The release's metsigdb subset holds them; a release without it has no answer.
-- param: resource VARCHAR[] = NULL
-- param: set_type VARCHAR[] = NULL
-- param: set_sub_type VARCHAR[] = NULL
-- param: organism BIGINT = NULL
-- param: sets VARCHAR[] = NULL
-- param: entity_given BOOLEAN = false
-- param: entity_keys VARCHAR[] = NULL
-- param: inchikeys VARCHAR[] = NULL
-- param: hmdb VARCHAR[] = NULL
-- param: chebi VARCHAR[] = NULL
-- param: kegg VARCHAR[] = NULL
-- param: pubchem VARCHAR[] = NULL
-- param: structure_keys VARCHAR[] = NULL
-- param: exclude_classes VARCHAR[] = NULL
-- param: fields VARCHAR[] = ["entity", "metabolite_label", "inchikey", "hmdb", "pubchem", "chebi", "kegg", "resource", "set", "set_label", "set_type", "set_sub_type", "set_size"]
-- param: limit INTEGER = 1000
-- param: offset INTEGER = 0
-- param: total BOOLEAN = false
-- output: json response
-- example: {"resource": ["Reactome"], "limit": 5}

-- step: present AS value
SELECT CASE WHEN {{include:metsigdb_present x}} THEN true
    ELSE error('invalid: this release has no MetSigDB subset') END

-- step: page AS table
SELECT m.*, row_number() OVER (ORDER BY m.resource, m.set_source_id, m.metabolite_entity_key) AS n
FROM (
    SELECT m.* FROM metsigdb_membership m
    WHERE {{include:metsigdb_conditions m}} AND {{include:metsigdb_kept m}}
    ORDER BY m.resource, m.set_source_id, m.metabolite_entity_key
    LIMIT {{limit}} + 1 OFFSET {{offset}}
) m

-- step: matched AS value when {{total}}
SELECT count(*) FROM metsigdb_membership m
WHERE {{include:metsigdb_conditions m}} AND {{include:metsigdb_kept m}}

-- step: response AS value
SELECT to_json(struct_pack(
    count := count(*) FILTER (WHERE n <= {{limit}}),
    has_more := count(*) > {{limit}},
    rows := coalesce(list(op_metsigdb_fields(p, {{fields}}) ORDER BY n) FILTER (WHERE n <= {{limit}}), []),
    total := {{matched}}))
FROM __page p
