-- query: metsigdb_groups
-- One page of MetSigDB sets, grouped by resource, each with its members. Every count is in
-- sets: a page of N sets never splits one across a boundary. A set is never cut: sets are
-- added whole until the next would pass `row_ceiling` members, and a set larger than the
-- ceiling is returned complete and alone. The exclusion thins `returned` and leaves the set
-- list alone: a set that loses every member is returned with none.
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
-- param: set_fields VARCHAR[] = ["set", "set_label", "set_type", "set_sub_type", "set_size"]
-- param: member_fields VARCHAR[] = ["entity", "metabolite_label", "inchikey", "hmdb", "pubchem", "chebi", "kegg"]
-- param: limit INTEGER = 1000
-- param: offset INTEGER = 0
-- param: row_ceiling INTEGER = 50000
-- param: total BOOLEAN = false
-- output: json response
-- example: {"resource": ["Reactome"], "limit": 2}

-- step: present AS value
SELECT CASE WHEN {{include:metsigdb_present x}} THEN true
    ELSE error('invalid: this release has no MetSigDB subset') END

-- step: candidates AS table
SELECT m.resource, m.set_source_id, any_value(m.set_entity_key) AS set_entity_key,
    any_value(m.set_label) AS set_label, any_value(m.set_type) AS set_type,
    any_value(m.set_sub_type) AS set_sub_type, any_value(m.organism) AS organism,
    any_value(m.set_size) AS set_size,
    count(*) FILTER (WHERE {{include:metsigdb_kept m}}) AS returned
FROM metsigdb_membership m
WHERE {{include:metsigdb_conditions m}}
GROUP BY m.resource, m.set_source_id
ORDER BY m.resource, m.set_source_id
LIMIT {{limit}} + 1 OFFSET {{offset}}

-- step: page AS table
SELECT c.*, row_number() OVER w AS n, sum(c.returned) OVER (w ROWS UNBOUNDED PRECEDING) AS running
FROM __candidates c
WINDOW w AS (ORDER BY c.resource, c.set_source_id)

-- step: kept AS table
SELECT * FROM __page WHERE n <= {{limit}} AND (n = 1 OR running <= {{row_ceiling}})

-- step: members AS table
SELECT m.* FROM metsigdb_membership m
JOIN __kept k ON k.resource = m.resource AND k.set_source_id = m.set_source_id
WHERE {{include:metsigdb_conditions m}} AND {{include:metsigdb_kept m}}

-- step: matched AS value when {{total}}
SELECT count(*) FROM (SELECT 1 FROM metsigdb_membership m
    WHERE {{include:metsigdb_conditions m}} GROUP BY m.resource, m.set_source_id)

-- step: response AS value
WITH described AS (
    SELECT k.resource, k.n,
        json_merge_patch(
            op_metsigdb_fields(struct_pack(
                metabolite_entity_key := NULL::VARCHAR, metabolite_label := NULL::VARCHAR,
                metabolite_entity_type := NULL::VARCHAR, metabolite_structure_key := NULL::VARCHAR,
                inchikey := NULL::VARCHAR, smiles := NULL::VARCHAR, hmdb := NULL::VARCHAR,
                pubchem := NULL::VARCHAR, chebi := NULL::VARCHAR, kegg := NULL::VARCHAR,
                has_structure := NULL::BOOLEAN, exclusion_classes := NULL::VARCHAR[],
                resource := k.resource, set_source_id := k.set_source_id,
                set_entity_key := k.set_entity_key, set_label := k.set_label,
                set_type := k.set_type, set_sub_type := k.set_sub_type, organism := k.organism,
                set_size := k.set_size, set_context := NULL::VARCHAR,
                provenance_source := NULL::VARCHAR, provenance_record := NULL::VARCHAR,
                build_id := NULL::VARCHAR), {{set_fields}}),
            json_object('returned', k.returned, 'members', coalesce((
                SELECT list(op_metsigdb_fields(m, {{member_fields}})
                    ORDER BY m.metabolite_entity_key)
                FROM __members m
                WHERE m.resource = k.resource AND m.set_source_id = k.set_source_id), []))) AS one_set
    FROM __kept k
), groups AS (
    SELECT resource, min(n) AS first_n,
        json_object('resource', resource, 'sets', list(one_set ORDER BY n)) AS one_group
    FROM described GROUP BY resource
)
SELECT to_json(struct_pack(
    count := (SELECT count(*) FROM __kept),
    has_more := (SELECT count(*) FROM __candidates) > {{limit}}
        OR (SELECT count(*) FROM __kept) < (SELECT count(*) FROM __page WHERE n <= {{limit}}),
    groups := coalesce((SELECT list(one_group ORDER BY first_n) FROM groups), []),
    total := {{matched}}))
