-- query: relation_evidence_batch
-- For up to 1,000 relations: per evidence source, in first-seen order, the
-- modification sites of either endpoint's molecular form and the PubMed ids its records
-- state (a record without its own annotations has its relation's); what
-- relation_evidence returns for each, from one read. An unknown key has none.
-- param: relation_keys VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- param: protein_entity_keys VARCHAR[] = NULL
-- param: transcript_entity_keys VARCHAR[] = NULL
-- param: isoform_identifiers VARCHAR[] = NULL
-- param: molecular_endpoint_mode VARCHAR = "any"
-- output: json response
-- example: {"relation_keys": ["missing"]}

-- step: checked AS value
SELECT CASE WHEN len(coalesce({{relation_keys}}, [])) <= 1000 THEN true
    ELSE error('invalid: At most 1000 relation_pks per request') END

-- steps: relation_by_keys

-- step: ids AS list
SELECT _rid FROM __found ORDER BY _rid

-- step: evidence_rows AS table
SELECT * FROM relation_evidence WHERE _rid IN ({{ids*}})

-- step: inheriting AS list
SELECT DISTINCT _rid FROM __evidence_rows WHERE annotations IS NULL ORDER BY _rid

-- step: inherited_rows AS table
SELECT * FROM relation_annotation WHERE _rid IN ({{inheriting*}})

-- step: response AS value
WITH inherited AS (
    SELECT _rid, list(value ORDER BY ordinal) AS publications FROM __inherited_rows
    WHERE term = 'publications' AND coalesce(lower(scope), 'relation') = 'relation' GROUP BY ALL
), records AS (
    SELECT f.relation_key, f.resource, f.relation_id, e.ordinal,
        coalesce(nullif(e.source, ''), 'unknown') AS source,
        CASE WHEN e.annotations IS NULL THEN coalesce(i.publications, [])
            ELSE [a.value FOR a IN e.annotations IF a.term = 'publications'
                AND coalesce(lower(a.scope), 'relation') = 'relation'] END AS publications,
        e.subject_molecular_form.modifications AS subject_modifications,
        e.object_molecular_form.modifications AS object_modifications,
        row_number() OVER (PARTITION BY f.relation_key ORDER BY f.resource, f.relation_id, e.ordinal) AS seq
    FROM __evidence_rows e JOIN __found f USING (_rid) LEFT JOIN inherited i USING (_rid)
    WHERE op_evidence_ok(e, {{protein_entity_keys}}, {{transcript_entity_keys}},
            {{isoform_identifiers}}, {{molecular_endpoint_mode}})
        OR (op_unset({{protein_entity_keys}}) AND op_unset({{transcript_entity_keys}})
            AND op_unset({{isoform_identifiers}}))
), sites AS (
    -- each source's sites, first-seen first, without repeats
    SELECT relation_key, source, list(site ORDER BY first) AS sites FROM (
        SELECT relation_key, source, site, min(struct_pack(a := seq, b := side_order, c := n)) AS first FROM (
            SELECT relation_key, source, seq, s.side_order, m.n, struct_pack(side := s.side,
                term := m.mod.term, residue := m.mod.residue, "position" := m.mod.position,
                end_position := m.mod.end_position, description := m.mod.description) AS site
            FROM records,
                (VALUES ('subject', 1), ('object', 2)) AS s(side, side_order),
                LATERAL (SELECT unnest(CASE s.side WHEN 'subject' THEN subject_modifications
                    ELSE object_modifications END) AS mod,
                    generate_subscripts(CASE s.side WHEN 'subject' THEN subject_modifications
                    ELSE object_modifications END, 1) AS n) AS m)
        GROUP BY ALL)
    GROUP BY ALL
), pmids AS (
    SELECT relation_key, source, list(pmid ORDER BY first) AS pmids FROM (
        SELECT relation_key, source, pmid, min(struct_pack(a := seq, b := p, c := q)) AS first FROM (
            SELECT relation_key, source, seq, v.p, w.q,
                regexp_replace(trim(w.part), '^PMID:\s*', '', 'i') AS pmid
            FROM records,
                LATERAL (SELECT unnest(publications) AS value, generate_subscripts(publications, 1) AS p) AS v,
                LATERAL (SELECT unnest(regexp_split_to_array(coalesce(v.value, ''), '[;,|]')) AS part,
                    generate_subscripts(regexp_split_to_array(coalesce(v.value, ''), '[;,|]'), 1) AS q) AS w)
        WHERE regexp_full_match(pmid, '[0-9]+')
        GROUP BY ALL)
    GROUP BY ALL
), sources AS (
    SELECT relation_key, source, min(seq) AS seq FROM records GROUP BY ALL
), summaries AS (
    SELECT s.relation_key, list(struct_pack(source := s.source,
            sites := coalesce(t.sites, []), pmids := coalesce(p.pmids, [])) ORDER BY s.seq) AS items
    FROM sources s LEFT JOIN sites t USING (relation_key, source) LEFT JOIN pmids p USING (relation_key, source)
    GROUP BY ALL
)
SELECT to_json(struct_pack(evidence := coalesce(
    (SELECT map_from_entries(list(struct_pack(k := k, v := coalesce(s.items, [])) ORDER BY i))
     FROM (SELECT k, min(i) AS i FROM (SELECT unnest({{relation_keys}}) AS k,
               generate_subscripts({{relation_keys}}, 1) AS i) GROUP BY k) keys
     LEFT JOIN summaries s ON s.relation_key = keys.k),
    map([], []))))
