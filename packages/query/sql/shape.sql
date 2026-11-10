-- Macros the step-file queries share (sql/queries/*.sql): written into every
-- database's query layer, so the API and both clients shape results the same way.

-- True when a list filter is not given: NULL or empty, as the API's "if values:".
CREATE OR REPLACE MACRO op_unset(x) AS (x IS NULL OR len(x) = 0);

-- One endpoint's molecular form matches every given identity constraint.
CREATE OR REPLACE MACRO op_form_ok(form, proteins, transcripts, isoforms) AS (
    (op_unset(proteins) OR list_contains(proteins, form.protein_entity_key))
    AND (op_unset(transcripts) OR list_contains(transcripts, form.transcript_entity_key))
    AND (op_unset(isoforms) OR list_contains(
        isoforms, form.isoform_identifier.ns || ':' || form.isoform_identifier.id))
);

-- An evidence row matches the molecular form filters at the endpoints the mode names.
CREATE OR REPLACE MACRO op_evidence_ok(ev, proteins, transcripts, isoforms, mode) AS (
    CASE mode
        WHEN 'source' THEN op_form_ok(ev.subject_molecular_form, proteins, transcripts, isoforms)
        WHEN 'target' THEN op_form_ok(ev.object_molecular_form, proteins, transcripts, isoforms)
        WHEN 'both' THEN op_form_ok(ev.subject_molecular_form, proteins, transcripts, isoforms)
            AND op_form_ok(ev.object_molecular_form, proteins, transcripts, isoforms)
        ELSE op_form_ok(ev.subject_molecular_form, proteins, transcripts, isoforms)
            OR op_form_ok(ev.object_molecular_form, proteins, transcripts, isoforms)
    END
);

-- Biolink names as stored: predicates and entity types by name, storage name or CURIE.
-- An unknown name is an error naming it.
CREATE OR REPLACE MACRO op_biolink(kind_, names) AS TABLE
SELECT CASE WHEN n.name IS NULL
    THEN error('invalid: Unknown Biolink ' || kind_ || ': ' || u.v) ELSE n.name END AS name
FROM unnest(names::VARCHAR[]) AS u(v)
LEFT JOIN op_names n ON n.kind = kind_ AND n.alias = lower(trim(u.v));


-- An entity row passes the search filters: its resource, type, taxon, reference and key.
CREATE OR REPLACE MACRO op_entity_ok(e, scope, type_names, taxa, references_, keys) AS (
    (op_unset(scope) OR list_contains(scope, e._resource))
    AND (op_unset(type_names) OR list_contains(type_names, e.entity_type))
    AND (op_unset(taxa) OR list_contains(taxa, e.taxon))
    AND (op_unset(references_) OR list_contains(references_, e.reference_entity_key))
    AND (keys IS NULL OR list_contains(keys, e.entity_key))
);

-- Values in first-seen order without repeats (positions are 1-based).
CREATE OR REPLACE MACRO op_ordered_distinct(xs) AS (
    list_transform(
        list_filter(list_transform(xs, lambda x, i: struct_pack(x := x, i := i)),
            lambda s: list_position(xs, s.x) = s.i),
        lambda s: s.x)
);

-- The last value of a relation-scope annotation term, as a dict built from the
-- annotations in order keeps it (interaction_label).
CREATE OR REPLACE MACRO op_last_value(anns, term_) AS (
    list_reverse(list_filter(anns, lambda a: a.scope = 'relation' AND a.term = term_))[1].value
);

-- One relation as the API's relation summary (shape_relation_summary).
CREATE OR REPLACE MACRO op_relation_json(r, anns) AS (
    struct_pack(
        relationPk := coalesce(r.relation_key, ''),
        subjectEntityPk := coalesce(r.subject_entity_key, ''),
        subjectReferenceEntityKey := r.subject_reference_entity_key,
        objectReferenceEntityKey := r.object_reference_entity_key,
        predicate := coalesce(r.predicate, ''),
        displayLabel := CASE
            WHEN r.predicate = 'interacts_with'
                AND op_last_value(anns, 'causal_mechanism_qualifier') = 'binding' THEN 'Binds'
            WHEN r.predicate = 'affects'
                AND op_last_value(anns, 'object_aspect_qualifier') = 'transport'
                AND op_last_value(anns, 'causal_mechanism_qualifier') = 'relocalization'
                AND coalesce(op_last_value(anns, 'object_direction_qualifier'), '') = ''
                THEN 'Transports'
        END,
        -- effect qualifiers, the first value of each
        qualifiers := to_json(map_from_entries(coalesce(list_transform(
            list_filter(list_transform(
                list_filter(anns, lambda a: a.scope = 'relation' AND coalesce(a.value, '') <> ''
                    AND a.term IN ('object_direction_qualifier', 'object_aspect_qualifier',
                        'causal_mechanism_qualifier')),
                lambda a, i: struct_pack(a := a, i := i)),
                lambda s: list_position(list_transform(
                    list_filter(anns, lambda a: a.scope = 'relation' AND coalesce(a.value, '') <> ''
                        AND a.term IN ('object_direction_qualifier', 'object_aspect_qualifier',
                            'causal_mechanism_qualifier')),
                    lambda a: a.term), s.a.term) = s.i),
            lambda s: struct_pack(k := s.a.term, v := s.a.value)), []))),
        objectEntityPk := coalesce(r.object_entity_key, ''),
        relationCategory := r.category,
        interactionClass := r.interaction_class,
        isDirected := coalesce(r.is_directed, false),
        sign := coalesce(r.sign, 0),
        participantTypes := op_ordered_distinct(list_filter(
            [r.subject_type, r.object_type], lambda t: coalesce(t, '') <> '')),
        evidenceCount := coalesce(r.evidence_count, 0),
        sources := list_filter(coalesce(r.sources, []), lambda x: coalesce(x, '') <> '')
    )
);

-- A quantity annotation's fields as one struct, NULL when all are empty.
CREATE OR REPLACE MACRO op_quantity(numeric_value, unit, unit_prefix, binary_relation, source_field, comparator) AS (
    CASE WHEN numeric_value IS NULL AND unit IS NULL AND unit_prefix IS NULL
        AND binary_relation IS NULL AND source_field IS NULL AND comparator IS NULL THEN NULL
    ELSE struct_pack(
        has_numeric_value := numeric_value,
        has_unit := unit,
        has_unit_prefix := unit_prefix,
        has_binary_relation := binary_relation,
        source_field := source_field,
        comparator := comparator) END
);

-- A relation's qualifier columns as relation-scope annotations, in the order of
-- RELATION_QUALIFIERS (what the API shapes a relation from when its annotations are
-- not read).
CREATE OR REPLACE MACRO op_qualifier_annotations(r) AS (
    flatten([
        list_transform(coalesce(r.object_aspect_qualifier, []), lambda v:
            struct_pack(term := 'object_aspect_qualifier', value := v, scope := 'relation')),
        list_transform(coalesce(r.object_direction_qualifier, []), lambda v:
            struct_pack(term := 'object_direction_qualifier', value := v, scope := 'relation')),
        list_transform(coalesce(r.causal_mechanism_qualifier, []), lambda v:
            struct_pack(term := 'causal_mechanism_qualifier', value := v, scope := 'relation'))
    ])
);

-- Annotations of one scope (subject, object, relation: the default) as the API lists
-- them (split_compiled_annotations): term and value (NULL when empty), and the
-- quantity only when there is one. A list of JSON values.
CREATE OR REPLACE MACRO op_scoped_annotations(anns, scope_) AS (
    list_transform(
        list_filter(coalesce(anns, []), lambda a: trim(coalesce(a.term, '')) <> ''
            AND CASE lower(coalesce(a.scope, 'relation'))
                WHEN 'subject' THEN 'subject' WHEN 'object' THEN 'object' ELSE 'relation' END = scope_),
        lambda a: CASE WHEN a.quantity IS NULL
            THEN to_json(struct_pack(term := trim(a.term), value := nullif(a.value, '')))
            ELSE to_json(struct_pack(term := trim(a.term), value := nullif(a.value, ''), quantity := a.quantity))
        END)
);

-- A name's rank for display, least first (omnipath_core.display_names.name_rank): the
-- penalties of the rules it matches and of its length, then its length, then itself.
CREATE OR REPLACE MACRO op_name_rank(name_) AS (
    struct_pack(
        score := coalesce((SELECT sum(penalty) FROM op_name_rules WHERE regexp_matches(trim(name_), pattern)), 0)
            + 25 * (length(trim(name_)) > 80)::INTEGER + 8 * (length(trim(name_)) > 40)::INTEGER
            + 4 * (length(trim(name_)) < 4)::INTEGER,
        length := length(trim(name_)),
        name := trim(name_))
);

-- A JSON object with the fields of another set over it, nulls included (Python's
-- dict.update; json_merge_patch would drop the fields set to null).
CREATE OR REPLACE MACRO op_json_patch(obj, patch) AS (
    to_json(map_concat(json_transform(obj::JSON, '"MAP(VARCHAR, JSON)"'),
        json_transform(patch::JSON, '"MAP(VARCHAR, JSON)"')))
);
