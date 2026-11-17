-- The MetSigDB filters other than the exclusion, as a condition on the membership row
-- {{input}}: the vocabularies (already in their published spelling), the organism, the sets,
-- the metabolite by any identifier it may be, and its structure key.
(op_unset({{resource}}) OR list_contains({{resource}}, {{input}}.resource))
    AND (op_unset({{set_type}}) OR list_contains({{set_type}}, {{input}}.set_type))
    AND (op_unset({{set_sub_type}}) OR list_contains({{set_sub_type}}, {{input}}.set_sub_type))
    AND ({{organism}} IS NULL OR {{input}}.organism = {{organism}})
    AND (op_unset({{sets}}) OR list_contains({{sets}}, {{input}}.set_source_id))
    AND (NOT {{entity_given}}
        OR {{input}}.metabolite_entity_key IN ({{entity_keys*}})
        OR {{input}}.inchikey IN ({{inchikeys*}})
        OR {{input}}.hmdb IN ({{hmdb*}})
        OR {{input}}.chebi IN ({{chebi*}})
        OR {{input}}.kegg IN ({{kegg*}})
        OR {{input}}.pubchem IN ({{pubchem*}}))
    AND (op_unset({{structure_keys}})
        OR list_contains({{structure_keys}}, {{input}}.metabolite_structure_key))
