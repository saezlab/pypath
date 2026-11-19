-- The exclusion table's filters, as a condition on the assignment row {{input}}: the class,
-- the source and kind of the assignment, and the metabolite by any identifier it may be.
-- Reads the step table `ident` (one row per excluded metabolite with its identifiers).
(op_unset({{classes}}) OR list_contains({{classes}}, {{input}}.exclusion_class))
    AND (op_unset({{sources}}) OR list_contains({{sources}}, {{input}}.assignment_source))
    AND (op_unset({{origins}}) OR list_contains({{origins}}, {{input}}.origin))
    AND (op_unset({{match_kinds}}) OR list_contains({{match_kinds}}, {{input}}.match_kind))
    AND (NOT {{entity_given}}
        OR {{input}}.entity_key IN ({{entity_keys*}})
        OR {{input}}.entity_key IN (SELECT i.entity_key FROM __ident i
            WHERE i.inchikey IN ({{inchikeys*}}) OR i.hmdb IN ({{hmdb*}})
                OR i.chebi IN ({{chebi*}}) OR i.kegg IN ({{kegg*}}) OR i.pubchem IN ({{pubchem*}})))
