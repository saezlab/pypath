-- Steps: fast, whether the relation filters select by entity keys or by reference keys
-- alone, without a scope or molecular forms: one range of the endpoint table (the
-- relation_range steps), in page order.

-- step: fast AS value
SELECT (NOT op_unset({{entities}})) <> (NOT op_unset({{reference_entity_keys}}))
    AND op_unset({{scope_entity_ids}}) AND op_unset({{protein_entity_keys}})
    AND op_unset({{transcript_entity_keys}}) AND op_unset({{isoform_identifiers}})
