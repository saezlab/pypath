DELETE FROM exclusion_metabolite x
USING exclusion_veto v
WHERE x.entity_key = v.entity_key AND x.exclusion_class = v.exclusion_class;
