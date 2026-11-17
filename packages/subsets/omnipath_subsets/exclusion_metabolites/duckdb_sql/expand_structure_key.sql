-- Variants of a listed member: charge, tautomer and stereo forms share the
-- connectivity-level structure key with the form a curator listed.
--
-- A variant is added to the classes of the member it shares a key with, as a
-- row of its own with `match_kind = 'structure_key'` and the member it came
-- from. Nothing is merged: the variant stays a different entity.
--
-- An entity that matched in its own right in a class gets no key row there. Where
-- several members share the key, the lowest entity key is named, so a rebuild
-- names the same one.

INSERT INTO exclusion_metabolite
SELECT DISTINCT ON (variant.entity_key, member.exclusion_class, member.assignment_source)
       variant.entity_key,
       member.exclusion_class,
       member.assignment_source,
       member.origin,
       'structure_key',
       NULL,
       NULL,
       'shares structure key ' || listed.group_key,
       member.entity_key,
       member.list_version,
       member.build_id
FROM exclusion_metabolite member
JOIN structure_group listed ON listed.entity_key = member.entity_key
JOIN structure_group variant
  ON variant.group_key = listed.group_key AND variant.entity_key <> member.entity_key
WHERE member.match_kind = 'entity'
  AND NOT EXISTS (
    SELECT 1 FROM exclusion_metabolite own
    WHERE own.entity_key = variant.entity_key
      AND own.exclusion_class = member.exclusion_class
      AND own.match_kind = 'entity'
  )
ORDER BY variant.entity_key, member.exclusion_class, member.assignment_source, member.entity_key;
