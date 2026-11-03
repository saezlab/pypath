-- Variants of a listed member: charge, tautomer and stereo forms share the
-- connectivity-level structure key (level 1) with the form a curator listed.
--
-- A variant is added to the classes of the member it shares a key with, as a
-- row of its own with `match_kind = 'structure_key'` and the member it came
-- from. Nothing is merged: the variant stays a different entity.
--
-- An entity that matched in its own right in a class gets no key row there. Where
-- several members share the key, the lowest entity id is named, so a rebuild
-- names the same one.

INSERT INTO exclusion_metabolite
  (entity_id, exclusion_class, assignment_source, origin, match_kind, seed_namespace,
   seed_identifier, evidence, via_entity_id, list_version, build_id)
SELECT DISTINCT ON (variant.entity_id, member.exclusion_class, member.assignment_source)
       variant.entity_id,
       member.exclusion_class,
       member.assignment_source,
       member.origin,
       'structure_key',
       NULL,
       NULL,
       'shares structure key ' || listed.group_key,
       member.entity_id,
       member.list_version,
       member.build_id
FROM exclusion_metabolite member
JOIN chemical_resolution_group_member listed
  ON listed.entity_id = member.entity_id AND listed.level_id = 1
JOIN chemical_resolution_group_member variant
  ON variant.group_key = listed.group_key
 AND variant.level_id = 1
 AND variant.entity_id <> member.entity_id
WHERE member.match_kind = 'entity'
  AND NOT EXISTS (
    SELECT 1 FROM exclusion_metabolite own
    WHERE own.entity_id = variant.entity_id
      AND own.exclusion_class = member.exclusion_class
      AND own.match_kind = 'entity'
  )
ORDER BY variant.entity_id, member.exclusion_class, member.assignment_source, member.entity_id;
