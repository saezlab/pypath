-- Write one resource into the membership table.
--
-- `set_size` is derived here, from the final staged population of each set, never from
-- an upstream total. A source that publishes no count still gets one.
INSERT INTO metsigdb_membership (
  resource, set_source_id, metabolite_entity_key, metabolite_label, metabolite_entity_type,
  inchikey, metabolite_structure_key, smiles, hmdb, pubchem, chebi, kegg,
  has_structure, exclusion_classes,
  set_entity_key, set_label, set_type, set_sub_type, organism, set_size, set_context,
  provenance_source, provenance_record, build_id
)
SELECT $resource,
       $set_id_prefix || s.set_source_id,
       s.metabolite_entity_key,
       p.metabolite_label,
       p.metabolite_entity_type,
       p.inchikey,
       p.metabolite_structure_key,
       p.smiles, p.hmdb, p.pubchem, p.chebi, p.kegg,
       p.inchikey IS NOT NULL OR p.smiles IS NOT NULL OR p.metabolite_structure_key IS NOT NULL,
       coalesce(x.classes, []::VARCHAR[]),
       s.set_entity_key,
       a.set_label,
       $set_type,
       s.set_sub_type,
       a.organism,
       sz.set_size,
       s.set_context,
       $provenance_source,
       s.provenance_record,
       $build_id
FROM metsigdb_stage s
JOIN metsigdb_projection p ON p.entity_key = s.metabolite_entity_key
JOIN metsigdb_set_attrs a ON a.entity_key = s.set_entity_key
LEFT JOIN (
  SELECT entity_key, list(DISTINCT exclusion_class ORDER BY exclusion_class) AS classes
  FROM exclusion_member GROUP BY entity_key
) x ON x.entity_key = s.metabolite_entity_key
JOIN (
  SELECT set_source_id, count(*)::INTEGER AS set_size
  FROM metsigdb_stage GROUP BY set_source_id
) sz ON sz.set_source_id = s.set_source_id;
