-- Write one resource into the substrate under the deterministic row identity.
--
-- `set_size` is derived here, from the final staged population of each set,
-- never from an upstream total. A source that publishes no count still gets
-- one, and a capped run gets the count of what it actually published.
--
-- The conflict clause is what makes a rebuild idempotent: the identity columns
-- are the primary key, so a second run updates the row in place instead of
-- adding one.
--
-- `set_label` and `organism` come from metsigdb_set_attrs, which reads the set
-- entity rather than the source string.

INSERT INTO metsigdb_membership (
  resource, set_source_id, metabolite_entity_id,
  metabolite_label, metabolite_entity_type, metabolite_structure_key,
  inchikey, smiles, hmdb, pubchem, chebi, kegg,
  set_entity_id, set_label, set_type, set_sub_type, organism, set_size, set_context,
  provenance_source, provenance_record, build_id
)
SELECT %(resource)s,
       s.set_source_id,
       s.metabolite_entity_id,
       p.metabolite_label,
       p.metabolite_entity_type,
       p.metabolite_structure_key,
       p.inchikey, p.smiles, p.hmdb, p.pubchem, p.chebi, p.kegg,
       s.set_entity_id,
       a.set_label,
       %(set_type)s,
       s.set_sub_type,
       a.organism,
       sz.set_size,
       s.set_context,
       %(provenance_source)s,
       s.provenance_record,
       %(build_id)s
FROM metsigdb_stage s
JOIN metsigdb_projection p ON p.entity_id = s.metabolite_entity_id
JOIN metsigdb_set_attrs a ON a.entity_id = s.set_entity_id
JOIN (
  SELECT set_source_id, count(*)::int AS set_size
  FROM metsigdb_stage GROUP BY set_source_id
) sz ON sz.set_source_id = s.set_source_id
ON CONFLICT (resource, set_source_id, metabolite_entity_id) DO UPDATE SET
  metabolite_label       = EXCLUDED.metabolite_label,
  metabolite_entity_type = EXCLUDED.metabolite_entity_type,
  metabolite_structure_key = EXCLUDED.metabolite_structure_key,
  inchikey               = EXCLUDED.inchikey,
  smiles                 = EXCLUDED.smiles,
  hmdb                   = EXCLUDED.hmdb,
  pubchem                = EXCLUDED.pubchem,
  chebi                  = EXCLUDED.chebi,
  kegg                   = EXCLUDED.kegg,
  set_entity_id          = EXCLUDED.set_entity_id,
  set_label              = EXCLUDED.set_label,
  set_type               = EXCLUDED.set_type,
  set_sub_type           = EXCLUDED.set_sub_type,
  organism               = EXCLUDED.organism,
  set_size               = EXCLUDED.set_size,
  set_context            = EXCLUDED.set_context,
  provenance_source      = EXCLUDED.provenance_source,
  provenance_record      = EXCLUDED.provenance_record,
  build_id               = EXCLUDED.build_id;
