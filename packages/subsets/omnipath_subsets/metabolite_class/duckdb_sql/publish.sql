-- The rows of `metabolite_class`, with the readable columns and the supporting resources.
--
-- The HMDB name comes first because HMDB uses the common name (`Cortisol`). The release
-- label is the one of the first resource in name order, which for many molecules is
-- ChEBI's systematic name.
CREATE TEMP TABLE mc_ident AS
SELECT i.entity_key, i.ns, i.id
FROM release_identifier i
WHERE i.ns IN ('inchikey', 'hmdb')
  AND i.entity_key IN (SELECT entity_key FROM mc_assigned);

INSERT INTO metabolite_class
WITH named AS (
  SELECT entity_key,
         min(label) FILTER (WHERE resource = 'hmdb') AS hmdb_label,
         min(label) FILTER (WHERE resource = 'chebi') AS chebi_label
  FROM resource_entity
  WHERE resource IN ('hmdb', 'chebi')
    AND label IS NOT NULL
    AND entity_key IN (SELECT entity_key FROM mc_assigned)
  GROUP BY entity_key
),
keyed AS (
  SELECT entity_key,
         min(id) FILTER (WHERE ns = 'inchikey') AS inchikey,
         min(CASE WHEN regexp_full_match(id, 'HMDB[0-9]{1,7}')
                  THEN 'HMDB' || lpad(substr(id, 5), 7, '0') END)
           FILTER (WHERE ns = 'hmdb') AS hmdb
  FROM mc_ident
  GROUP BY entity_key
),
support AS (
  SELECT entity_key, list(DISTINCT resource ORDER BY resource) AS supported_by
  FROM mc_claim
  GROUP BY entity_key
)
SELECT a.entity_key,
       coalesce(n.hmdb_label, n.chebi_label, e.label),
       CASE WHEN e.namespace = 'inchikey' THEN e.identifier ELSE k.inchikey END,
       k.hmdb,
       a.metabolite_class,
       a.match_kind,
       a.via_entity_key,
       coalesce(s.supported_by, []::VARCHAR[]),
       $build_id
FROM mc_assigned a
JOIN release_entity e USING (entity_key)
LEFT JOIN named n USING (entity_key)
LEFT JOIN keyed k USING (entity_key)
LEFT JOIN support s USING (entity_key);
