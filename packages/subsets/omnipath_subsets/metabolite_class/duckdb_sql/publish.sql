-- The rows of `metabolite_class`, with the readable columns and the supporting resources.
--
-- The HMDB name comes first because HMDB uses the common name (`Cortisol`). The release
-- label is the one of the first resource in name order, which for many molecules is
-- ChEBI's systematic name.
-- The rows of the published molecules only, from each resource: the merged views of the
-- release would aggregate every entity and identifier first.
CREATE TEMP TABLE mc_entity AS
SELECT e.resource, e.entity_id, e.entity_key, e.namespace, e.identifier, e.label
FROM resource_entity e
WHERE e.entity_key IN (SELECT entity_key FROM mc_assigned);

CREATE TEMP TABLE mc_ident AS
SELECT e.entity_key, i.ns, i.id
FROM resource_identifier i
JOIN mc_entity e ON e.resource = i.resource AND e.entity_id = i.entity_id
WHERE i.ns IN ('inchikey', 'hmdb')
UNION
SELECT entity_key, namespace, identifier
FROM mc_entity
WHERE namespace IN ('inchikey', 'hmdb');

INSERT INTO metabolite_class
WITH named AS (
  SELECT entity_key,
         min(label) FILTER (WHERE resource = 'hmdb') AS hmdb_label,
         min(label) FILTER (WHERE resource = 'chebi') AS chebi_label,
         -- The release label: the first resource in name order that has one.
         arg_min(label, resource) FILTER (WHERE label IS NOT NULL) AS release_label,
         arg_min(namespace, resource) FILTER (WHERE namespace IS NOT NULL) AS namespace,
         arg_min(identifier, resource) FILTER (WHERE identifier IS NOT NULL) AS identifier
  FROM mc_entity
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
       coalesce(n.hmdb_label, n.chebi_label, n.release_label),
       CASE WHEN n.namespace = 'inchikey' THEN n.identifier ELSE k.inchikey END,
       k.hmdb,
       a.metabolite_class,
       a.match_kind,
       a.via_entity_key,
       coalesce(s.supported_by, []::VARCHAR[]),
       $build_id
FROM mc_assigned a
LEFT JOIN named n USING (entity_key)
LEFT JOIN keyed k USING (entity_key)
LEFT JOIN support s USING (entity_key);
