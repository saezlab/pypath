-- ClassyFire membership extraction: direct classes, plus every ancestor.
--
-- ClassyFire is two sources. HMDB assigns a metabolite to a ChemOnt class. ChemOnt
-- supplies the `subclass_of` hierarchy that assignment expands over.
--
-- A metabolite belongs to its direct class and to every broader class above it, so a
-- query for a broad chemical class matches every metabolite classified beneath it.
--
-- The direct assignment wins where a metabolite reaches one class both directly and
-- through a longer path, and the shallowest ancestor path wins otherwise. `set_context`
-- records which of the two a row is.

CREATE OR REPLACE TEMP TABLE metsigdb_chemont_edge AS
SELECT DISTINCT subject_entity_key AS child, object_entity_key AS parent
FROM release_relation
WHERE resource = $hierarchy_source
  AND statement_kind IN ('ontology', 'ontology_axiom')
  AND predicate = 'subclass_of';

-- The transitive closure of the hierarchy, keeping the shortest path to each ancestor.
-- The depth guard stops a cycle in the source from running away.
CREATE OR REPLACE TEMP TABLE metsigdb_chemont_ancestor AS
WITH RECURSIVE walk (node, ancestor, depth) AS (
  SELECT child, parent, 1 FROM metsigdb_chemont_edge
  UNION ALL
  SELECT w.node, e.parent, w.depth + 1
  FROM walk w
  JOIN metsigdb_chemont_edge e ON e.child = w.ancestor
  WHERE w.depth < 20
)
SELECT node, ancestor, min(depth) AS depth
FROM walk GROUP BY node, ancestor;

-- The direct assignments, one row per class and metabolite.
CREATE OR REPLACE TEMP TABLE metsigdb_classyfire_direct AS
SELECT DISTINCT ON (cls.entity_key, met.entity_key)
       cls.entity_key AS class_entity_key,
       met.entity_key AS metabolite_entity_key,
       re.resource, re.dataset, re.row_id
FROM release_relation_evidence re
JOIN release_entity cls
  ON cls.entity_key IN (re.subject_entity_key, re.object_entity_key)
 AND cls.entity_type = $set_entity_type
 -- HMDB also associates metabolites with disease terms. Only the terms ChemOnt itself
 -- publishes are chemical classes.
 AND cls.entity_key IN (SELECT entity_key FROM resource_entity WHERE resource = $hierarchy_source)
JOIN release_entity met
  ON met.entity_key = CASE WHEN cls.entity_key = re.subject_entity_key
                          THEN re.object_entity_key ELSE re.subject_entity_key END
 AND list_contains($chemical_entity_types, met.entity_type)
WHERE re.resource = $source
  AND re.statement_kind NOT IN ('ontology', 'ontology_axiom')
  AND re.predicate = 'associated_with'
ORDER BY cls.entity_key, met.entity_key, re.row_id;

-- Direct assignments and ancestors together, before the fold.
CREATE OR REPLACE TEMP TABLE metsigdb_classyfire_pair AS
SELECT DISTINCT ON (set_entity_key, metabolite_entity_key)
       set_entity_key, metabolite_entity_key, depth, via_entity_key, resource, dataset, row_id
FROM (
  SELECT d.class_entity_key AS set_entity_key,
         d.metabolite_entity_key,
         0                  AS depth,
         d.class_entity_key AS via_entity_key,
         d.resource, d.dataset, d.row_id
  FROM metsigdb_classyfire_direct d
  UNION ALL
  SELECT a.ancestor,
         d.metabolite_entity_key,
         a.depth,
         d.class_entity_key,
         d.resource, d.dataset, d.row_id
  FROM metsigdb_classyfire_direct d
  JOIN metsigdb_chemont_ancestor a ON a.node = d.class_entity_key
) combined
ORDER BY set_entity_key, metabolite_entity_key, depth, via_entity_key, row_id;

CREATE OR REPLACE TEMP TABLE metsigdb_stage AS
SELECT se.identifier AS set_source_id,
       se.entity_key AS set_entity_key,
       pr.metabolite_entity_key,
       json_object('source', pr.resource, 'dataset', pr.dataset,
                   'row_id', pr.row_id)::VARCHAR AS provenance_record,
       NULL::VARCHAR AS set_sub_type,
       CASE WHEN pr.depth = 0
            THEN json_object('assignment', 'direct', 'depth', 0)
            ELSE json_object('assignment', 'ancestor', 'depth', pr.depth, 'via', ve.identifier)
       END::VARCHAR AS set_context
FROM metsigdb_classyfire_pair pr
JOIN release_entity se ON se.entity_key = pr.set_entity_key
JOIN release_entity ve ON ve.entity_key = pr.via_entity_key;

DROP TABLE metsigdb_chemont_edge;
DROP TABLE metsigdb_chemont_ancestor;
DROP TABLE metsigdb_classyfire_direct;
DROP TABLE metsigdb_classyfire_pair;
