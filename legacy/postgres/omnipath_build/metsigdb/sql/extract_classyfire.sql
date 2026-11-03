-- ClassyFire membership extraction: direct classes, plus every ancestor.
--
-- ClassyFire is two sources. HMDB assigns a metabolite to a ChemOnt class.
-- ChemOnt supplies the `is_a` hierarchy that assignment expands over, in
-- `entity_ontology_relation`.
--
-- A metabolite belongs to its direct class and to every broader class above
-- it, so a query for a broad chemical class matches every metabolite
-- classified beneath it. The expansion is the whole point of the resource, and
-- it is also the dominant cost: 1,535,588 direct assignments become 3,453,875
-- published rows, which is 98 parts in a hundred of the substrate.
--
-- The direct assignment wins where a metabolite reaches one class both
-- directly and through a longer path, and the shallowest ancestor path wins
-- otherwise. `set_context` records which of the two a row is.

DROP TABLE IF EXISTS metsigdb_chemont_edge;
DROP TABLE IF EXISTS metsigdb_chemont_ancestor;
DROP TABLE IF EXISTS metsigdb_classyfire_direct;
DROP TABLE IF EXISTS metsigdb_classyfire_pair;
DROP TABLE IF EXISTS metsigdb_stage;

CREATE TEMP TABLE metsigdb_chemont_edge AS
SELECT r.subject_entity_id AS child, r.object_entity_id AS parent
FROM entity_ontology_relation r
JOIN vocab_relation_predicate p
  ON p.relation_predicate_id = r.predicate_id
 AND p.name = 'is_a'
WHERE r.source_id = %(hierarchy_source_id)s;

CREATE INDEX ON metsigdb_chemont_edge (child);
ANALYZE metsigdb_chemont_edge;

-- The transitive closure of the hierarchy, keeping the shortest path to each
-- ancestor. The depth guard stops a cycle in the source from running away;
-- the measured maximum is 11.
CREATE TEMP TABLE metsigdb_chemont_ancestor AS
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

CREATE INDEX ON metsigdb_chemont_ancestor (node);
ANALYZE metsigdb_chemont_ancestor;

-- The direct assignments, one row per class and metabolite.
CREATE TEMP TABLE metsigdb_classyfire_direct AS
SELECT DISTINCT ON (cls.entity_id, res.entity_id)
       cls.entity_id AS class_entity_id,
       res.entity_id AS metabolite_entity_id,
       re.source_id,
       re.dataset_id,
       re.row_id
FROM relation_evidence re
JOIN vocab_relation_predicate p
  ON p.relation_predicate_id = re.predicate_id
 AND p.name = 'associated_with'
JOIN entity cls
  ON cls.entity_id = re.subject_entity_id
 AND cls.entity_type_id = %(set_entity_type_id)s
JOIN entity_evidence_resolution res
  ON res.source_id = re.source_id
 AND res.entity_evidence_id = re.object_entity_evidence_id
 AND res.status_id = 1
JOIN entity met
  ON met.entity_id = res.entity_id
 AND met.entity_type_id = %(chemical_entity_type_id)s
WHERE re.source_id = %(source_id)s
ORDER BY cls.entity_id, res.entity_id, re.row_id;

ANALYZE metsigdb_classyfire_direct;

-- Direct assignments and ancestors together, before the fold. The identifiers
-- and the JSON are built after it, so the sort carries narrow rows.
CREATE TEMP TABLE metsigdb_classyfire_pair AS
SELECT DISTINCT ON (set_entity_id, metabolite_entity_id)
       set_entity_id, metabolite_entity_id, depth, via_entity_id,
       source_id, dataset_id, row_id
FROM (
  SELECT d.class_entity_id AS set_entity_id,
         d.metabolite_entity_id,
         0                  AS depth,
         d.class_entity_id  AS via_entity_id,
         d.source_id, d.dataset_id, d.row_id
  FROM metsigdb_classyfire_direct d
  UNION ALL
  SELECT a.ancestor,
         d.metabolite_entity_id,
         a.depth,
         d.class_entity_id,
         d.source_id, d.dataset_id, d.row_id
  FROM metsigdb_classyfire_direct d
  JOIN metsigdb_chemont_ancestor a ON a.node = d.class_entity_id
) combined
ORDER BY set_entity_id, metabolite_entity_id, depth;

ANALYZE metsigdb_classyfire_pair;

CREATE TEMP TABLE metsigdb_stage AS
SELECT se.canonical_identifier AS set_source_id,
       se.entity_id            AS set_entity_id,
       pr.metabolite_entity_id,
       jsonb_build_object(
         'source_id',  pr.source_id,
         'dataset_id', pr.dataset_id,
         'row_id',     pr.row_id
       ) AS provenance_record,
       NULL::text AS set_sub_type,
       CASE WHEN pr.depth = 0
            THEN jsonb_build_object('assignment', 'direct', 'depth', 0)
            ELSE jsonb_build_object(
                   'assignment', 'ancestor',
                   'depth', pr.depth,
                   'via', ve.canonical_identifier)
       END AS set_context
FROM metsigdb_classyfire_pair pr
JOIN entity se ON se.entity_id = pr.set_entity_id
JOIN entity ve ON ve.entity_id = pr.via_entity_id
LIMIT %(max_records)s;

DROP TABLE metsigdb_chemont_edge;
DROP TABLE metsigdb_chemont_ancestor;
DROP TABLE metsigdb_classyfire_direct;
DROP TABLE metsigdb_classyfire_pair;
