-- The working tables of the exclusion metabolite product. The same names and
-- columns are written as Parquet. Entities are named by their release
-- `entity_key`.

CREATE TABLE exclusion_metabolite (
  entity_key        VARCHAR NOT NULL,
  exclusion_class   VARCHAR NOT NULL
    CHECK (exclusion_class IN ('ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors',
                               'cofactors_extended')),
  assignment_source VARCHAR NOT NULL
    CHECK (assignment_source IN
      ('metaproviz_seed', 'project_seed', 'override_add', 'chebi_hierarchy', 'chebi_role')),
  -- `curated` for the three seed sources, `derived` for the two ontology ones.
  origin            VARCHAR NOT NULL CHECK (origin IN ('curated', 'derived')),
  -- `structure_key` marks a variant of a listed member (a charge, tautomer or
  -- stereo form) that shares its structure key. The variant is a different
  -- entity and is never merged with the member.
  match_kind        VARCHAR NOT NULL CHECK (match_kind IN ('entity', 'structure_key')),
  -- The identifier the assignment came from. Null for a derived row.
  seed_namespace    VARCHAR,
  seed_identifier   VARCHAR,
  -- The seed line for a curated row; the ontology term and path for a derived one.
  evidence          VARCHAR NOT NULL,
  -- For a structure_key row, the listed member whose key it shares.
  via_entity_key    VARCHAR,
  -- Commit of the seed and override files, or uncommitted:<sha256>.
  list_version      VARCHAR NOT NULL,
  -- The digest of the release manifest the build ran against.
  build_id          VARCHAR NOT NULL,
  CHECK (match_kind <> 'structure_key' OR via_entity_key IS NOT NULL),
  CHECK ((origin = 'curated') =
         (assignment_source IN ('metaproviz_seed', 'project_seed', 'override_add')))
);

CREATE TABLE exclusion_metabolite_report (
  kind            VARCHAR NOT NULL
    CHECK (kind IN ('seed_coverage', 'class_size', 'seed_not_derived', 'derived_not_seeded',
                    'shared_key', 'no_key', 'veto_applied', 'uncommitted_inputs',
                    'cofactor_reactions', 'cofactor_cross_check')),
  exclusion_class VARCHAR,
  namespace       VARCHAR,
  resolved        INTEGER,
  fallback_only   INTEGER,
  unresolved      INTEGER,
  -- The identifiers, labels or keys behind the count, as JSON text.
  detail          VARCHAR,
  build_id        VARCHAR NOT NULL
);

-- The ChEBI terms the closures start from, resolved to entities.
CREATE TABLE exclusion_root (exclusion_class VARCHAR, entity_key VARCHAR, curie VARCHAR);
CREATE TABLE exclusion_role_root (exclusion_class VARCHAR, entity_key VARCHAR, curie VARCHAR);

-- The vetoes of the override file, resolved to entities, and a record of what
-- each one removed.
CREATE TABLE exclusion_veto (
  exclusion_class VARCHAR NOT NULL, namespace VARCHAR NOT NULL, identifier VARCHAR NOT NULL,
  entity_key VARCHAR NOT NULL
);
CREATE TABLE exclusion_veto_removed (
  veto_class VARCHAR NOT NULL, namespace VARCHAR NOT NULL, identifier VARCHAR NOT NULL,
  entity_key VARCHAR NOT NULL, assignment_source VARCHAR NOT NULL, match_kind VARCHAR NOT NULL
);
-- Each veto of the file with its note and what it resolved to, for the report.
CREATE TABLE exclusion_veto_listing (
  exclusion_class VARCHAR, namespace VARCHAR, identifier VARCHAR, note VARCHAR,
  resolved BOOLEAN, entities VARCHAR
);
