-- Exclusion metabolites: chemical entities that mass spectrometry does not
-- measure, in classes a consumer can leave out of a set.
--
-- One row is one entity in one class from one source of the assignment.
-- `omnipath-build` writes both tables and rewrites them whole on every build,
-- so a seed that was removed leaves no row behind. `omnipath-metabo` and the
-- COSMOS subset read them.
--
-- The names are unprefixed, as metsigdb_membership is.
--
-- The DDL is idempotent: applying it twice changes nothing.

CREATE TABLE IF NOT EXISTS exclusion_metabolite (
  entity_id         uuid NOT NULL,
  exclusion_class   text NOT NULL
    CONSTRAINT exclusion_metabolite_class_check
    CHECK (exclusion_class IN ('ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors',
            'cofactors_extended')),
  assignment_source text NOT NULL
    CONSTRAINT exclusion_metabolite_source_check
    CHECK (assignment_source IN
      ('metaproviz_seed', 'project_seed', 'override_add', 'chebi_hierarchy', 'chebi_role')),
  -- `curated` for the three seed sources, `derived` for the two ontology ones.
  origin            text NOT NULL
    CONSTRAINT exclusion_metabolite_origin_check
    CHECK (origin IN ('curated', 'derived')),
  -- `structure_key` marks a variant of a listed member (a charge, tautomer or
  -- stereo form) that shares its structure key. The variant is a different
  -- entity and is never merged with the member.
  match_kind        text NOT NULL
    CONSTRAINT exclusion_metabolite_match_check
    CHECK (match_kind IN ('entity', 'structure_key')),
  -- The identifier the assignment came from. Null for a derived row.
  seed_namespace    text,
  seed_identifier   text,
  -- The seed line for a curated row; the ontology term and path for a derived one.
  evidence          text NOT NULL,
  -- For a structure_key row, the listed member whose key it shares.
  via_entity_id     uuid,
  -- Commit of the seed and override files, or uncommitted:<sha256>.
  list_version      text NOT NULL,
  -- The digest of the release manifest the build ran against.
  build_id          text NOT NULL,
  PRIMARY KEY (entity_id, exclusion_class, assignment_source),
  CONSTRAINT exclusion_metabolite_key_member_check
    CHECK (match_kind <> 'structure_key' OR via_entity_id IS NOT NULL),
  CONSTRAINT exclusion_metabolite_origin_source_check
    CHECK ((origin = 'curated') =
           (assignment_source IN ('metaproviz_seed', 'project_seed', 'override_add')))
);

-- A table made before `cofactors_extended` existed carries the narrower check.
DO $$
BEGIN
  IF EXISTS (SELECT 1 FROM pg_constraint
             WHERE conrelid = 'exclusion_metabolite'::regclass
               AND conname = 'exclusion_metabolite_class_check'
               AND pg_get_constraintdef(oid) NOT LIKE '%cofactors_extended%') THEN
    ALTER TABLE exclusion_metabolite DROP CONSTRAINT exclusion_metabolite_class_check;
    ALTER TABLE exclusion_metabolite ADD CONSTRAINT exclusion_metabolite_class_check
      CHECK (exclusion_class IN ('ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors',
                                 'cofactors_extended'));
  END IF;
END $$;

-- The access path of the route's anti-join: one class at a time, by entity.
CREATE INDEX IF NOT EXISTS exclusion_metabolite_class_entity_idx
  ON exclusion_metabolite (exclusion_class, entity_id);

CREATE TABLE IF NOT EXISTS exclusion_metabolite_report (
  kind            text NOT NULL
    CONSTRAINT exclusion_metabolite_report_kind_check
    CHECK (kind IN ('seed_coverage', 'class_size', 'seed_not_derived', 'derived_not_seeded',
                    'shared_key', 'no_key', 'veto_applied', 'uncommitted_inputs',
                    'cofactor_reactions', 'cofactor_cross_check')),
  exclusion_class text,
  namespace       text,
  resolved        integer,
  fallback_only   integer,
  unresolved      integer,
  -- The identifiers, labels or keys behind the count.
  detail          jsonb,
  build_id        text NOT NULL
);

-- For a person reading the list: every row with the entity's label and its
-- common identifiers.
CREATE OR REPLACE VIEW exclusion_metabolite_members AS
SELECT x.exclusion_class,
       x.entity_id,
       coalesce(e.label, e.canonical_identifier) AS label,
       ids.hmdb,
       ids.pubchem,
       ids.chebi,
       ids.kegg,
       x.assignment_source,
       x.origin,
       x.match_kind,
       x.seed_namespace,
       x.seed_identifier,
       x.via_entity_id
FROM exclusion_metabolite x
JOIN entity e ON e.entity_id = x.entity_id
LEFT JOIN LATERAL (
  SELECT
    -- HMDB ships two accession widths; show the seven-digit form.
    max(CASE WHEN v.value ~ '^HMDB[0-9]+$'
             THEN 'HMDB' || lpad(substring(v.value FROM 5), 7, '0')
             ELSE v.value END) FILTER (WHERE it.name = 'Hmdb:OM:0004')         AS hmdb,
    max(v.value) FILTER (WHERE it.name = 'Pubchem Compound:OM:0002')           AS pubchem,
    -- ChEBI is shown as the namespace-local accession.
    max(CASE WHEN v.value ~* '^CHEBI:[0-9]+$'
             THEN substring(v.value FROM 7) ELSE v.value END)
      FILTER (WHERE it.name = 'Chebi:MI:0474')                                 AS chebi,
    max(v.value) FILTER (WHERE it.name = 'Kegg Compound:MI:2012')              AS kegg
  FROM entity_identifier_lookup l
  JOIN identifier_evidence v ON v.identifier_id = l.identifier_id
  JOIN vocab_identifier_type it ON it.identifier_type_id = v.identifier_type_id
  WHERE l.entity_id = x.entity_id
) ids ON true;
