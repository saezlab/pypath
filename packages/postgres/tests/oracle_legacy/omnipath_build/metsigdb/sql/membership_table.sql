-- MetSigDB unified membership substrate (cycle 010).
--
-- One row is one canonical metabolite entity in one resource-native set.
-- `omnipath-build` writes this table. `omnipath-metabo` reads it and nothing
-- else for /sets/metsigdb.
--
-- The name is unprefixed, following the derived-table convention of this
-- package (canonical_entity, gene_protein_representative). The metabo_* prefix
-- marks tables written by the omnipath-metabo post-build step, and this table
-- is written by the core build.
--
-- The DDL is idempotent: applying it twice changes nothing.

CREATE TABLE IF NOT EXISTS metsigdb_membership (
  -- Row identity. A rebuild refreshes in place, and never duplicates.
  resource              text  NOT NULL,
  set_source_id         text  NOT NULL,
  metabolite_entity_id  uuid  NOT NULL,

  -- Metabolite side, projected from the canonical entity layer.
  metabolite_label       text NOT NULL,
  metabolite_entity_type text NOT NULL,
  inchikey               text,
  -- The connectivity-level InChIKey block: charge, stereo and tautomer
  -- variants of one molecule share it. The anchor stays the entity, so
  -- identity is unchanged; this is what lets a consumer join a membership
  -- from one resource to the same molecule in another when the two picked
  -- different protonation states.
  metabolite_structure_key text,
  smiles                 text,
  hmdb                   text,
  pubchem                text,
  chebi                  text,
  kegg                   text,

  -- Set side. `set_entity_id` is the canonical entity behind the set, and it
  -- is what a name joins on: matching `set_source_id` against
  -- entity.canonical_identifier is unsafe, because MACdb trait ids are bare
  -- integers that collide with ChEBI ids in that same column.
  --
  -- `set_label` comes from entity_ontology_term. Reactome, MACdb and
  -- ClassyFire name every set. KEGG and WikiPathways carry no name in this
  -- build, so theirs stay null.
  --
  -- `organism` is read from the taxonomy the source recorded on its set
  -- evidence, never guessed from the identifier.
  set_entity_id uuid,
  set_label     text,
  set_type    text   NOT NULL,
  -- A finer semantic than `set_type`, where the source publishes one. MACdb
  -- takes it from the trait type it records, KEGG marks its whole-metabolism
  -- overview maps. A resource with no structured answer leaves it null rather
  -- than inferring one from a label.
  set_sub_type text,
  organism    bigint,
  set_size    integer NOT NULL,
  set_context jsonb,

  -- Provenance. `provenance_source` cites the pypath.inputs_v2 module and the
  -- commit that parsed the resource, read from the `resources` table.
  provenance_source text NOT NULL,
  provenance_record jsonb,
  build_id          text NOT NULL,

  CONSTRAINT metsigdb_membership_pkey
    PRIMARY KEY (resource, set_source_id, metabolite_entity_id),
  CONSTRAINT metsigdb_membership_set_type_check
    CHECK (set_type IN ('disease', 'pathway', 'chemical_class')),
  CONSTRAINT metsigdb_membership_resource_check
    CHECK (resource IN ('KEGG', 'Reactome', 'WikiPathways', 'MACdb', 'ClassyFire')),
  CONSTRAINT metsigdb_membership_set_size_check
    CHECK (set_size > 0)
);

-- Columns added after the first publication, for a database that already
-- carries the table. CREATE TABLE IF NOT EXISTS is a no-op there, so the new
-- columns arrive here, and they arrive before the indexes that read them.
ALTER TABLE metsigdb_membership
  ADD COLUMN IF NOT EXISTS set_entity_id uuid,
  ADD COLUMN IF NOT EXISTS metabolite_structure_key text,
  ADD COLUMN IF NOT EXISTS set_sub_type text;

-- The five shared filter columns of the API contract, and nothing else. The
-- serving layer filters on these alone, so these alone are indexed.
--
-- `resource` and `set_source_id` lead the primary key, so a filter on either
-- already has an index. The three below cover the rest.
CREATE INDEX IF NOT EXISTS metsigdb_membership_metabolite_idx
  ON metsigdb_membership (metabolite_entity_id);

CREATE INDEX IF NOT EXISTS metsigdb_membership_set_type_idx
  ON metsigdb_membership (set_type);

-- Partial: two resources fill the sub-type, so the index stays small and a
-- sub-type filter never scans the ClassyFire bulk.
CREATE INDEX IF NOT EXISTS metsigdb_membership_set_sub_type_idx
  ON metsigdb_membership (set_sub_type) WHERE set_sub_type IS NOT NULL;

-- Cross-resource grouping runs on the structure key, which is the point of
-- publishing it.
CREATE INDEX IF NOT EXISTS metsigdb_membership_structure_idx
  ON metsigdb_membership (metabolite_structure_key);

-- Partial: Reactome and WikiPathways rows carry an organism, so the index stays
-- small and an `organism` filter never scans the null-organism bulk.
CREATE INDEX IF NOT EXISTS metsigdb_membership_organism_idx
  ON metsigdb_membership (organism) WHERE organism IS NOT NULL;

-- Paging needs no index of its own. The API contract orders by the row
-- identity, and the primary key already indexes exactly those three columns in
-- that order.

COMMENT ON TABLE metsigdb_membership IS
  'MetSigDB v1: canonical metabolite memberships in resource-native sets from '
  'KEGG, Reactome, WikiPathways, MACdb and ClassyFire (cycle 010).';

-- Cycle 012: identifier lookup and case-insensitive vocabulary filters.
--
-- Partial or full is decided per column by measurement, not by habit. The
-- identifier columns are sparse **per metabolite** and populated **per row**,
-- and an index is built per row: ClassyFire is 98.4 per cent of the substrate
-- carries HMDB and InChIKey for nearly every metabolite it publishes. So
-- `WHERE ... IS NOT NULL` excludes one row in a hundred of `inchikey` and
-- `hmdb`, and a little over a quarter of `pubchem` — harmless but pointless —
-- while it excludes four fifths of `chebi`
-- and nineteen twentieths of `kegg`.
CREATE INDEX IF NOT EXISTS metsigdb_membership_inchikey_idx
  ON metsigdb_membership (inchikey);

CREATE INDEX IF NOT EXISTS metsigdb_membership_hmdb_idx
  ON metsigdb_membership (hmdb);

CREATE INDEX IF NOT EXISTS metsigdb_membership_pubchem_idx
  ON metsigdb_membership (pubchem);

-- Partial: one row in five carries a ChEBI id, one in twenty a KEGG id.
CREATE INDEX IF NOT EXISTS metsigdb_membership_chebi_idx
  ON metsigdb_membership (chebi) WHERE chebi IS NOT NULL;

CREATE INDEX IF NOT EXISTS metsigdb_membership_kegg_idx
  ON metsigdb_membership (kegg) WHERE kegg IS NOT NULL;

-- `smiles` is deliberately unindexed and is not an accepted `entity` value. It
-- is 111 bytes per row, the widest column in the table, and indexing it would
-- cost more than every other identifier index combined — to serve a lookup
-- nobody performs, because a SMILES string is not canonical across producers
-- and has no shape a recognition rule could claim.

-- Case-insensitive vocabulary matching needs **no index of its own**. The three
-- closed vocabularies are canonicalised to their published spelling before the
-- value reaches SQL, so the predicate stays exact equality and `resource` keeps
-- the primary key's index condition.
--
-- Lowering the column instead was tried and measured: one five-row WikiPathways
-- page went from 6 buffers to 91,442, because `lower(resource) = …` demotes the
-- primary key from an index condition to a filter and the scan discards
-- 3,503,370 ClassyFire rows before reaching a match.
