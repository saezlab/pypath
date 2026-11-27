-- The tables of the metabolite class product. The two published tables are written as
-- Parquet under the same names. Entities are named by their release `entity_key`.

CREATE TABLE metabolite_class (
  entity_key       VARCHAR NOT NULL PRIMARY KEY,
  -- A readable name: the HMDB name, else the ChEBI name, else the release label.
  label            VARCHAR,
  -- The entity identifier when it is an InChIKey, else the lowest InChIKey it carries.
  inchikey         VARCHAR,
  -- The lowest HMDB accession it carries, in the seven-digit form.
  hmdb             VARCHAR,
  -- The best class of the connectivity group, or of the molecule when it has no group.
  metabolite_class VARCHAR NOT NULL,
  -- `entity` when the molecule's own first claim gives the class; `structure_key` when
  -- another member of its connectivity group does. The two are never merged.
  match_kind       VARCHAR NOT NULL CHECK (match_kind IN ('entity', 'structure_key')),
  -- For a structure_key row, the lowest entity key of the members that give the class.
  via_entity_key   VARCHAR,
  -- Every resource that asserts the molecule itself, sorted; empty for a variant that
  -- no source claims.
  supported_by     VARCHAR[] NOT NULL,
  -- The digest of the release manifest the build ran against.
  build_id         VARCHAR NOT NULL,
  CHECK ((match_kind = 'structure_key') = (via_entity_key IS NOT NULL))
);

CREATE TABLE metabolite_class_report (
  kind             VARCHAR NOT NULL
    CHECK (kind IN ('class_size', 'overlap', 'shared_key', 'no_key')),
  metabolite_class VARCHAR,
  n                INTEGER NOT NULL,
  -- The classes, keys or labels behind the count, as JSON text.
  detail           VARCHAR,
  build_id         VARCHAR NOT NULL
);

-- The classes of the rule file, and one row for each claim a source makes.
CREATE TEMP TABLE mc_class (name VARCHAR NOT NULL, precedence INTEGER NOT NULL);
CREATE TEMP TABLE mc_claim (
  entity_key VARCHAR NOT NULL, metabolite_class VARCHAR NOT NULL, priority INTEGER NOT NULL,
  resource VARCHAR NOT NULL
);
