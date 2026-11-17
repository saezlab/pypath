-- The MetSigDB membership table: one row is one canonical metabolite entity in
-- one resource-native set. It is written as Parquet under the same name.
-- Entities are named by their release `entity_key`.

CREATE TABLE metsigdb_membership (
  -- Row identity: (resource, set_source_id, metabolite_entity_key).
  resource                 VARCHAR NOT NULL
    CHECK (resource IN ('KEGG', 'Reactome', 'WikiPathways', 'MACdb', 'ClassyFire')),
  set_source_id            VARCHAR NOT NULL,
  metabolite_entity_key    VARCHAR NOT NULL,

  -- Metabolite side, projected from the release's entities and identifiers.
  metabolite_label         VARCHAR NOT NULL,
  metabolite_entity_type   VARCHAR NOT NULL,
  inchikey                 VARCHAR,
  -- The connectivity-level InChIKey block: charge, stereo and tautomer variants of
  -- one molecule share it. The anchor stays the entity.
  metabolite_structure_key VARCHAR,
  smiles                   VARCHAR,
  hmdb                     VARCHAR,
  pubchem                  VARCHAR,
  chebi                    VARCHAR,
  kegg                     VARCHAR,
  -- Whether the metabolite has an InChIKey, a SMILES or a structure key. A member without
  -- one is a generic class or a name the build could not anchor to a structure.
  has_structure            BOOLEAN NOT NULL,
  -- The exclusion classes of the metabolite (cofactors, ions, ...), sorted; empty when it
  -- is in none. Leaving a class out of a set needs no join.
  exclusion_classes        VARCHAR[] NOT NULL,

  -- Set side. `set_entity_key` is the entity behind the set.
  set_entity_key           VARCHAR,
  set_label                VARCHAR,
  set_type                 VARCHAR NOT NULL
    CHECK (set_type IN ('disease', 'pathway', 'chemical_class')),
  -- A finer semantic than `set_type`, where the source publishes one.
  set_sub_type             VARCHAR,
  organism                 BIGINT,
  set_size                 INTEGER NOT NULL CHECK (set_size > 0),
  set_context              VARCHAR,

  -- Provenance, as JSON text for the record.
  provenance_source        VARCHAR NOT NULL,
  provenance_record        VARCHAR,
  build_id                 VARCHAR NOT NULL
);
