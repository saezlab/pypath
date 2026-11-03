-- The vetoes of the override file, resolved to entities, and a record of what
-- each one removed. Both live for the product transaction only.

DROP TABLE IF EXISTS exclusion_veto;
DROP TABLE IF EXISTS exclusion_veto_removed;

CREATE TEMP TABLE exclusion_veto (
  exclusion_class text NOT NULL,
  namespace       text NOT NULL,
  identifier      text NOT NULL,
  entity_id       uuid NOT NULL
);

CREATE TEMP TABLE exclusion_veto_removed (
  veto_class        text NOT NULL,
  namespace         text NOT NULL,
  identifier        text NOT NULL,
  entity_id         uuid NOT NULL,
  assignment_source text NOT NULL,
  match_kind        text NOT NULL
);
