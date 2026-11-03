-- The COSMOS prior-knowledge network, as binary edges over the reaction stars.
--
-- One row is one COSMOS edge: a metabolite feeding an enzyme, an enzyme
-- producing a metabolite, or a connector tying a bare identifier to the
-- per-reaction gene node that carries it.
--
-- `omnipath-build` writes this table. It is not the interaction fact table and
-- it never becomes one: a COSMOS gene node is one enzyme *in one reaction*, so
-- an enzyme catalysing five hundred reactions is five hundred nodes. Those are
-- modelling artifacts of the formalism, not entities of the graph, and giving
-- them entity rows would move every count the rest of the build reports.
--
-- The labels are the node identifiers the COSMOS R package expects, and they
-- are the join key a consumer uses. The entity ids beside them are the handle
-- back into the canonical layer, and they are null on the pseudo-node that
-- stands in for an unknown enzyme, because no entity is the pseudo-node.
--
-- There is deliberately **no foreign key to `interaction`**. The interaction
-- projection truncates that table on every derive, and Postgres refuses a
-- TRUNCATE whose dependants it was not asked to name, so a key here would
-- break the core build rather than protect this one.
--
-- The DDL is idempotent: applying it twice changes nothing.

CREATE TABLE IF NOT EXISTS cosmos_edge (
  cosmos_edge_id bigserial PRIMARY KEY,
  -- Which build produced the row, so a consumer can tell a stale edge from a
  -- current one. Read from the manifest, as the MetSigDB substrate reads it.
  build_id text NOT NULL,

  -- The two endpoints, as COSMOS node strings.
  --   Metab__<id>_<compartment>, or Metab__<id> where no compartment is known
  --   Gene<N>__<id>
  --   Gene<N>__orphanReac<reaction id>
  -- with `_rev` appended to the gene node, and only to the gene node, on the
  -- reverse half of a reversible reaction.
  --
  -- N is the reaction's index, and it is part of the label, so the same enzyme
  -- in two reactions carries two labels by design. Both halves of a reversible
  -- reaction share it.
  source_label text NOT NULL,
  target_label text NOT NULL,

  -- The canonical entities behind the labels. Null where the endpoint is the
  -- pseudo-enzyme of a reaction nobody named a catalyst for.
  source_entity_id uuid,
  target_entity_id uuid,

  -- metabolite / protein / reaction.
  source_type text NOT NULL,
  target_type text NOT NULL,

  -- A compartment per endpoint rather than one per edge: a transport reaction
  -- moves the same molecule between two of them, and a single column could not
  -- say that.
  source_compartment text,
  target_compartment text,

  -- The mode of regulation. A reaction edge carries no sign: the direction is
  -- the chemistry, not an activation or an inhibition, so this is always 1.
  mor smallint NOT NULL DEFAULT 1,

  -- catalysis / transport / connector.
  interaction_type text NOT NULL,

  -- Where the edge came from: the event entity, the N-ary header it was read
  -- off, and the index that the gene labels of this reaction carry.
  reaction_entity_id uuid,
  interaction_id uuid,
  reaction_index integer,

  -- The enzyme was unknown and the event stands in for it.
  orphan boolean NOT NULL DEFAULT false,

  -- The reverse half of a reversible reaction: the products feed the enzyme
  -- and the reactants come out of it. The gene node of such a row carries the
  -- `_rev` suffix, so the two halves are distinct nodes of the network rather
  -- than one node with two edge sets.
  reverse boolean NOT NULL DEFAULT false,

  -- What the resources published about the direction: `reversible`,
  -- `left_to_right`, or null where nobody said. Null is a statement of its
  -- own and is never filled in with a default — Rhea, the largest contributor
  -- of reactions, states no direction at all, and a reaction that is
  -- directionless in its source must stay distinguishable from one known to
  -- run one way.
  direction text,

  -- The namespace each label's identifier actually came from. COSMOS wants
  -- ChEBI on the metabolite side and UniProt on the gene side, and the
  -- projection translates into both where a mapping answers; where none does,
  -- the label keeps the identifier the build canonicalised it to and this
  -- column names that namespace instead. It is the whole point of the pair: a
  -- consumer tells a translated label from a fallen-back one without parsing
  -- the string, and a row that is not in the wanted namespace cannot look like
  -- a row that is.
  source_id_type text,
  target_id_type text,

  -- The resources that contributed the reaction.
  sources text[] NOT NULL,

  CONSTRAINT cosmos_edge_mor_check CHECK (mor IN (-1, 1)),
  CONSTRAINT cosmos_edge_interaction_type_check
    CHECK (interaction_type IN ('catalysis', 'transport', 'connector')),
  CONSTRAINT cosmos_edge_node_type_check
    CHECK (
      source_type IN ('metabolite', 'protein', 'reaction')
      AND target_type IN ('metabolite', 'protein', 'reaction')
    ),
  CONSTRAINT cosmos_edge_direction_check
    CHECK (direction IN ('reversible', 'left_to_right')),
  -- Only a reaction stated reversible has a reverse half.
  CONSTRAINT cosmos_edge_reverse_check
    CHECK (NOT reverse OR direction = 'reversible')
);

-- A consumer walks the network by label, which is what makes both label
-- indexes load bearing rather than one of them redundant.
CREATE INDEX IF NOT EXISTS cosmos_edge_source_label_idx
  ON cosmos_edge (source_label);

CREATE INDEX IF NOT EXISTS cosmos_edge_target_label_idx
  ON cosmos_edge (target_label);

CREATE INDEX IF NOT EXISTS cosmos_edge_reaction_idx
  ON cosmos_edge (reaction_entity_id);

CREATE INDEX IF NOT EXISTS cosmos_edge_interaction_type_idx
  ON cosmos_edge (interaction_type);

-- Beyond the four an edge query needs: the build's own counting step asks
-- which headers reached the projection, and that question keys on the header.
CREATE INDEX IF NOT EXISTS cosmos_edge_interaction_idx
  ON cosmos_edge (interaction_id);

COMMENT ON TABLE cosmos_edge IS
  'COSMOS prior-knowledge network: reaction stars projected to binary '
  'metabolite/enzyme edges, with per-reaction gene pseudo-nodes.';
