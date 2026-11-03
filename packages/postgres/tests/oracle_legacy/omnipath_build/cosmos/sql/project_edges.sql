-- The reaction edges: every (enzyme, reactant) and every (enzyme, product).
--
-- The reaction star is already assembled. The interaction projection wrote one
-- N-ary header per reaction and one `interaction_party` row per participant,
-- each in the role its resources published, so this step reads roles rather
-- than re-deriving them from annotations.
--
-- Only `reactant`, `product` and `enzyme` are projected. `cofactor`,
-- `regulator` and `member` are counted and left out: the binarisation COSMOS
-- was built from has no such roles, every species in it is a reactant or a
-- product by the sign of its coefficient, and projecting a regulator as a
-- substrate would state something no resource said.
--
-- A reaction the resources call reversible is emitted twice: once as written,
-- and once with the arrow turned round, so that the products feed the enzyme
-- and the reactants come out of it. The gene node of the second half carries
-- the `_rev` suffix, which is what makes it a node of its own; the metabolite
-- labels are untouched, because a molecule is the same molecule whichever way
-- the reaction runs. Both halves keep the reaction's index.
--
-- A reaction nobody stated a direction for is emitted forward only, and
-- `direction` stays null rather than being filled in. Rhea publishes no
-- direction on any of its reactions, so the null block is large, and it means
-- "no resource said" — which is not the same claim as "runs one way".

INSERT INTO cosmos_edge (
  build_id,
  source_label, target_label,
  source_entity_id, target_entity_id,
  source_type, target_type,
  source_compartment, target_compartment,
  mor, interaction_type,
  reaction_entity_id, interaction_id, reaction_index,
  orphan, reverse, direction,
  source_id_type, target_id_type,
  sources
)
WITH reaction_event AS (
  -- The event entity behind a header, taken from the record rather than from
  -- the projection's staging tables: those are unlogged scratch and the derive
  -- drops them the moment the interactions phase ends, so nothing outside that
  -- phase can read them. The record survives, and a reaction's record rows are
  -- the star's own `(event, member)` pairs, whose subject is the event.
  --
  -- A header can carry more than one event: it is keyed on the chemistry, so
  -- two resources describing the same reaction through two event entities
  -- reach one header. `min` picks a representative, which is what keeps the
  -- reaction index stable across rebuilds.
  --
  -- The pick is an ordered aggregate rather than `min`, which Postgres does
  -- not define over `uuid`, and it orders on the same expression the reaction
  -- index ranks on, so the representative and the numbering cannot disagree.
  SELECT
    record.interaction_id,
    (array_agg(
      record.subject_entity_id ORDER BY record.subject_entity_id
    ))[1] AS reaction_entity_id
  FROM interaction_fact_resource record
  JOIN entity event ON event.entity_id = record.subject_entity_id
  JOIN vocab_entity_type event_type
    ON event_type.entity_type_id = event.entity_type_id
  WHERE event_type.name = ANY(%(reaction_entity_types)s)
  GROUP BY record.interaction_id
),
reaction_header AS (
  -- The role gate is not decoration. An event entity is the subject of plenty
  -- of ordinary pairs too — a reaction is part of a pathway — and those reach
  -- an arity-2 header that carries `subject`/`object` parties and no chemistry.
  -- Requiring a side role is what tells the two apart.
  SELECT
    header.interaction_id,
    header.sources,
    event.reaction_entity_id
  FROM interaction header
  JOIN reaction_event event ON event.interaction_id = header.interaction_id
  -- The resource gate keeps this a metabolic network. A general reaction
  -- database states its events in proteins and complexes rather than small
  -- molecules, and one of those reaching the output would be labelled a
  -- metabolite by the role it plays rather than by what it is. The overlap
  -- test rather than a containment one is deliberate: a reaction two
  -- resources describe stays on the word of whichever of them is in scope.
  WHERE header.sources && %(reaction_sources)s
    AND EXISTS (
    SELECT 1
    FROM interaction_party party
    JOIN vocab_relation_role role
      ON role.relation_role_id = party.role_id
    WHERE party.interaction_id = header.interaction_id
      AND role.name IN ('reactant', 'product')
  )
),
reaction_direction AS (
  -- Reversibility is an annotation of the **event entity**, not of any one of
  -- its participant relations, so it is read off the entity's evidence rather
  -- than off the relation evidence the roles came from.
  --
  -- Two resources can disagree, and 539 reactions on the current build do.
  -- Reversible wins: stating that a reaction runs left to right is often just
  -- the direction a resource writes its equations in, while stating that it is
  -- reversible is a positive claim about the reverse half, and dropping that
  -- claim would silently delete edges a resource published.
  --
  -- The two spellings of each value — `REVERSIBLE` against `reversible`,
  -- `LEFT-TO-RIGHT` against `left_to_right` — come from the resources
  -- themselves, so the comparison folds case and reads the hyphen and the
  -- underscore alike.
  SELECT
    resolution.entity_id AS reaction_entity_id,
    CASE
      WHEN bool_or(replace(lower(direction.value), '-', '_') = 'reversible')
      THEN 'reversible'
      WHEN bool_or(replace(lower(direction.value), '-', '_') = 'left_to_right')
      THEN 'left_to_right'
    END AS direction
  FROM annotation direction
  JOIN entity_evidence_annotation stated
    ON stated.annotation_key = direction.annotation_key
  JOIN entity_evidence_resolution resolution
    ON resolution.source_id = stated.source_id
   AND resolution.entity_evidence_id = stated.entity_evidence_id
  JOIN reaction_header reaction
    ON reaction.reaction_entity_id = resolution.entity_id
  WHERE direction.term = %(direction_term)s
  GROUP BY resolution.entity_id
),
indexed_reaction AS (
  -- N, the number every gene label of this reaction carries. Ranked by the
  -- event id, so a rebuild that reads the same graph mints the same labels.
  -- The Python original numbers by first occurrence in a DataFrame, which
  -- depends on row order and is not stable; this is deliberately not that.
  SELECT
    reaction.interaction_id,
    reaction.sources,
    reaction.reaction_entity_id,
    event.canonical_identifier AS reaction_identifier,
    -- The pseudo-enzyme's label is the reaction's own identifier, which no
    -- translation touches: a reaction is neither a protein nor a metabolite,
    -- and neither namespace has anything to say about it. Its namespace is
    -- still reported in the same vocabulary as every other endpoint's.
    coalesce(event_namespace.namespace, 'unknown') AS reaction_id_type,
    -- The event's own type decides the edge kind. The header's interaction
    -- class would be the other candidate and is the worse one: the class
    -- derivation resolves every reaction to `other` on the current build, so
    -- reading it would type every transport event as a catalysis.
    CASE
      WHEN event_type.name = %(transport_entity_type)s THEN 'transport'
      ELSE 'catalysis'
    END AS interaction_type,
    stated.direction,
    dense_rank() OVER (
      ORDER BY reaction.reaction_entity_id
    )::integer AS reaction_index
  FROM reaction_header reaction
  LEFT JOIN reaction_direction stated
    ON stated.reaction_entity_id = reaction.reaction_entity_id
  JOIN entity event ON event.entity_id = reaction.reaction_entity_id
  JOIN vocab_entity_type event_type
    ON event_type.entity_type_id = event.entity_type_id
  LEFT JOIN _cos_namespace event_namespace
    ON event_namespace.identifier_type_id
       = event.canonical_identifier_type_id
),
party AS (
  -- `_cos_label` is the translation stage's answer for this entity: the
  -- identifier the label carries and the namespace that identifier belongs to.
  -- It is a left join and both columns are coalesced, because the stage is
  -- allowed to say nothing — a build with no reachable mapping database stages
  -- a fallback for every entity, and a build whose entity reached no party row
  -- of a projected role is not in the table at all.
  SELECT
    reaction.interaction_id,
    role.name AS role_name,
    party.entity_id,
    nullif(party.compartment, '') AS compartment,
    coalesce(label.identifier, participant.canonical_identifier)
      AS canonical_identifier,
    coalesce(label.id_namespace, participant_namespace.namespace, 'unknown')
      AS id_type
  FROM indexed_reaction reaction
  JOIN interaction_party party
    ON party.interaction_id = reaction.interaction_id
  JOIN vocab_relation_role role
    ON role.relation_role_id = party.role_id
  JOIN entity participant ON participant.entity_id = party.entity_id
  LEFT JOIN _cos_label label ON label.entity_id = party.entity_id
  LEFT JOIN _cos_namespace participant_namespace
    ON participant_namespace.identifier_type_id
       = participant.canonical_identifier_type_id
  WHERE role.name IN ('reactant', 'product', 'enzyme')
),
metabolite AS (
  SELECT
    party.interaction_id,
    party.role_name,
    party.entity_id,
    party.compartment,
    party.id_type,
    'Metab__' || party.canonical_identifier
      || coalesce('_' || party.compartment, '') AS label
  FROM party
  WHERE party.role_name IN ('reactant', 'product')
),
gene AS (
  SELECT
    party.interaction_id,
    party.entity_id,
    'Gene' || reaction.reaction_index || '__' || party.canonical_identifier
      AS label,
    party.compartment,
    party.id_type,
    'protein'::text AS node_type,
    false AS orphan
  FROM party
  JOIN indexed_reaction reaction
    ON reaction.interaction_id = party.interaction_id
  WHERE party.role_name = 'enzyme'
  UNION ALL
  -- A reaction nobody named a catalyst for keeps its edges, with the event
  -- standing in for the enzyme. Dropping it would throw away chemistry the
  -- resources did state in order to hide one thing they did not.
  SELECT
    reaction.interaction_id,
    NULL::uuid,
    'Gene' || reaction.reaction_index || '__orphanReac'
      || reaction.reaction_identifier,
    NULL::text,
    reaction.reaction_id_type,
    'reaction',
    true
  FROM indexed_reaction reaction
  WHERE NOT EXISTS (
    SELECT 1
    FROM party
    WHERE party.interaction_id = reaction.interaction_id
      AND party.role_name = 'enzyme'
  )
)
SELECT
  %(build_id)s,
  edge.source_label,
  edge.target_label,
  edge.source_entity_id,
  edge.target_entity_id,
  edge.source_type,
  edge.target_type,
  edge.source_compartment,
  edge.target_compartment,
  1,
  reaction.interaction_type,
  reaction.reaction_entity_id,
  reaction.interaction_id,
  reaction.reaction_index,
  gene.orphan,
  half.reverse,
  reaction.direction,
  edge.source_id_type,
  edge.target_id_type,
  coalesce(reaction.sources, ARRAY[]::text[])
FROM indexed_reaction reaction
JOIN metabolite ON metabolite.interaction_id = reaction.interaction_id
JOIN gene ON gene.interaction_id = reaction.interaction_id
-- One half for an irreversible reaction, two for a reversible one. The second
-- half is a distinct node of the network, not a second edge on the same one,
-- which is why the suffix lands on the label here rather than on the row.
CROSS JOIN LATERAL (
  SELECT
    pass.reverse,
    gene.label || CASE WHEN pass.reverse THEN '_rev' ELSE '' END AS label
  FROM (VALUES (false), (true)) AS pass (reverse)
  WHERE NOT pass.reverse OR reaction.direction = 'reversible'
) AS half
-- The metabolite's role is the whole of the orientation: in the forward half
-- what is consumed points at the enzyme and what is produced is pointed at by
-- it, and the reverse half turns both round. Laying both orderings out and
-- keeping one is what stops eighteen columns turning into eighteen CASE
-- expressions.
CROSS JOIN LATERAL (
  SELECT *
  FROM (
    VALUES
      (
        true,
        metabolite.label, metabolite.entity_id, 'metabolite'::text,
        metabolite.compartment, metabolite.id_type,
        half.label, gene.entity_id, gene.node_type,
        gene.compartment, gene.id_type
      ),
      (
        false,
        half.label, gene.entity_id, gene.node_type,
        gene.compartment, gene.id_type,
        metabolite.label, metabolite.entity_id, 'metabolite'::text,
        metabolite.compartment, metabolite.id_type
      )
  ) AS ordering (
    from_metabolite,
    source_label, source_entity_id, source_type,
    source_compartment, source_id_type,
    target_label, target_entity_id, target_type,
    target_compartment, target_id_type
  )
  WHERE ordering.from_metabolite
        = ((metabolite.role_name = 'reactant') <> half.reverse)
) AS edge
