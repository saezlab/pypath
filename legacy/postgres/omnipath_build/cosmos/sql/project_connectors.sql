-- The connectors: one edge per distinct gene node, from the bare identifier to
-- the node that carries it.
--
-- A gene node is an enzyme in one reaction, so its label is unusable as a join
-- key by anything upstream: transcriptomics arrives keyed on a gene, not on a
-- gene-in-a-reaction. The connector is the attachment point, and it is what
-- lets a measurement reach every instance of its enzyme at once.
--
-- The bare side of an orphan node is the reaction's own identifier, not the
-- `orphanReac` spelling of it, because that is the identifier a caller holds.
-- The reverse half of a reversible reaction is a node of its own, so it gets
-- its own connector from the same bare identifier.
--
-- Metabolite nodes get no connector. Their labels are the identifier plus a
-- compartment suffix, which a consumer can strip; a gene label cannot be
-- stripped without knowing the reaction index.

INSERT INTO cosmos_edge (
  build_id,
  source_label, target_label,
  source_entity_id, target_entity_id,
  source_type, target_type,
  mor, interaction_type,
  reaction_entity_id, interaction_id, reaction_index,
  orphan, reverse, direction,
  source_id_type, target_id_type,
  sources
)
WITH gene_node AS (
  SELECT DISTINCT
    node.label,
    node.entity_id,
    node.id_type,
    edge.reaction_entity_id,
    edge.interaction_id,
    edge.reaction_index,
    edge.orphan,
    edge.reverse,
    edge.direction,
    edge.sources
  FROM cosmos_edge edge
  CROSS JOIN LATERAL (
    VALUES
      (edge.source_label, edge.source_entity_id, edge.source_type,
       edge.source_id_type),
      (edge.target_label, edge.target_entity_id, edge.target_type,
       edge.target_id_type)
  ) AS node (label, entity_id, node_type, id_type)
  WHERE edge.interaction_type <> 'connector'
    AND node.node_type <> 'metabolite'
)
SELECT
  %(build_id)s,
  -- The bare side is the same identifier the node label carries, translation
  -- and all. Reading the entity's canonical identifier here instead would
  -- point a connector from an Entrez gene at a node named after a UniProt
  -- accession, and the connector exists precisely so that a caller holding one
  -- identifier reaches every node that carries it.
  coalesce(
    enzyme_label.identifier,
    enzyme.canonical_identifier,
    event.canonical_identifier
  ),
  gene_node.label,
  gene_node.entity_id,
  gene_node.entity_id,
  -- Both ends stand in the gene position of the network, so both are typed as
  -- the gene side. The pseudo-enzyme's own `reaction` type is recorded on the
  -- reaction edges, where it describes an endpoint of the chemistry rather
  -- than the identity of one node stated twice.
  'protein',
  'protein',
  1,
  'connector',
  gene_node.reaction_entity_id,
  gene_node.interaction_id,
  gene_node.reaction_index,
  gene_node.orphan,
  gene_node.reverse,
  gene_node.direction,
  gene_node.id_type,
  gene_node.id_type,
  gene_node.sources
FROM gene_node
LEFT JOIN entity enzyme ON enzyme.entity_id = gene_node.entity_id
LEFT JOIN _cos_label enzyme_label ON enzyme_label.entity_id = gene_node.entity_id
LEFT JOIN entity event ON event.entity_id = gene_node.reaction_entity_id
WHERE coalesce(
        enzyme_label.identifier,
        enzyme.canonical_identifier,
        event.canonical_identifier
      ) IS NOT NULL
