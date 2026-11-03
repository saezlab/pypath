/** Vocabulary namespaces describe record kinds, independently of distributing resource. */
const namespaceKinds: Record<string, string> = {
  rhea: 'Reaction',
  kegg_reaction: 'Reaction',
  bigg_reaction: 'Reaction',
  human_gem_reaction: 'Reaction',
  'drugcentral.target_group': 'Target group',
  pfocr: 'Figure',
  chemont: 'Chemical class',
  ec: 'Enzyme class',
  uniprot_keyword: 'Keyword',
  'human_gem.gpr_and': 'Required gene group',
  'recon3d.gpr_and': 'Required gene group',
  kegg_orthology: 'Orthology group',
  macdb_trait: 'Study trait / context',
};
export function entityKind(entity: {
  entityType?: string | null;
  canonicalIdentifierType?: string | null;
}): string | undefined {
  const type = entity.entityType?.replace(/^biolink:/i, '').toLowerCase();
  const ns = entity.canonicalIdentifierType?.toLowerCase() || '';
  if (type === 'molecular_activity' && ns === 'reactome') return 'Reaction';
  if (namespaceKinds[ns]) return namespaceKinds[ns];
  if (type === 'organism_taxon') return 'Organism';
  if (type === 'named_thing') return 'Unspecified entity';
  if (type === 'information_content_entity') return 'Information record';
  if (type === 'ontology_class') return 'Ontology term';
  return undefined;
}
