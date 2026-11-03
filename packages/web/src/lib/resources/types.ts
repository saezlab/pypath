import type { components } from '../api/types.generated';
import type { Declared, RequiredValues } from '../api/contracts';
import { BIOLINK_ENTITIES, ENTITY_LOOKUP } from '../domain/biolink.generated';

export type ResourceFile = Declared<components['schemas']['ResourceFile']>;
type DeclaredTag = Declared<components['schemas']['ResourceTag']>;
/** A tag the input module declares, or an entity type its build holds. */
export type ResourceTag = Omit<DeclaredTag, 'dimension'> & {
  dimension: DeclaredTag['dimension'] | 'entity_type';
};
export type LicenseUseStatus = components['schemas']['ResourceLicenseUse']['academic'];
type WireResource = components['schemas']['ResourceRecord'];
export type ResourceRecord = Partial<
  Omit<
    Declared<WireResource>,
    | 'resource_id'
    | 'resource_name'
    | 'version'
    | 'entity_count'
    | 'interaction_count'
    | 'total_size_bytes'
    | 'files'
    | 'resolution_stats'
    | 'tags'
  >
> &
  RequiredValues<
    WireResource,
    'resource_id' | 'resource_name' | 'entity_count' | 'interaction_count' | 'total_size_bytes'
  > & {
    version: string | null;
    files: ResourceFile[];
    tags?: ResourceTag[];
    resolution_stats?: {
      input_entities: number;
      resolved_entities: number;
      unresolved_entities: number;
      not_applicable_entities: number;
    } | null;
  };
/** The entity types a build holds, most frequent first: not declared, counted by the build. */
function entityTypeTags(counts: unknown): ResourceTag[] {
  return Object.entries((counts ?? {}) as Record<string, { input?: number }>)
    .sort(([, a], [, b]) => (b.input ?? 0) - (a.input ?? 0))
    .map(([type]) => {
      const term = BIOLINK_ENTITIES[ENTITY_LOOKUP[type.toLowerCase()]];
      return {
        id: `entity:${type}`,
        label: term?.label ?? type,
        dimension: 'entity_type',
        description: term?.description ?? '',
      };
    });
}
export function resourceFromWire(resource: WireResource): ResourceRecord {
  const id = resource.resource_id || resource.key || resource.resource;
  if (!id) throw new Error('Resource catalog entry has no identity');
  const stats = resource.resolution_stats;
  return {
    ...resource,
    tags: [...(resource.tags ?? []), ...entityTypeTags(stats?.by_entity_type)],
    resource_id: id,
    // Curated display names ("BindingDB") replace bare identifiers ("bindingdb").
    resource_name: resource.display_name || resource.resource_name || resource.resource || id,
    version: resource.version ?? null,
    files: resource.files ?? [],
    entity_count: resource.entity_count ?? resource.entities_count ?? 0,
    interaction_count: resource.interaction_count ?? resource.relations_count ?? 0,
    total_size_bytes: resource.total_size_bytes ?? 0,
    resolution_stats: stats
      ? {
          input_entities: Number(stats.input_entities) || 0,
          resolved_entities: Number(stats.resolved_entities) || 0,
          unresolved_entities: Number(stats.unresolved_entities) || 0,
          not_applicable_entities: Number(stats.not_applicable_entities) || 0,
        }
      : null,
  };
}
