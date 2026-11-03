import type { ResourceRecord, ResourceTag } from './types';

export const dimensions = [
  { id: 'topic', label: 'Biological topic' },
  { id: 'content', label: 'Content' },
  { id: 'knowledge_origin', label: 'Knowledge origin' },
  { id: 'entity_type', label: 'Entities' },
] as const;
export type LicenseAudience = 'academic' | 'commercial';
export type ResourceSort = 'name' | 'entities' | 'relations' | 'size';

export function resourceTags(resources: ResourceRecord[]): ResourceTag[] {
  const tags = new Map<string, ResourceTag>();
  for (const resource of resources) for (const tag of resource.tags ?? []) tags.set(tag.id, tag);
  return [...tags.values()].sort((a, b) => a.label.localeCompare(b.label));
}

/** OR within a dimension, AND across dimensions and license audiences. */
export function filterResources(
  resources: ResourceRecord[],
  query: string,
  selected: string[],
  licenses: string[],
): ResourceRecord[] {
  const vocabulary = new Map(resourceTags(resources).map((tag) => [tag.id, tag]));
  const groups = new Map<string, string[]>();
  for (const id of selected) {
    const dimension = vocabulary.get(id)?.dimension ?? 'unknown';
    groups.set(dimension, [...(groups.get(dimension) ?? []), id]);
  }
  const needle = query.trim().toLowerCase();
  return resources.filter((resource) => {
    const tags = resource.tags ?? [];
    const haystack = [
      resource.resource_id,
      resource.resource_name,
      resource.description,
      ...tags.map((tag) => tag.label),
    ]
      .join(' ')
      .toLowerCase();
    return (
      (!needle || haystack.includes(needle)) &&
      [...groups.values()].every((ids) => tags.some((tag) => ids.includes(tag.id))) &&
      licenses.every(
        (audience) =>
          (audience === 'academic' || audience === 'commercial') &&
          resource.license_use?.[audience] === 'allowed',
      )
    );
  });
}

export type SortDirection = 'asc' | 'desc';

/** Names read A–Z; counts and sizes read largest first. */
export function defaultDirection(sort: string): SortDirection {
  return sort === 'entities' || sort === 'relations' || sort === 'size' ? 'desc' : 'asc';
}

export function sortResources(
  resources: ResourceRecord[],
  sort: string,
  direction: SortDirection = defaultDirection(sort),
): ResourceRecord[] {
  const fields = {
    entities: 'entity_count',
    relations: 'interaction_count',
    size: 'total_size_bytes',
  } as const;
  const field = fields[sort as keyof typeof fields];
  const sign = direction === 'asc' ? 1 : -1;
  return [...resources].sort((a, b) => {
    const primary = field
      ? ((a[field] || 0) - (b[field] || 0)) * sign
      : a.resource_name.localeCompare(b.resource_name) * sign;
    return primary || a.resource_name.localeCompare(b.resource_name);
  });
}

export function formatBytes(bytes = 0): string {
  if (bytes <= 0) return '0 B';
  const units = ['B', 'KB', 'MB', 'GB', 'TB'];
  const index = Math.min(Math.floor(Math.log(bytes) / Math.log(1024)), units.length - 1);
  return `${(bytes / 1024 ** index).toFixed(index ? 1 : 0)} ${units[index]}`;
}
