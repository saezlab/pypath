import { API_SERVICE_URL } from '$lib/server/api-config';

import type { components } from '$lib/api/types.generated';
import { resourceFromWire, type ResourceRecord } from '$lib/resources/types';
export type {
  ResourceRecord,
  ResourceFile,
  ResourceTag,
  LicenseUseStatus,
} from '$lib/resources/types';

export interface ResourcesSummary {
  totalResources: number;
  totalEntities: number;
  totalInteractions: number;
  totalBytes: number;
}

export async function listResources(release = 'latest'): Promise<ResourceRecord[]> {
  try {
    const base = API_SERVICE_URL.endsWith('/') ? API_SERVICE_URL : `${API_SERVICE_URL}/`;
    const res = await fetch(
      new URL(`resources?shape=svelte&release=${encodeURIComponent(release)}`, base),
    );
    if (!res.ok) throw new Error(`Resource catalog returned HTTP ${res.status}`);
    const data = (await res.json()) as components['schemas']['ResourcesResponse'];
    return (data.resources || []).map(resourceFromWire);
  } catch (error) {
    console.error('Resource catalog unavailable', error);
    return [];
  }
}

export function summarizeResources(resources: ResourceRecord[]): ResourcesSummary {
  return {
    totalResources: resources.length,
    totalEntities: resources.reduce((sum, resource) => sum + (resource.entity_count || 0), 0),
    totalInteractions: resources.reduce(
      (sum, resource) => sum + (resource.interaction_count || 0),
      0,
    ),
    totalBytes: resources.reduce((sum, resource) => sum + (resource.total_size_bytes || 0), 0),
  };
}
