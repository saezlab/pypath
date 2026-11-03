import { releaseFetch } from '$lib/api/release';
import type { SearchFilters } from '$lib/types/search';
import type { InteractionDetailsData } from '$lib/types/interactions';
import type { components } from './types.generated';
import {
  entityFromWire,
  entitySearchFromWire,
  relationSearchFromWire,
  ontologyTermFromWire,
} from './adapters';
import type {
  Declared,
  WireEntities,
  WireEntity,
  WireEntitySearch,
  WireRelationSearch,
  WireFacet,
  WireFacetOption,
  WireCursor,
  WireOntologyTerm,
} from './contracts';

export type SelectionScopeResponse = components['schemas']['SelectionScopeResponse'];
export type EntityDetailsResponse = components['schemas']['EntityDetailsResponse'];
export type EntityFilterOptions = components['schemas']['EntityFilterOptions'];
export type RelationFilterOptions = components['schemas']['RelationFilterOptions'];
export type ResourcesResponse = components['schemas']['ResourcesResponse'];
export type StatsSource = components['schemas']['StatsSource'];
export type StatsEntityType = components['schemas']['StatsEntityType'];
export type StatsInteractionType = components['schemas']['StatsInteractionType'];
export type TermsResponse = components['schemas']['TermsResponse'];
export type BuildManifest = components['schemas']['BuildManifest'];

export type EntitySearchCursor = WireCursor;
export type SelectionScopeMode = 'union' | 'intersection';
export type SelectionScopeRequest = Partial<
  Declared<components['schemas']['SelectionScopeRequest']>
>;

export async function fetchSelectionScope(
  params: SelectionScopeRequest,
): Promise<SelectionScopeResponse> {
  const res = await releaseFetch('/app-api/selection/scope', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!res.ok) throw new Error('Failed to resolve selection scope');
  return res.json() as Promise<SelectionScopeResponse>;
}

/** Landing examples: entities, and the group cards the grouped view shows instead. */
export async function fetchEntityExamples(signal?: AbortSignal) {
  const response = await releaseFetch('/app-api/entities/examples', { signal });
  if (!response.ok) throw new Error('Failed to load examples');
  const body = (await response.json()) as WireEntitySearch & {
    groups?: Array<{ entity: Pick<WireEntity, 'entityPk'> & Partial<WireEntity> }>;
  };
  return {
    ...entitySearchFromWire(body),
    groups: (body.groups ?? []).map((group) => ({
      ...group,
      entity: entityFromWire(group.entity),
    })),
  };
}

export type TopHit = { key: string; label: string; detail: string };

/**
 * The top hit of an entity search, which a relations search shows the relations of: the
 * top group (its group key, e.g. 'gene:entrez:7157') when grouping, else the top entity.
 */
export async function fetchTopHit(
  query: string,
  filters: SearchFilters,
  grouped: boolean,
  signal?: AbortSignal,
): Promise<TopHit | null> {
  if (!grouped) {
    const entity = (await fetchEntitiesSearch({ query, limit: 1, filters }, signal)).entities[0];
    return entity
      ? {
          key: entity.entityPk,
          label: entity.displayName || entity.canonicalIdentifier || entity.entityPk,
          detail: entity.entityType ?? '',
        }
      : null;
  }
  const response = await releaseFetch('/app-api/entities/groups', {
    method: 'POST',
    signal,
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ strategy: 'auto', query, filters, limit: 1, member_limit: 1 }),
  });
  if (!response.ok) throw new Error('Failed to resolve the search');
  const group = (
    (await response.json()) as {
      groups: Array<{
        group_key: string;
        is_group: boolean;
        member_count: number;
        entity: { displayName?: string | null; entityType?: string | null };
      }>;
    }
  ).groups[0];
  if (!group) return null;
  return {
    key: group.group_key,
    label: group.entity.displayName || group.group_key,
    detail: group.is_group
      ? `${group.group_key.startsWith('gene:') ? 'Gene group' : 'Structure group'} · ${group.member_count} entities`
      : (group.entity.entityType ?? ''),
  };
}

export async function fetchEntitiesSearch(
  params: {
    query?: string;
    limit?: number;
    cursor?: EntitySearchCursor | null;
    filters?: SearchFilters;
  },
  signal?: AbortSignal,
) {
  const hasLargeFilters =
    params.filters &&
    ((params.filters.entity_pks?.length ?? 0) > 50 ||
      (params.filters.annotation_term_ids?.length ?? 0) > 10);

  if (hasLargeFilters) {
    const res = await releaseFetch('/app-api/entities/search', {
      method: 'POST',
      signal,
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        query: params.query || '',
        limit: params.limit ?? 20,
        cursor: params.cursor,
        filters: params.filters,
      }),
    });
    if (!res.ok) throw new Error('Failed to fetch entities');
    return entitySearchFromWire((await res.json()) as WireEntitySearch);
  }

  const url = new URL('/app-api/entities/search', window.location.origin);
  if (params.query) url.searchParams.set('q', params.query);
  if (params.limit) url.searchParams.set('limit', String(params.limit));
  if (params.cursor != null) url.searchParams.set('cursor', JSON.stringify(params.cursor));
  if (params.filters && Object.keys(params.filters).length > 0) {
    url.searchParams.set('filters', JSON.stringify(params.filters));
  }

  const res = await releaseFetch(url, { signal });
  if (!res.ok) throw new Error('Failed to fetch entities');
  return entitySearchFromWire((await res.json()) as WireEntitySearch);
}

export async function fetchEntityFilterOptions() {
  const res = await releaseFetch('/app-api/entities/filter-options');
  if (!res.ok) throw new Error('Failed to fetch entity filter options');
  return res.json() as Promise<EntityFilterOptions>;
}

export async function fetchRelationsSearch(
  params: {
    filters?: SearchFilters;
    limit?: number;
    offset?: number;
  },
  signal?: AbortSignal,
) {
  const res = await releaseFetch('/app-api/relations/search', {
    method: 'POST',
    signal,
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!res.ok) throw new Error('Failed to fetch relations');
  return relationSearchFromWire((await res.json()) as WireRelationSearch);
}

export async function fetchRelationFilterOptions() {
  const res = await releaseFetch('/app-api/relations/filter-options');
  if (!res.ok) throw new Error('Failed to fetch relation filter options');
  return res.json() as Promise<RelationFilterOptions>;
}

export async function fetchRelationEvidence(relationPk: string | number, filters?: SearchFilters) {
  const query =
    filters && Object.keys(filters).length
      ? `?filters=${encodeURIComponent(JSON.stringify(filters))}`
      : '';
  const res = await releaseFetch(
    `/app-api/relations/${encodeURIComponent(String(relationPk))}/evidence${query}`,
  );
  if (!res.ok) throw new Error('Failed to fetch relation evidence');
  return res.json() as Promise<{ evidence: InteractionDetailsData['evidence'] }>;
}

export async function fetchEntitiesByPublicIds(publicIds: string[]) {
  const res = await releaseFetch('/app-api/entities/by-public-ids', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ public_ids: publicIds }),
  });
  if (!res.ok) throw new Error('Failed to fetch entities by public ids');
  const data = (await res.json()) as WireEntities;
  return { ...data, entities: (data.entities ?? []).map(entityFromWire) };
}

export async function fetchEntitiesByPks(pks: Array<string | number>) {
  const res = await releaseFetch('/app-api/entities/by-pks', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ pks }),
  });
  if (!res.ok) throw new Error('Failed to fetch entities by pks');
  const data = (await res.json()) as WireEntities;
  return { ...data, entities: (data.entities ?? []).map(entityFromWire) };
}

export async function fetchOntologySearch(
  params: {
    query?: string;
    prefixes?: string[];
    ontologyIds?: string[];
    limit?: number;
    offset?: number;
  },
  signal?: AbortSignal,
) {
  const url = new URL('/app-api/ontology/search', window.location.origin);
  if (params.query) url.searchParams.set('q', params.query);
  if (params.limit != null) url.searchParams.set('limit', String(params.limit));
  if (params.offset != null) url.searchParams.set('offset', String(params.offset));
  if (params.prefixes?.length) url.searchParams.set('prefixes', params.prefixes.join(','));
  if (params.ontologyIds?.length) url.searchParams.set('ontologyIds', params.ontologyIds.join(','));

  const res = await releaseFetch(url, { signal });
  if (!res.ok) throw new Error('Failed to fetch ontology terms');
  return ((await res.json()) as WireOntologyTerm[]).map(ontologyTermFromWire);
}

export async function fetchScopedOntologySearch(
  params: {
    entityPks?: Array<string | number>;
    termIds?: string[];
    selectionScope?: SelectionScopeRequest;
    query?: string;
    prefixes?: string[];
    ontologyIds?: string[];
    limit?: number;
    offset?: number;
  },
  signal?: AbortSignal,
) {
  const res = await releaseFetch('/app-api/ontology/scoped-search', {
    method: 'POST',
    signal,
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      entityPks: params.entityPks || [],
      termIds: params.termIds || [],
      selectionScope: params.selectionScope,
      query: params.query || '',
      prefixes: params.prefixes,
      ontologyIds: params.ontologyIds,
      limit: params.limit ?? 24,
      offset: params.offset ?? 0,
    }),
  });
  if (!res.ok) throw new Error('Failed to fetch scoped ontology terms');
  return ((await res.json()) as WireOntologyTerm[]).map(ontologyTermFromWire);
}

export async function fetchOntologyPrefixes() {
  const res = await releaseFetch('/app-api/ontology/prefixes');
  if (!res.ok) throw new Error('Failed to fetch ontology prefixes');
  return res.json() as Promise<{ prefixes: string[] }>;
}

export async function fetchScopedOntologyPrefixCounts(params: {
  entityPks?: Array<string | number>;
  annotationTermIds?: string[];
  query?: string;
}) {
  const res = await releaseFetch('/app-api/ontology/prefix-counts', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!res.ok) throw new Error('Failed to fetch scoped ontology prefix counts');
  return res.json() as Promise<Array<{ prefix: string; scopedCount: number }>>;
}

export async function fetchScopedOntologyIdCounts(
  params: {
    entityPks?: Array<string | number>;
    annotationTermIds?: string[];
    selectionScope?: SelectionScopeRequest;
    query?: string;
  },
  signal?: AbortSignal,
) {
  const res = await releaseFetch('/app-api/ontology/ontology-id-counts', {
    method: 'POST',
    signal,
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!res.ok) throw new Error('Failed to fetch scoped ontology id counts');
  return res.json() as Promise<Array<{ ontologyId: string; scopedCount: number }>>;
}

export async function fetchTerms(termIds: string[]) {
  const res = await releaseFetch('/app-api/terms', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ term_ids: termIds }),
  });
  if (!res.ok) throw new Error('Failed to fetch terms');
  return res.json() as Promise<{ terms: Record<string, unknown | null> }>;
}

export async function fetchScopedEntityFacetCounts(
  params: {
    entityIds?: Array<string | number>;
    annotationTermIds?: string[];
    entityTypes?: string[];
    sources?: string[];
    ncbi_tax_id?: string[];
    query?: string;
    facetLimit?: number;
  },
  signal?: AbortSignal,
) {
  const res = await releaseFetch('/app-api/entities/scoped-facets', {
    method: 'POST',
    signal,
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!res.ok) throw new Error('Failed to fetch scoped entity facet counts');
  return res.json() as Promise<WireFacet[]>;
}

/** The relation filter options the scope and the other groups' filters leave. */
export async function fetchScopedRelationFacetOptions(
  params: {
    entityIds?: Array<string | number>;
    endpointMode?: 'any' | 'both';
    mode?: SelectionScopeMode;
    annotationTermIds?: string[];
    predicates?: string[];
    relation_categories?: string[];
    object_aspect_qualifier?: string[];
    object_direction_qualifier?: string[];
    causal_mechanism_qualifier?: string[];
    interactionTypes?: string[];
    partnerTypes?: string[];
    sources?: string[];
    taxonomyIds?: string[];
    taxonomyLimit?: number;
    taxonomyQuery?: string;
  },
  signal?: AbortSignal,
) {
  const res = await releaseFetch('/app-api/relations/scoped-facets', {
    method: 'POST',
    signal,
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  });
  if (!res.ok) throw new Error('Failed to fetch relation filter options');
  return res.json() as Promise<WireFacetOption[]>;
}

// Statistics transport aliases come from generated OpenAPI schemas.

async function fetchStats<T>(path: string): Promise<T> {
  const res = await releaseFetch(`/api/stats/${path}`);
  if (!res.ok) throw new Error(`Failed to fetch stats/${path}`);
  return res.json() as Promise<T>;
}

export const fetchStatsSources = () => fetchStats<StatsSource[]>('sources');
export const fetchStatsEntityTypes = () => fetchStats<StatsEntityType[]>('entity-types');
export const fetchStatsInteractionTypes = () =>
  fetchStats<StatsInteractionType[]>('interaction-types');
export const fetchStatsIdentifierTypes = () =>
  fetchStats<Array<{ identifierType: string; count: number }>>('identifier-types');
export const fetchStatsChemicalClasses = () =>
  fetchStats<Array<{ chemicalClass: string; count: number }>>('chemical-classes');
export const fetchStatsMetabolicDomains = () =>
  fetchStats<Array<{ metabolicDomain: string; count: number }>>('metabolic-domains');
export const fetchStatsBuildManifest = () => fetchStats<BuildManifest>('build-manifest');
export const fetchStatsCoverageProfile = () =>
  fetchStats<Array<{ nResources: number; nEntities: number }>>('coverage-profile');
export const fetchStatsResourceOverlap = (contentKind = 'entity') =>
  fetchStats<Array<{ sourceA: string; sourceB: string; overlap: number }>>(
    `resource-overlap?contentKind=${encodeURIComponent(contentKind)}`,
  );
