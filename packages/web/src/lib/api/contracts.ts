import type { components } from './types.generated';
/** Generated transport fields, excluding OpenAPI's arbitrary-extra-field index. */
export type Declared<T> = {
  [K in keyof T as string extends K ? never : number extends K ? never : K]: T[K];
};
export type RequiredValues<T, K extends keyof T> = { [P in K]-?: NonNullable<T[P]> };
export type WireEntity = components['schemas']['EntitySummary'];
export type WireIdentifier = components['schemas']['EntityIdentifier'];
export type WireRelation = components['schemas']['RelationSummary'];
export type WireEntitySearch = components['schemas']['EntitySearchResponse'];
export type WireRelationSearch = components['schemas']['RelationSearchResponse'];
export type WireEntities = components['schemas']['EntitiesResponse'];
export type WireFilters = components['schemas']['SearchFilters'];
export type WireFacet = components['schemas']['FacetCount'];
export type WireFacetOption = components['schemas']['FacetOption'];
export type WireCursor = components['schemas']['EntitySearchCursor'];

export type WireOntologyTerm = components['schemas']['OntologyTerm'];
export type WireOntologyChildren = components['schemas']['OntologyChildrenResponse'];
export type WireAdminJob = components['schemas']['AdminJobResponse'];
export type WireAdminStage = components['schemas']['AdminStage'];
export type WireAdminJobs = components['schemas']['AdminJobsResponse'];
