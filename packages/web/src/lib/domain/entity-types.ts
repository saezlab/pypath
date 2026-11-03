import type {
  Declared,
  RequiredValues,
  WireEntity,
  WireIdentifier,
  WireRelation,
} from '../api/contracts';

/**
 * Domain entity and relation types for the OmniPath web application.
 * Replaces the legacy Drizzle schema type definitions.
 */

export type EntityFacetHints = {
  chemicalClasses: string[];
  metabolicDomains: string[];
  structuralSpecificities: string[];
};

/** Transport defaults are normalized by api/adapters; these additions are UI state. */
export type Entity = Partial<
  Omit<
    Declared<WireEntity>,
    | 'identifiers'
    | 'canonicalIdentifier'
    | 'canonicalIdentifierType'
    | 'entityType'
    | 'taxonomyId'
    | 'sources'
    | 'entityAttributes'
    | 'entityFacetHints'
    | 'relationCount'
  >
> &
  Pick<WireEntity, 'entityPk'> &
  RequiredValues<WireEntity, 'canonicalIdentifier' | 'canonicalIdentifierType' | 'sources'> & {
    displayName?: string;
    identifiersNextCursor?: string | null;
    groupStrategy?: 'chemical_connectivity' | 'gene_reference';
    memberEntityTypes?: string[];
    groupMemberKeys?: string[];
    groupMemberCount?: number;
    groupResources?: string[];
    groupDetailsLoaded?: boolean;
    groupMemberCursor?: string | null;
    detailNextCursor?: string | null;
    groupQuery?: string;
    groupFilters?: import('$lib/types/search').SearchFilters;
    label?: string | null;
    resolutionStatus?: string | null;
    entityType: string | null;
    taxonomyId: string | null;
    entityAttributes: unknown;
    sources: string[];
    relationCount?: number;
    entityFacetHints?: EntityFacetHints;
  };

export type EntityOntologyHierarchy = {
  termId: string;
  ontologyPrefix: string | null;
  label: string | null;
  definition: string | null;
  ontologyId: string | null;
  childCount: number;
  parentCount: number;
};

export type EntityInsert = Entity;

export type EntityIdentifier = Omit<
  Declared<WireIdentifier>,
  'entityPk' | 'identifier' | 'identifierType' | 'id'
> &
  RequiredValues<WireIdentifier, 'entityPk' | 'identifier' | 'identifierType'> & { id?: string };

export type Identifier = {
  key: string;
  value: string;
  id?: string;
  entityPk?: string;
  identifier?: string;
  identifierType?: string;
};

export type EntityRelation = Omit<
  Declared<WireRelation>,
  | 'subjectEntityPk'
  | 'objectEntityPk'
  | 'predicate'
  | 'participantTypes'
  | 'evidenceCount'
  | 'sources'
  | 'sign'
  | 'isDirected'
> &
  RequiredValues<
    WireRelation,
    | 'subjectEntityPk'
    | 'objectEntityPk'
    | 'predicate'
    | 'participantTypes'
    | 'evidenceCount'
    | 'sources'
  > & {
    sign?: number;
    isDirected?: boolean;
  };

export type EntityRelationEvidence = {
  source: string;
  relationEvidencePk: string;
  relationPk: string;
  recordAttributes: unknown;
  subjectAttributes: unknown;
  objectAttributes: unknown;
  evidence: unknown;
};
