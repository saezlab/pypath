import type {
  WireEntity,
  WireEntitySearch,
  WireRelation,
  WireRelationSearch,
  WireOntologyTerm,
} from './contracts';
import type { EntityWithIdentifiers } from '../types/entities';
import type { EntityRelation, EntityOntologyHierarchy } from '../domain/entity-types';
import type { InteractionListRow } from '../types/interactions';

const record = (value: unknown): value is Record<string, unknown> =>
  !!value && typeof value === 'object' && !Array.isArray(value);
const text = (value: unknown): string | undefined =>
  typeof value === 'string' ? value : undefined;

/** Normalize optional transport values once, before presentation or UI enrichment. */
export function entityFromWire(
  entity: Pick<WireEntity, 'entityPk'> & Partial<WireEntity>,
): EntityWithIdentifiers {
  if (typeof entity.entityPk !== 'string' || !entity.entityPk)
    throw new Error('API entity is missing its exact key');
  const hierarchy = entity.ontologyHierarchy;
  return {
    ...entity,
    entityPk: entity.entityPk,
    canonicalIdentifier: entity.canonicalIdentifier ?? '',
    canonicalIdentifierType: entity.canonicalIdentifierType ?? '',
    entityType: entity.entityType ?? null,
    taxonomyId: entity.taxonomyId ?? null,
    sources: entity.sources ?? [],
    entityAttributes: entity.entityAttributes ?? null,
    entityFacetHints: entity.entityFacetHints
      ? {
          chemicalClasses: entity.entityFacetHints.chemicalClasses ?? [],
          metabolicDomains: entity.entityFacetHints.metabolicDomains ?? [],
          structuralSpecificities: entity.entityFacetHints.structuralSpecificities ?? [],
        }
      : undefined,
    displayName: text(entity.displayName),
    relationCount: typeof entity.relationCount === 'number' ? entity.relationCount : undefined,
    identifiers: (entity.identifiers ?? []).map((identifier) => ({
      ...identifier,
      id: identifier.id ?? undefined,
      entityPk: identifier.entityPk ?? entity.entityPk,
      identifier: identifier.identifier ?? identifier.id ?? '',
      identifierType: identifier.identifierType ?? identifier.type ?? '',
    })),
    ontologyHierarchy: record(hierarchy)
      ? ({
          termId: text(hierarchy.termId) ?? entity.canonicalIdentifier ?? entity.entityPk,
          ontologyPrefix: text(hierarchy.ontologyPrefix) ?? null,
          ontologyId: text(hierarchy.ontologyId) ?? null,
          label: text(hierarchy.label) ?? null,
          definition: text(hierarchy.definition) ?? null,
          childCount: Number(hierarchy.childCount) || 0,
          parentCount: Number(hierarchy.parentCount) || 0,
        } satisfies EntityOntologyHierarchy)
      : null,
  };
}
export function relationFromWire(relation: WireRelation): EntityRelation {
  return {
    ...relation,
    subjectEntityPk: relation.subjectEntityPk ?? '',
    objectEntityPk: relation.objectEntityPk ?? '',
    predicate: relation.predicate ?? '',
    participantTypes: relation.participantTypes ?? [],
    evidenceCount: relation.evidenceCount ?? 0,
    sources: relation.sources ?? [],
    sign: typeof relation.sign === 'number' ? relation.sign : undefined,
    isDirected: typeof relation.isDirected === 'boolean' ? relation.isDirected : undefined,
  };
}
export function entitySearchFromWire(response: WireEntitySearch) {
  return {
    ...response,
    entities: (response.entities ?? []).map(entityFromWire),
    nextCursor: response.nextCursor ?? null,
  };
}
export function relationSearchFromWire(response: WireRelationSearch) {
  const rows: InteractionListRow[] = [];
  for (const row of response.rows ?? []) {
    if (!record(row.relation) || !record(row.subjectEntity) || !record(row.objectEntity)) continue;
    if (
      typeof row.relation.relationPk !== 'string' ||
      typeof row.subjectEntity.entityPk !== 'string' ||
      typeof row.objectEntity.entityPk !== 'string'
    )
      throw new Error('API relation row contains invalid entity keys');
    rows.push({
      relation: relationFromWire(row.relation as WireRelation),
      subjectEntity: entityFromWire(row.subjectEntity as WireEntity),
      objectEntity: entityFromWire(row.objectEntity as WireEntity),
    });
  }
  return {
    ...response,
    rows,
    relations: (response.relations ?? []).map(relationFromWire),
    nextCursor: response.nextCursor ?? null,
  };
}

export function ontologyTermFromWire(term: WireOntologyTerm) {
  return {
    ...term,
    label: term.label ?? null,
    definition: term.definition ?? null,
    sources: term.sources ?? [],
    synonyms: term.synonyms ?? [],
  };
}
