import type {
  Entity,
  EntityIdentifier,
  EntityOntologyHierarchy,
  Identifier as DrizzleIdentifier,
} from '$lib/domain/entity-types';

export type Identifier = DrizzleIdentifier;
export type EntityRecord = Entity;
export type EntityWithIdentifiers = Entity & {
  identifiers: EntityIdentifier[];
  ontologyHierarchy?: EntityOntologyHierarchy | null;
  identifiersTotal?: number;
  relationCount?: number;
};
export type HydratedEntityRecord = EntityWithIdentifiers;
