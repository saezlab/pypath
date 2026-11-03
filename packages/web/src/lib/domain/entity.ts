import type { Entity } from './entity-types';

export function publicId(
  entity: Entity | { canonicalIdentifier?: string | null; entityPk?: string | null },
): string {
  return entity.entityPk || '';
}

export function displayLabel(entity: {
  label?: string | null;
  canonicalIdentifier?: string | null;
  entityPk?: string | null;
}): string {
  return entity.label || entity.canonicalIdentifier || entity.entityPk || '';
}
