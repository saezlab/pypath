export interface SelectedEntity {
  id: string;
  entityId?: string | number;
  entityPk?: string | number;
  name: string;
  type?: string;
  cv_terms?: string[];
  references?: string[];
  associated_entity_ids?: Array<string | number>;
  fullResult?: unknown;
  groupKey?: string;
}
export interface SelectedAnnotation {
  id: string;
  label: string;
  namespace?: string;
  definition?: string | null;
}

const EXACT_KEY =
  /^(?:[0-9a-f]{64}|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$/i;
export function entityPrimaryKeyFromSelection(entity: SelectedEntity): string | number | null {
  if (entity.entityPk != null) return entity.entityPk;
  if (typeof entity.entityId === 'number') return entity.entityId;
  const value = String(entity.entityId ?? entity.id).trim();
  return EXACT_KEY.test(value) ? value : null;
}

export function selectionFromUrl(
  url: URL | null,
  fallbackEntities: string[],
  fallbackAnnotations: string[],
) {
  const ownsSelection =
    !!url && ['selection', 'entities', 'annotations'].some((key) => url.searchParams.has(key));
  const parse = (value: string | null) => [
    ...new Set(
      (value || '')
        .split(',')
        .map((id) => id.trim())
        .filter(Boolean),
    ),
  ];
  return {
    entityIds: ownsSelection ? parse(url!.searchParams.get('entities')) : fallbackEntities,
    annotationIds: ownsSelection
      ? parse(url!.searchParams.get('annotations'))
      : fallbackAnnotations,
  };
}

export function isEntityCache(value: unknown): value is Record<string, SelectedEntity> {
  return (
    !!value &&
    typeof value === 'object' &&
    !Array.isArray(value) &&
    Object.entries(value).every(
      ([id, item]) =>
        !!item &&
        typeof item === 'object' &&
        item.id === id &&
        typeof item.name === 'string' &&
        (item.entityPk == null ||
          typeof item.entityPk === 'string' ||
          typeof item.entityPk === 'number'),
    )
  );
}
export function isAnnotationCache(value: unknown): value is Record<string, SelectedAnnotation> {
  return (
    !!value &&
    typeof value === 'object' &&
    !Array.isArray(value) &&
    Object.entries(value).every(
      ([id, item]) =>
        !!item && typeof item === 'object' && item.id === id && typeof item.label === 'string',
    )
  );
}
