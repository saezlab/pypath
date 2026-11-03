import type { Entity, Identifier as LegacyIdentifier } from '$lib/domain/entity-types';
import type { EntityWithIdentifiers } from '$lib/types/entities';
import { entityKind } from './../domain/presentation';
import { publicId } from '$lib/domain/entity';
import { IDENTIFIER_SLUGS } from '$lib/domain/identifier-types';
import {
  BIOLINK_PREDICATES,
  PREDICATE_LOOKUP,
  BIOLINK_ENTITIES,
  ENTITY_LOOKUP,
  BIOLINK_SLOTS,
  SLOT_LOOKUP,
} from '$lib/domain/biolink.generated';

export type EntityLike = Entity | EntityWithIdentifiers;
export type EntityIdentifierLike = LegacyIdentifier;

function isObject(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null;
}

export function isEntityRecord(entity: EntityLike): entity is Entity {
  return 'entityPk' in entity && 'canonicalIdentifier' in entity;
}

export function getEntityPublicId(entity: EntityLike): string {
  return publicId(entity);
}

export function getEntityTypeValue(entity: EntityLike): string | undefined | null {
  return entity.entityType;
}

export function getEntityTypeLabel(entity: EntityLike): string {
  const entityType = getEntityTypeValue(entity);
  return entityKind(entity) || (entityType ? getIdentifierTypeLabel(entityType) : 'Entity');
}

function isFallbackIdentifierType(identifierType: string | null | undefined): boolean {
  const rawType = identifierType?.trim().toLowerCase() || '';
  const label = getIdentifierTypeLabel(identifierType || '').toLowerCase();
  const normalized = [rawType, label].join(' ').replace(/[_-]/g, ' ');
  return label === 'fallback' || normalized.includes('unresolved entity key');
}

function isFallbackEntity(entity: EntityLike): boolean {
  return isFallbackIdentifierType(entity.canonicalIdentifierType);
}

export function getEntityResolutionStatus(entity: EntityLike): string | undefined {
  const status = (entity as { resolutionStatus?: unknown }).resolutionStatus;
  return typeof status === 'string' && status.trim() ? status.trim().toLowerCase() : undefined;
}

function normalizedEntityTypeLabel(entity: EntityLike): string {
  return (entity.entityType || '').toLowerCase().replace(/[\s_]/g, '');
}

export function isResolverBackedEntity(entity: EntityLike): boolean {
  const typeLabel = normalizedEntityTypeLabel(entity);
  return typeLabel === 'protein' || isChemicalEntity(entity);
}

export function isUnresolvedEntity(entity: EntityLike): boolean {
  if (!isResolverBackedEntity(entity)) return false;
  const status = getEntityResolutionStatus(entity);
  if (status) return status === 'unresolved';
  return isFallbackEntity(entity);
}

function isHashLikeIdentifier(value: string | null | undefined): boolean {
  const text = value?.trim() || '';
  return /^[a-f0-9-]{24,}$/i.test(text) && /[a-f]/i.test(text);
}

export function isChemicalEntity(entity: EntityLike): boolean {
  const typeLabel = normalizedEntityTypeLabel(entity);
  return (
    typeLabel === 'chemicalentity' ||
    typeLabel === 'chemical' ||
    typeLabel === 'smallmolecule' ||
    typeLabel === 'compound' ||
    typeLabel === 'metabolite' ||
    typeLabel === 'drug' ||
    typeLabel === 'lipid'
  );
}

export function isCvTermEntity(entity: EntityLike): boolean {
  const typeLabel = (entity.entityType || '').toLowerCase().replace(/[\s_]/g, '');
  return typeLabel === 'ontologyclass' || typeLabel === 'cvterm';
}

export function getEntityIdentifiers(entity: EntityLike): EntityIdentifierLike[] {
  const typedEntity = entity as { identifiers?: unknown };
  const raw = Array.isArray(typedEntity.identifiers) ? typedEntity.identifiers : [];
  return raw
    .map((item: unknown) => {
      if (!isObject(item)) return null;
      const key =
        typeof item.key === 'string'
          ? item.key
          : typeof item.identifierType === 'string'
            ? item.identifierType
            : typeof item.identifier_type === 'string'
              ? item.identifier_type
              : null;
      const value =
        typeof item.value === 'string'
          ? item.value
          : typeof item.identifier === 'string'
            ? item.identifier
            : null;

      if (!key || !value) return null;
      return { key, value } satisfies EntityIdentifierLike;
    })
    .filter((identifier): identifier is EntityIdentifierLike => Boolean(identifier));
}

function getShortestAvailableIdentifier(
  entity: EntityLike,
  excludedValues: string[] = [],
): EntityIdentifierLike | undefined {
  const excluded = new Set(excludedValues.map((value) => value.trim()).filter(Boolean));
  const identifiers = getEntityIdentifiers(entity)
    .map((identifier) => ({
      key: identifier.key.trim(),
      value: identifier.value.trim(),
    }))
    .filter((identifier) => identifier.key && identifier.value && !excluded.has(identifier.value));
  const nonFallbackIdentifiers = identifiers.filter(
    (identifier) => !isFallbackIdentifierType(identifier.key),
  );
  const candidates = nonFallbackIdentifiers.length > 0 ? nonFallbackIdentifiers : identifiers;

  return candidates.sort(
    (a, b) =>
      a.value.length - b.value.length ||
      a.value.localeCompare(b.value) ||
      a.key.localeCompare(b.key),
  )[0];
}

export function getEntityPrimaryIdentifierBadge(entity: EntityLike): EntityIdentifierLike {
  const fallbackIdentifier = isFallbackEntity(entity)
    ? getShortestAvailableIdentifier(entity)
    : undefined;
  if (fallbackIdentifier) return fallbackIdentifier;
  if (isFallbackEntity(entity)) return { key: 'Identifier', value: 'Unavailable' };
  return { key: entity.canonicalIdentifierType, value: entity.canonicalIdentifier };
}

export function getRelationPredicateLabel(predicate: string | null | undefined, sign = 0): string {
  if (predicate?.replace(/^biolink:/i, '').toLowerCase() === 'affects') {
    if (sign === 1) return 'Affects (increased)';
    if (sign === -1) return 'Affects (decreased)';
  }
  if (
    [
      'physically_interacts_with',
      'directly_physically_interacts_with',
      'indirectly_physically_interacts_with',
    ].includes(predicate?.replace(/^biolink:/i, '') || '')
  )
    return 'Interacts with';
  const text = predicate?.trim() || '';
  const key = PREDICATE_LOOKUP[text.toLowerCase()];
  return key ? BIOLINK_PREDICATES[key].label : text;
}

export function getIdentifierTypeLabel(identifierType: string | null | undefined): string {
  const text = identifierType?.trim() || '';
  const key = text.toLowerCase();
  const entity = ENTITY_LOOKUP[key];
  if (entity) return BIOLINK_ENTITIES[entity].label;
  const slot = SLOT_LOOKUP[key];
  if (slot) return BIOLINK_SLOTS[slot].label;
  return IDENTIFIER_SLUGS[key] ?? text;
}

export function getIdentifierDisplayTypeForValue(
  entity: EntityLike,
  value: string,
): string | undefined {
  const trimmed = value.trim();
  if (!trimmed) return undefined;

  const identifier = getEntityIdentifiers(entity).find((item) => item.value.trim() === trimmed);
  if (identifier) return identifier.key;

  if (entity.canonicalIdentifier?.trim() === trimmed) {
    return entity.canonicalIdentifierType;
  }

  return undefined;
}

function identifierLabel(identifierType: string): string {
  return getIdentifierTypeLabel(identifierType).toLowerCase();
}

function normalizedIdentifierTypeText(identifierType: string): string {
  return [identifierType, getIdentifierTypeLabel(identifierType)]
    .join(' ')
    .toLowerCase()
    .replace(/[_-]/g, ' ');
}

function identifierTypeMatches(identifierType: string, types: string[]): boolean {
  const normalized = normalizedIdentifierTypeText(identifierType);
  return types.some((type) => normalized.includes(type.toLowerCase().replace(/[_-]/g, ' ')));
}

function uniqueStrings(values: Array<string | null | undefined>): string[] {
  return Array.from(new Set(values.map((value) => value?.trim()).filter(Boolean) as string[]));
}

export function classifyEntityIdentifiers(entity: EntityLike): {
  names: string[];
  synonyms: string[];
  geneSymbols: string[];
} {
  const names: string[] = [];
  const synonyms: string[] = [];
  const geneSymbols: string[] = [];

  const rawLabel = (entity as { label?: unknown }).label;
  if (typeof rawLabel === 'string' && rawLabel.trim()) {
    const trimmedLabel = rawLabel.trim();
    const entityType = getEntityTypeLabel(entity).toLowerCase();
    if (entityType === 'protein' || entityType === 'gene') {
      geneSymbols.push(trimmedLabel);
    } else {
      names.push(trimmedLabel);
    }
  }

  for (const identifier of getEntityIdentifiers(entity)) {
    const rawKey = identifier.key.trim().toLowerCase();
    const label = identifierLabel(identifier.key);
    const value = identifier.value.trim();
    if (!value) continue;

    if (
      rawKey === 'genesymbol' ||
      rawKey === 'gene_symbol' ||
      rawKey === 'gene_name_primary' ||
      rawKey === 'symbol' ||
      rawKey === 'om:0200' ||
      label.includes('gene symbol') ||
      label.includes('gene name primary') ||
      label.includes('gene name')
    ) {
      geneSymbols.push(value);
      continue;
    }
    if (
      rawKey === 'genesymbol-syn' ||
      rawKey === 'gene_name_synonym' ||
      rawKey === 'gene_symbol_synonym' ||
      rawKey === 'om:0201' ||
      label.includes('gene name synonym') ||
      label.includes('gene symbol synonym')
    ) {
      synonyms.push(value);
      continue;
    }
    if (
      label === 'name' ||
      label.endsWith(':name') ||
      label.includes(' entry name') ||
      rawKey === 'name' ||
      rawKey === 'om:0202'
    ) {
      names.push(value);
      continue;
    }
    if (label.includes('synonym') || rawKey.includes('synonym') || rawKey === 'om:0203') {
      synonyms.push(value);
    }
  }

  return {
    names: uniqueStrings(names),
    synonyms: uniqueStrings(synonyms),
    geneSymbols: uniqueStrings(geneSymbols),
  };
}

function mapEntityAttributesToDescriptions(attributes: unknown): string[] {
  if (!Array.isArray(attributes)) return [];

  const preferredKeywords = [
    'function',
    'description',
    'disease',
    'subcellular location',
    'pathway',
    'activity regulation',
    'tissue specificity',
    'developmental stage',
    'note',
  ];

  const preferred: string[] = [];
  const fallback: string[] = [];

  for (const attribute of attributes) {
    if (!isObject(attribute)) continue;
    const value = typeof attribute.value === 'string' ? attribute.value.trim() : '';
    if (!value) continue;
    const term = typeof attribute.term === 'string' ? attribute.term : '';
    const label = identifierLabel(term);
    if (preferredKeywords.some((keyword) => label.includes(keyword))) {
      preferred.push(value);
    } else {
      fallback.push(value);
    }
  }

  return uniqueStrings([...preferred, ...fallback]).slice(0, 20);
}

export function getEntityDescriptions(entity: EntityLike): string[] {
  return mapEntityAttributesToDescriptions(entity.entityAttributes);
}

const DESCRIPTION_ATTRIBUTE_TERM_KEYS = [
  'OM:0603:Function',
  'OM:0605:Disease',
  'OM:0604:Subcellular Location',
  'function',
  'disease_involvement',
  'subcellular_location',
] as const;

export function getAllowedEntityDescriptions(entity: EntityLike): string[] {
  if (!Array.isArray(entity.entityAttributes)) return [];

  const allowedTerms = new Set(DESCRIPTION_ATTRIBUTE_TERM_KEYS.map((term) => term.toLowerCase()));
  const allowedLabels = new Set(
    DESCRIPTION_ATTRIBUTE_TERM_KEYS.map((term) => getIdentifierTypeLabel(term).toLowerCase()),
  );

  return entity.entityAttributes.flatMap((attribute) => {
    if (!isObject(attribute)) return [];
    const term = typeof attribute.term === 'string' ? attribute.term.trim().toLowerCase() : '';
    const value = typeof attribute.value === 'string' ? attribute.value.trim() : '';
    const label = getIdentifierTypeLabel(term).toLowerCase();
    if (!term || !value || (!allowedTerms.has(term) && !allowedLabels.has(label))) return [];
    return [value];
  });
}

function getPreferredName(names: string[]): string | undefined {
  const scored = uniqueStrings(names).map((name, index) => {
    const trimmed = name.trim();
    let score = 0;
    if (/^https?:\/\//i.test(trimmed)) score += 100;
    if (/^(MLS|SMR|cid_|ZINC|SID_|CID_|InChI=)/i.test(trimmed)) score += 100;
    if (/^[A-Z0-9._-]{1,3}$/.test(trimmed)) score += 10;
    if (/^[A-Z][0-9A-Z]{3,}$/.test(trimmed)) score += 8;
    if (trimmed.length > 80) score += 25;
    if (trimmed.length > 40) score += 8;
    if (trimmed.length < 4) score += 4;
    if (/^[A-Z][a-z]/.test(trimmed)) score -= 3;
    return { name: trimmed, score, index };
  });

  scored.sort((a, b) => a.score - b.score || a.name.length - b.name.length || a.index - b.index);
  return scored[0]?.name;
}

export function getIdentifiersByType(entity: EntityLike, types: string[]): string[] {
  const values: string[] = [];
  for (const identifier of getEntityIdentifiers(entity)) {
    const value = identifier.value.trim();
    if (identifierTypeMatches(identifier.key, types) && value) {
      values.push(value);
    }
  }
  return uniqueStrings(values);
}

function getIdentifierByType(entity: EntityLike, types: string[]): string | undefined {
  return getIdentifiersByType(entity, types)[0];
}

function getFallbackEntityName(names: string[], geneSymbols: string[]): string | undefined {
  return geneSymbols[0] || getPreferredName(names) || names[0];
}

export function getEntityDisplayName(entity: EntityLike): string {
  if (entity.displayName?.trim()) return entity.displayName;
  const { names, geneSymbols } = classifyEntityIdentifiers(entity);
  const publicId = getEntityPublicId(entity);
  const entityTypeLabel = getEntityTypeLabel(entity).toLowerCase();

  if (isFallbackEntity(entity)) {
    return (
      getFallbackEntityName(names, geneSymbols) ||
      getShortestAvailableIdentifier(entity)?.value ||
      'Entity'
    );
  }

  if (isChemicalEntity(entity)) {
    const preferredName = getPreferredName(names);
    const chebiId = getIdentifierByType(entity, ['chebi']);
    const pubchemId = getIdentifierByType(entity, ['pubchem compound', 'pubchem']);
    const hmdbId = getIdentifierByType(entity, ['hmdb']);
    const chemblId = getIdentifierByType(entity, ['chembl']);
    if (preferredName && !/^\d+$/.test(preferredName) && !/^InChI=/i.test(preferredName))
      return preferredName;
    const identifierFallbacks = [chebiId, chemblId, hmdbId, pubchemId];
    return (
      identifierFallbacks.find((identifier) => identifier && !/^\d+$/.test(identifier)) ||
      identifierFallbacks.find(Boolean) ||
      preferredName ||
      entity.canonicalIdentifier ||
      publicId
    );
  }

  if (entityTypeLabel === 'protein' || entityTypeLabel === 'gene') {
    const uniprotId = getIdentifierByType(entity, ['uniprot', 'uniprotkb']);
    return geneSymbols[0] || uniprotId || names[0] || entity.canonicalIdentifier || publicId;
  }

  return geneSymbols[0] || names[0] || entity.canonicalIdentifier || publicId;
}

export function getEntitySecondaryName(entity: EntityLike): string | undefined {
  const { names, geneSymbols } = classifyEntityIdentifiers(entity);
  const entityTypeLabel = getEntityTypeLabel(entity).toLowerCase();

  if (isFallbackEntity(entity)) {
    const primaryName = getFallbackEntityName(names, geneSymbols);
    return getShortestAvailableIdentifier(entity, primaryName ? [primaryName] : [])?.value;
  }

  if (isCvTermEntity(entity)) {
    const primary = names[0];
    const canonicalIdentifier = entity.canonicalIdentifier || undefined;
    return primary && canonicalIdentifier && primary !== canonicalIdentifier
      ? canonicalIdentifier
      : undefined;
  }

  if (entityTypeLabel === 'protein' || entityTypeLabel === 'gene') {
    const canonicalIdentifier = entity.canonicalIdentifier || undefined;
    const primary = geneSymbols[0];
    if (primary && canonicalIdentifier && primary !== canonicalIdentifier) {
      return canonicalIdentifier;
    }
    const uniprotId = getIdentifierByType(entity, ['uniprot', 'uniprotkb']);
    if (uniprotId && uniprotId !== primary) return uniprotId;
    return names[0] && names[0] !== primary ? names[0] : undefined;
  }

  if (isChemicalEntity(entity)) {
    return (
      getIdentifierByType(entity, ['chebi']) ||
      getIdentifierByType(entity, ['pubchem compound', 'pubchem']) ||
      getIdentifierByType(entity, ['hmdb']) ||
      getIdentifierByType(entity, ['chembl']) ||
      undefined
    );
  }

  if (geneSymbols[0] && names[0] && geneSymbols[0] !== names[0]) {
    return names[0];
  }

  return undefined;
}

export function getEntityBadgeIdentifier(entity: EntityLike): string {
  const secondaryName = getEntitySecondaryName(entity);
  if (secondaryName) return secondaryName;

  const canonicalIdentifier = entity.canonicalIdentifier?.trim() || '';
  if (isFallbackEntity(entity) || isHashLikeIdentifier(canonicalIdentifier)) {
    return getShortestAvailableIdentifier(entity)?.value || getEntityDisplayName(entity);
  }

  return canonicalIdentifier || getEntityDisplayName(entity);
}

export function getEntitySmiles(entity: EntityLike): string | null {
  return (
    getIdentifierByType(entity, [
      'smiles',
      'canonical smiles',
      'canonical_smiles',
      'biotin tag',
      'biotin',
    ]) || null
  );
}
