import type { SearchFilters } from '../types/search';

export type IdentifierLike = { key: string; value: string };

export type ScopeEntity = {
  entityPk: string;
  entityType?: string | null;
  label?: string | null;
  displayName?: string | null;
  canonicalIdentifier?: string | null;
};

/** One way to narrow an entity's relations: the whole entity, or one product or group member. */
export type RelationScope = {
  id: string;
  kind: 'all' | 'product' | 'member';
  /** Short text for the scope button. */
  label: string;
  /** Human name and identifier of the product or member. */
  title?: string;
  identifier?: string;
  filters: SearchFilters;
};

export type MolecularContextLike = {
  referenceEntityKey?: string | null;
  filters?: SearchFilters | null;
  catalogueProducts?: ScopeEntity[];
};

export const ALL_SCOPE = 'all';

function entityName(entity: ScopeEntity): string {
  return (
    entity.displayName?.trim() ||
    entity.label?.trim() ||
    entity.canonicalIdentifier?.trim() ||
    entity.entityPk
  );
}

/**
 * Scopes for the relations tab. A gene reference narrows by the product its evidence names; a
 * group narrows by member; anything else is a single scope over the entity itself.
 */
export function relationScopes(
  entity: ScopeEntity,
  options: {
    context?: MolecularContextLike | null;
    memberKeys?: readonly string[];
    members?: readonly ScopeEntity[];
  } = {},
): { scopes: RelationScope[]; initial: string } {
  const { context, memberKeys = [], members = [] } = options;
  const name = entityName(entity);
  const base = context?.filters;
  if (context?.referenceEntityKey && base && Object.keys(base).length) {
    const products = (context.catalogueProducts ?? []).map((product) => ({
      id: product.entityPk,
      kind: 'product' as const,
      label: product.canonicalIdentifier?.trim() || entityName(product),
      title: entityName(product),
      identifier: product.canonicalIdentifier?.trim() || undefined,
      filters: {
        ...base,
        [product.entityType === 'protein' ? 'protein_entity_keys' : 'transcript_entity_keys']: [
          product.entityPk,
        ],
      },
    }));
    return {
      scopes: [{ id: ALL_SCOPE, kind: 'all', label: `All of ${name}`, filters: base }, ...products],
      initial: products.some((product) => product.id === entity.entityPk)
        ? entity.entityPk
        : ALL_SCOPE,
    };
  }
  if (memberKeys.length > 1) {
    const byKey = new Map(members.map((member) => [member.entityPk, member]));
    return {
      scopes: [
        {
          id: ALL_SCOPE,
          kind: 'all',
          label: `All of ${name}`,
          filters: { scope_entity_ids: [...memberKeys] },
        },
        ...memberKeys.flatMap((key) => {
          const member = byKey.get(key);
          return member
            ? [
                {
                  id: key,
                  kind: 'member' as const,
                  label: entityName(member),
                  title: entityName(member),
                  identifier: member.canonicalIdentifier?.trim() || undefined,
                  filters: { scope_entity_ids: [key] },
                },
              ]
            : [];
        }),
      ],
      initial: ALL_SCOPE,
    };
  }
  return {
    scopes: [
      { id: ALL_SCOPE, kind: 'all', label: name, filters: { scope_entity_ids: [entity.entityPk] } },
    ],
    initial: ALL_SCOPE,
  };
}

// Identifier types worth showing next to the name, most specific first.
const KEY_IDENTIFIER_TYPES = [
  'entrez',
  'ncbigene',
  'hgnc',
  'ensg',
  'ensembl',
  'uniprot',
  'mirbase',
  'chebi',
  'hmdb',
  'kegg',
  'pubchem',
  'lipidmaps',
  'swisslipids',
  'rhea',
  'reactome',
  'inchikey',
];

function keyTypeRank(key: string): number {
  const normalized = key.toLowerCase().replace(/[^a-z]/g, '');
  return KEY_IDENTIFIER_TYPES.findIndex((type) => normalized.startsWith(type));
}

/**
 * The primary identifier plus well-known identifier types that have exactly one value. A type
 * with several values (a gene's many UniProt accessions) would show an arbitrary one, so it is
 * left to the identifiers tab.
 */
export function keyIdentifiers(
  primary: IdentifierLike | null,
  identifiers: readonly IdentifierLike[],
  limit = 4,
): IdentifierLike[] {
  const valuesByRank = new Map<number, Set<string>>();
  for (const identifier of identifiers) {
    const rank = keyTypeRank(identifier.key);
    const value = identifier.value.trim();
    if (rank < 0 || !value) continue;
    valuesByRank.set(rank, (valuesByRank.get(rank) ?? new Set()).add(value));
  }
  const picked: IdentifierLike[] = [];
  const ranks = new Set<number>();
  if (primary?.value && primary.value !== 'Unavailable') {
    picked.push(primary);
    ranks.add(keyTypeRank(primary.key));
  }
  for (const [rank, values] of [...valuesByRank].sort(([a], [b]) => a - b)) {
    if (picked.length >= limit) break;
    const [value] = values;
    if (values.size !== 1 || ranks.has(rank) || picked.some((item) => item.value === value))
      continue;
    picked.push(identifiers.find((item) => item.value.trim() === value)!);
    ranks.add(rank);
  }
  return picked;
}

// First match wins: specific types precede their family; '' means no link.
const IDENTIFIERS_ORG_NAMESPACES: Array<[string, string]> = [
  ['uniprot entry', ''],
  ['uniprot keyword', ''],
  ['uniprot', 'uniprot'],
  ['chebi', 'chebi'],
  ['chembl target component', ''],
  ['chembl target', 'chembl.target'],
  ['chembl variant', ''],
  ['chembl', 'chembl.compound'],
  ['pubchem substance', 'pubchem.substance'],
  ['hmdb', 'hmdb'],
  ['pubchem', 'pubchem.compound'],
  ['ensembl', 'ensembl'],
  ['hgnc', 'hgnc'],
  ['entrez', 'ncbigene'],
  ['ncbi gene', 'ncbigene'],
  ['taxonomy', 'taxonomy'],
  ['tax id', 'taxonomy'],
  ['reactome', 'reactome'],
  ['interpro', 'interpro'],
  ['rhea', 'rhea'],
];

/** identifiers.org link for a typed identifier, when its namespace is known. */
export function identifiersOrgHref(typeText: string, value: string): string | null {
  const trimmed = value.trim();
  if (!trimmed) return null;
  const text = typeText.toLowerCase();
  const namespace = IDENTIFIERS_ORG_NAMESPACES.find(([needle]) => text.includes(needle))?.[1];
  const compact = /^[A-Za-z][A-Za-z0-9_.-]*:\S+$/.test(trimmed)
    ? trimmed
    : namespace
      ? `${namespace}:${trimmed}`
      : null;
  if (!compact) return null;
  const separator = compact.indexOf(':');
  return `https://identifiers.org/${encodeURIComponent(compact.slice(0, separator))}:${encodeURIComponent(compact.slice(separator + 1))}`;
}
