export type ProductSummary = {
  entityPk: string;
  displayName?: string | null;
  label?: string | null;
  canonicalIdentifier?: string | null;
  canonicalIdentifierType?: string | null;
};

export function productLabel(key: string, products: readonly ProductSummary[]): string | undefined {
  const product = products.find((entry) => entry.entityPk === key);
  if (!product) return undefined;
  const label = product.displayName || product.label || '';
  const accession = product.canonicalIdentifier || '';
  if (label && accession && label !== accession) return `${label} · ${accession}`;
  return label || accession || undefined;
}

export function observationSummary(
  relation: { evidenceCount: number; sources: string[] },
  evidence?: readonly { source?: string | null }[],
): { evidenceCount: number; sources: string[] } {
  if (evidence === undefined) return relation;
  return {
    evidenceCount: evidence.length,
    sources: [...new Set(evidence.map((row) => row.source || 'Unknown'))].sort(),
  };
}

/** A sequence position or range; an unknown endpoint shows as "?". */
export function featureRange(position: unknown, end?: unknown): string {
  const start = typeof position === 'number' ? String(position) : '?';
  const stop = typeof end === 'number' ? String(end) : end == null ? start : '?';
  if (start === '?' && stop === '?') return '';
  return start === stop ? start : `${start}–${stop}`;
}
