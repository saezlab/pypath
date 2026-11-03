/** Display quantities without discarding their source field or comparison bound. */
export function measurementPresentation(record: Record<string, unknown>):
  | {
      value: string;
      unit?: string;
      sourceField?: string;
    }
  | undefined {
  let candidate = record.quantity;
  if (
    !candidate &&
    typeof record.value === 'string' &&
    record.value.startsWith('{"has_numeric_value":')
  ) {
    try {
      candidate = JSON.parse(record.value);
    } catch {
      return undefined;
    }
  }
  if (!candidate || typeof candidate !== 'object') return undefined;
  const quantity = candidate as Record<string, unknown>;
  const value = quantity.has_numeric_value;
  if (typeof value !== 'number' || !Number.isFinite(value)) return undefined;
  const operators: Record<string, string> = { less_than: '<', equal_to: '=', greater_than: '>' };
  const comparator =
    typeof quantity.comparator === 'string'
      ? quantity.comparator
      : operators[String(quantity.has_binary_relation)] || '';
  const unit = [quantity.has_unit_prefix, quantity.has_unit]
    .filter((part): part is string => typeof part === 'string' && part.length > 0)
    .join(' ');
  return {
    value: `${comparator ? comparator + ' ' : ''}${Number(value.toPrecision(7))}`,
    unit: unit || undefined,
    sourceField: typeof quantity.source_field === 'string' ? quantity.source_field : undefined,
  };
}

export function publicationPmid(term: string, value: string | undefined): string | undefined {
  if (!value) return undefined;
  if (term.toLowerCase() === 'publications') return /^PMID:(\d+)$/i.exec(value.trim())?.[1];
  if (['pmid', 'pubmed'].includes(term.toLowerCase()))
    return /^(?:PMID:)?(\d+)$/i.exec(value.trim())?.[1];
  return undefined;
}
