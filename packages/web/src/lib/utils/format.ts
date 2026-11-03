/**
 * Format a number for display with appropriate separators and units
 * @param value - The number to format
 * @param options - Optional formatting options
 * @returns Formatted string
 */
export function formatNumber(
  value: number | undefined | null,
  options: {
    /** Whether to use compact notation (e.g., 1.2M instead of 1,200,000) */
    compact?: boolean;
    /** Minimum number of fraction digits (default: 0) */
    minimumFractionDigits?: number;
    /** Maximum number of fraction digits (default: 2 for decimal numbers, 0 for integers) */
    maximumFractionDigits?: number;
    /** Locale for formatting (default: 'en-US') */
    locale?: string;
    /** Custom suffix to append (e.g., ' interactions') */
    suffix?: string;
    /** Custom prefix to prepend */
    prefix?: string;
    /** Return "0" for zero values instead of the formatted zero */
    zeroValue?: string;
  } = {},
): string {
  if (value === undefined || value === null || isNaN(value)) {
    return options.zeroValue ?? '0';
  }

  if (value === 0 && options.zeroValue !== undefined) {
    return options.zeroValue;
  }

  const {
    compact = true,
    minimumFractionDigits = 0,
    maximumFractionDigits = Number.isInteger(value) ? 0 : 2,
    locale = 'en-US',
    suffix = '',
    prefix = '',
  } = options;

  const formatter = new Intl.NumberFormat(locale, {
    notation: compact ? 'compact' : 'standard',
    minimumFractionDigits,
    maximumFractionDigits,
    compactDisplay: 'short',
  });

  const formatted = formatter.format(value);
  return `${prefix}${formatted}${suffix}`;
}

/**
 * Format a count with appropriate singular/plural suffix
 * @param count - The count to format
 * @param singular - The singular form of the item (e.g., 'interaction')
 * @param plural - The plural form of the item (default: singular + 's')
 * @param options - Additional formatting options
 * @returns Formatted string with count and appropriate singular/plural form
 */
export function formatCount(
  count: number | undefined | null,
  singular: string,
  plural?: string,
  options: Omit<Parameters<typeof formatNumber>[1], 'suffix'> = {},
): string {
  const actualCount = count ?? 0;
  const pluralForm = plural ?? `${singular}s`;
  const itemForm = actualCount === 1 ? singular : pluralForm;

  return formatNumber(actualCount, {
    ...options,
    suffix: ` ${itemForm}`,
  });
}

/**
 * Format facet counts for display (e.g., "Protein (1,234)")
 * @param label - The facet label
 * @param count - The count to display
 * @param options - Formatting options
 * @returns Formatted facet string
 */
export function formatFacetCount(
  label: string,
  count: number,
  options: {
    /** Whether to use compact notation for large numbers */
    compact?: boolean;
    /** Whether to show count in parentheses (default: true) */
    showParentheses?: boolean;
  } = {},
): string {
  const { compact = count > 9999, showParentheses = true } = options;

  const formattedCount = formatNumber(count, { compact });

  if (showParentheses) {
    return `${label} (${formattedCount})`;
  }

  return `${label}: ${formattedCount}`;
}

/**
 * Format large numbers with abbreviated units (K, M, B)
 * @param value - The number to format
 * @returns Formatted string with unit
 */
export function formatCompactNumber(value: number | undefined | null): string {
  return formatNumber(value, { compact: true });
}

/**
 * Format percentage values
 * @param value - The decimal value (0-1) or percentage value (0-100)
 * @param options - Formatting options
 * @returns Formatted percentage string
 */
export function formatPercentage(
  value: number | undefined | null,
  options: {
    /** Whether the input is already a percentage (0-100) or decimal (0-1) */
    isPercentage?: boolean;
    /** Number of decimal places to show */
    decimals?: number;
  } = {},
): string {
  if (value === undefined || value === null || isNaN(value)) {
    return '0%';
  }

  const { isPercentage = false, decimals = 0 } = options;
  const percentValue = isPercentage ? value : value * 100;

  return formatNumber(percentValue, {
    minimumFractionDigits: decimals,
    maximumFractionDigits: decimals,
    suffix: '%',
  });
}
