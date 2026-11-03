/** The options with the selected ones first, each group in its given order. */
export function selectedFirst<T>(options: readonly T[], isSelected: (option: T) => boolean): T[] {
  return [...options.filter(isSelected), ...options.filter((option) => !isSelected(option))];
}

/**
 * Splits facet options into those with results and those without. Options whose count is
 * unknown stay visible, and selected options stay visible even at zero so they can be cleared.
 * The selected options come first.
 */
export function splitEmptyOptions<T>(
  options: readonly T[],
  count: (option: T) => number | null | undefined,
  isSelected: (option: T) => boolean,
): { shown: T[]; empty: T[] } {
  const shown: T[] = [];
  const empty: T[] = [];
  for (const option of options) {
    (count(option) === 0 && !isSelected(option) ? empty : shown).push(option);
  }
  return { shown: selectedFirst(shown, isSelected), empty };
}
