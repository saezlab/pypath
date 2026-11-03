import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import { selectedFirst, splitEmptyOptions } from '../src/lib/utils/facets.ts';

test('zero-count facet options fold away unless selected or not yet counted', () => {
  const counts: Record<string, number | undefined> = {
    kegg: 20,
    brenda: 0,
    chebi: 0,
    new: undefined,
  };
  const { shown, empty } = splitEmptyOptions(
    ['kegg', 'brenda', 'chebi', 'new'],
    (option) => counts[option],
    (option) => option === 'chebi',
  );
  // The selected option leads the list.
  assert.deepEqual(shown, ['chebi', 'kegg', 'new']);
  assert.deepEqual(empty, ['brenda']);
});

test('selected options move to the top and both groups keep their order', () => {
  const chosen = new Set(['c', 'a']);
  const options = ['a', 'b', 'c', 'd'];
  assert.deepEqual(
    selectedFirst(options, (option) => chosen.has(option)),
    ['a', 'c', 'b', 'd'],
  );
  assert.deepEqual(options, ['a', 'b', 'c', 'd']);
  assert.deepEqual(
    selectedFirst(options, () => false),
    options,
  );
});
