import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import { measurementPresentation, publicationPmid } from '../src/lib/utils/measurements.ts';

test('quantity presentation preserves bounds, units, source fields and zero', () => {
  const quantity = { has_numeric_value: 0, has_unit: 'nM', comparator: '<=', source_field: 'IC50' };
  assert.deepEqual(measurementPresentation({ quantity }), {
    value: '<= 0',
    unit: 'nM',
    sourceField: 'IC50',
  });
  assert.deepEqual(
    measurementPresentation({ value: JSON.stringify(quantity) }),
    measurementPresentation({ quantity }),
  );
  assert.equal(measurementPresentation({ value: 'ordinary prose' }), undefined);
  assert.equal(measurementPresentation({ quantity: { has_numeric_value: NaN } }), undefined);
  assert.equal(
    measurementPresentation({
      quantity: { has_numeric_value: 4, has_binary_relation: 'greater_than' },
    })?.value,
    '> 4',
  );
});

test('native publication PMIDs are linked without interpreting DOI digits as PMIDs', () => {
  assert.equal(publicationPmid('publications', 'PMID:8784449'), '8784449');
  assert.equal(publicationPmid('publications', 'doi:10.1021/jm9602571'), undefined);
  assert.equal(publicationPmid('pubmed', '8784449'), '8784449');
});
