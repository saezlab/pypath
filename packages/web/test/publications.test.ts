import assert from 'node:assert/strict';
import test from 'node:test';
import { groupPublications, isPublicationTerm } from '../src/lib/utils/publications';

test('native and legacy publications are grouped with their sources', () => {
  const grouped = groupPublications([
    { term: 'publications', value: 'PMID:12345', source: 'a' },
    { term: 'pmid', value: '12345; 67890', source: 'b' },
    { term: 'publications', value: 'doi:10.1234/56789', source: 'a' },
    { term: 'publications', value: 'https://doi.org/10.1234/56789', source: 'b' },
    { term: 'publications', value: 'Unknown reference', source: 'a' },
    { term: 'description', value: '12345', source: 'c' },
  ]);
  assert.equal(grouped.length, 4);
  assert.deepEqual(grouped.find((p) => p.id === 'PMID:12345')?.sources, ['a', 'b']);
  assert.equal(
    grouped.find((p) => p.id === 'DOI:10.1234/56789')?.url,
    'https://doi.org/10.1234/56789',
  );
  assert.equal(grouped.find((p) => p.id === 'Unknown reference')?.url, undefined);
  assert.ok(!grouped.some((p) => p.id === 'PMID:56789'));
  assert.ok(isPublicationTerm('biolink:publications'));
  assert.ok(!isPublicationTerm('description'));
});
