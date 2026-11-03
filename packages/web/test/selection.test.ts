import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  entityPrimaryKeyFromSelection,
  selectionFromUrl,
} from '../src/lib/features/selection/model.ts';
import { buildSelectionUrl } from '../src/lib/features/selection/url-codecs.ts';
import { createUrlUpdater } from '../src/lib/features/selection/navigation.ts';
import {
  readStored,
  isStringList,
  writeStored,
} from '../src/lib/features/selection/persistence.ts';

const key = 'a'.repeat(64);
test('fresh shared links restore exact SHA-256 keys without a metadata cache', () => {
  assert.equal(entityPrimaryKeyFromSelection({ id: key, entityId: key, name: key }), key);
  assert.equal(entityPrimaryKeyFromSelection({ id: '1', entityId: '1', name: '1' }), null);
  const url = new URL(
    buildSelectionUrl({ entityIds: [key], annotationIds: [] }),
    'http://web.test',
  );
  assert.deepEqual(selectionFromUrl(url, ['unrelated'], ['stale-term']), {
    entityIds: [key],
    annotationIds: [],
  });
});
test('explicit empty selections own both lists; URLs without selection use local persistence', () => {
  for (const query of ['?entities=&annotations=', '?annotations=GO:1', '?selection=1']) {
    const value = selectionFromUrl(
      new URL(`http://web.test/selection${query}`),
      ['unrelated'],
      ['stale'],
    );
    assert.deepEqual(value.entityIds, []);
    assert.deepEqual(value.annotationIds, query.includes('GO:1') ? ['GO:1'] : []);
  }
  assert.deepEqual(selectionFromUrl(new URL('http://web.test/explore'), [key], []), {
    entityIds: [key],
    annotationIds: [],
  });
});
test('clear and synchronous edits use one complete navigation and preserve other URL state', async () => {
  let current = new URL('http://web.test/selection?entities=A&annotations=B&release=2');
  const calls: URL[] = [];
  const update = createUrlUpdater(
    () => current,
    async (url) => {
      calls.push(url);
    },
    (url) => {
      current = url;
    },
  );
  update({ entities: '', selection: '1' });
  update({ annotations: '' });
  await new Promise((resolve) => setImmediate(resolve));
  assert.equal(calls.length, 1);
  assert.equal(calls[0].searchParams.get('entities'), '');
  assert.equal(calls[0].searchParams.get('annotations'), '');
  assert.equal(calls[0].searchParams.get('release'), '2');
});
test('storage validation excludes malformed persisted IDs', () => {
  assert.deepEqual(readStored({ getItem: () => '[1,"x"]' }, 'key', [], isStringList), []);
  const stored: string[] = [];
  writeStored(
    {
      setItem: (_key, value) => {
        stored.push(value);
      },
    },
    'key',
    [key],
  );
  assert.deepEqual(JSON.parse(stored[0]), [key]);
});
