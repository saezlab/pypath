import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import {
  filterResources,
  sortResources,
  resourceTags,
  formatBytes,
} from '../src/lib/resources/catalog.ts';
import { resourceFromWire } from '../src/lib/resources/types.ts';
import type { ResourceRecord, ResourceTag } from '../src/lib/resources/types.ts';

const signaling: ResourceTag = {
  id: 'signaling',
  label: 'Signaling',
  dimension: 'topic',
  description: '',
};
const metabolism: ResourceTag = {
  id: 'metabolism',
  label: 'Metabolism',
  dimension: 'topic',
  description: '',
};
const proteins: ResourceTag = {
  id: 'entity:protein',
  label: 'Protein',
  dimension: 'entity_type',
  description: '',
};
const make = (
  id: string,
  tags: ResourceTag[],
  license_use?: ResourceRecord['license_use'],
): ResourceRecord => ({
  resource_id: id,
  resource_name: id,
  tags,
  license_use,
  version: '1',
  entity_count: 10,
  interaction_count: 20,
  total_size_bytes: 100,
  files: [],
  description: 'Curated biological knowledge',
});
const resources = [
  make('Alpha', [signaling, proteins], { academic: 'allowed', commercial: 'requires_permission' }),
  make('Beta', [metabolism, proteins], { academic: 'allowed', commercial: 'allowed' }),
  make('Gamma', [metabolism]),
];

test('resource filters combine alternatives within groups and intersect different groups', () => {
  assert.equal(filterResources(resources, '', ['signaling', 'metabolism'], []).length, 3);
  assert.deepEqual(
    filterResources(resources, '', ['signaling', 'metabolism', 'entity:protein'], []).map(
      (r) => r.resource_id,
    ),
    ['Alpha', 'Beta'],
  );
  assert.equal(filterResources(resources, '', ['missing_tag'], []).length, 0);
});
test('resource search includes descriptions and tags, and combines with license filters', () => {
  assert.equal(filterResources(resources, ' BIOLOGICAL ', [], []).length, 3);
  assert.equal(filterResources(resources, 'signaling', [], []).length, 1);
  assert.deepEqual(
    filterResources(resources, '', [], ['academic']).map((r) => r.resource_id),
    ['Alpha', 'Beta'],
  );
  assert.deepEqual(
    filterResources(resources, '', [], ['academic', 'commercial']).map((r) => r.resource_id),
    ['Beta'],
  );
  assert.equal(filterResources(resources, '', [], ['invalid']).length, 0);
  assert.equal(filterResources(resources, 'Alpha', [], ['commercial']).length, 0);
});
test('vocabulary is deduplicated and sorting does not mutate API data', () => {
  assert.equal(resourceTags(resources).length, 3);
  const input = [resources[1], { ...resources[0], entity_count: 99 }];
  assert.equal(sortResources(input, 'entities')[0].resource_id, 'Alpha');
  assert.equal(sortResources(input, 'unknown')[0].resource_id, 'Alpha');
  const byEntities = sortResources(input, 'entities').map((r) => r.resource_id);
  assert.deepEqual(
    sortResources(input, 'entities', 'asc').map((r) => r.resource_id),
    [...byEntities].reverse(),
  );
  assert.equal(input[0].resource_id, 'Beta');
  assert.equal(formatBytes(1024 ** 3), '1.0 GB');
  assert.equal(formatBytes(0), '0 B');
});
test('the entity types a build holds become tags of their own group', () => {
  const resource = resourceFromWire({
    resource_id: 'signor',
    description: '',
    queryable: true,
    tags: [{ id: 'signaling', label: 'Signaling', dimension: 'topic', description: '' }],
    resolution_stats: {
      by_entity_type: { protein: { input: 5 }, chemical_entity: { input: 9 }, novel_type: {} },
    },
  });
  assert.deepEqual(
    resource.tags?.map((tag) => [tag.id, tag.label, tag.dimension]),
    [
      ['signaling', 'Signaling', 'topic'],
      ['entity:chemical_entity', 'Chemical entity', 'entity_type'],
      ['entity:protein', 'Protein', 'entity_type'],
      ['entity:novel_type', 'novel_type', 'entity_type'],
    ],
  );
  const legacy = resourceFromWire({ resource_id: 'legacy', description: '', queryable: true });
  assert.deepEqual(legacy.tags, []);
});
