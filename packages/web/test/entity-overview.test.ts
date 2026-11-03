import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import {
  identifiersOrgHref,
  keyIdentifiers,
  relationScopes,
} from '../src/lib/utils/entity-overview.ts';

const P04637 = 'b'.repeat(64);
const K7PPA8 = 'c'.repeat(64);

test('gene references scope relations by the product their evidence names', () => {
  const context = {
    referenceEntityKey: 'entrez:7157',
    filters: { reference_entity_keys: ['entrez:7157'] },
    catalogueProducts: [
      { entityPk: P04637, entityType: 'protein', label: 'TP53', canonicalIdentifier: 'P04637' },
      { entityPk: K7PPA8, entityType: 'protein', label: 'TP53', canonicalIdentifier: 'K7PPA8' },
    ],
  };
  const group = relationScopes({ entityPk: 'gene:entrez:7157', label: 'TP53' }, { context });
  assert.equal(group.initial, 'all');
  assert.deepEqual(
    group.scopes.map((scope) => [scope.kind, scope.label]),
    [
      ['all', 'All of TP53'],
      ['product', 'P04637'],
      ['product', 'K7PPA8'],
    ],
  );
  assert.deepEqual(group.scopes[1].filters, {
    reference_entity_keys: ['entrez:7157'],
    protein_entity_keys: [P04637],
  });

  // Opening one product starts on that product instead of the whole gene.
  const product = relationScopes({ entityPk: P04637, label: 'TP53' }, { context });
  assert.equal(product.initial, P04637);
});

test('groups without a gene reference scope relations by member', () => {
  const { scopes, initial } = relationScopes(
    { entityPk: 'connectivity:WQZGKKKJIJFFOK', label: 'Glucose' },
    {
      memberKeys: ['m1', 'm2'],
      members: [{ entityPk: 'm2', label: 'Beta-D-Glucose', canonicalIdentifier: 'WQZG-VFU' }],
    },
  );
  assert.equal(initial, 'all');
  assert.deepEqual(scopes[0].filters, { scope_entity_ids: ['m1', 'm2'] });
  // Members whose details did not load are left out rather than shown as hashes.
  assert.deepEqual(
    scopes.slice(1).map((scope) => [scope.kind, scope.label, scope.filters]),
    [['member', 'Beta-D-Glucose', { scope_entity_ids: ['m2'] }]],
  );
});

test('a single entity has one scope over itself', () => {
  const { scopes } = relationScopes({ entityPk: 'x', label: 'ATP' });
  assert.deepEqual(scopes, [
    { id: 'all', kind: 'all', label: 'ATP', filters: { scope_entity_ids: ['x'] } },
  ]);
});

test('key identifiers keep the primary first and only unambiguous known types', () => {
  const picked = keyIdentifiers({ key: 'entrez', value: '7157' }, [
    { key: 'genesymbol', value: 'TP53' },
    { key: 'uniprot', value: 'P04637' },
    { key: 'uniprot', value: 'K7PPA8' },
    { key: 'ensg', value: 'ENSG00000141510' },
    { key: 'hgnc', value: 'HGNC:11998' },
    { key: 'entrez', value: '7157' },
  ]);
  assert.deepEqual(
    picked.map((identifier) => identifier.value),
    ['7157', 'HGNC:11998', 'ENSG00000141510'],
  );
  assert.deepEqual(keyIdentifiers({ key: 'Identifier', value: 'Unavailable' }, []), []);
});

test('identifiers.org links use known namespaces or existing CURIEs', () => {
  assert.equal(
    identifiersOrgHref('entrez Entrez Gene', '7157'),
    'https://identifiers.org/ncbigene:7157',
  );
  assert.equal(identifiersOrgHref('hgnc', 'HGNC:11998'), 'https://identifiers.org/HGNC:11998');
  assert.equal(identifiersOrgHref('genesymbol', 'TP53'), null);
  assert.equal(
    identifiersOrgHref('chembl target chembl target', 'CHEMBL203'),
    'https://identifiers.org/chembl.target:CHEMBL203',
  );
  assert.equal(identifiersOrgHref('chembl component chembl target component', '147'), null);
});
