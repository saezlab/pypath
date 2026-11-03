import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import { publicId } from '../src/lib/domain/entity.ts';
import {
  plainText,
  annotationLabel,
  safeSourceUrl,
} from '../src/lib/utils/annotation-presentation.ts';
import { entityKind } from '../src/lib/domain/presentation.ts';

test('detail and selection identity uses only the exact key, never a numeric accession', () => {
  assert.equal(publicId({ entityPk: 'complex-key', canonicalIdentifier: '1' }), 'complex-key');
  assert.equal(publicId({ canonicalIdentifier: '1' }), '');
});
test('source text keeps reaction arrows and evidence qualification', () => {
  assert.equal(plainText('<b>ATP</b> <=> ADP &lt;=&gt; AMP'), 'ATP <=> ADP <=> AMP');
  assert.equal(plainText('F->L /evidence="ECO:123"'), 'F->L /evidence="ECO:123"');
  assert.equal(safeSourceUrl('javascript:alert(1)'), undefined);
  assert.equal(annotationLabel('member_sd'), 'Standard deviation');
});
test('record kind follows vocabulary semantics, not arbitrary titles or resource names', () => {
  assert.equal(
    entityKind({ entityType: 'molecular_activity', canonicalIdentifierType: 'rhea' }),
    'Reaction',
  );
  assert.equal(
    entityKind({ entityType: 'ontology_class', canonicalIdentifierType: 'go' }),
    'Ontology term',
  );
  assert.equal(
    entityKind({ entityType: 'information_content_entity', canonicalIdentifierType: 'pfocr' }),
    'Figure',
  );
  assert.equal(
    entityKind({ entityType: 'gene_family', canonicalIdentifierType: 'kegg_orthology' }),
    'Orthology group',
  );
  assert.equal(
    entityKind({ entityType: 'protein', canonicalIdentifierType: 'uniprot' }),
    undefined,
  );
});

test('wire entity adapters normalize optional API values and retain exact identity', async () => {
  const { entityFromWire, relationSearchFromWire } = await import('../src/lib/api/adapters.ts');
  const entity = entityFromWire({
    entityPk: 'exact-key',
    canonicalIdentifier: null,
    sources: undefined,
    identifiers: [{ id: 'P1', type: 'uniprot' }],
  });
  assert.equal(entity.entityPk, 'exact-key');
  assert.equal(entity.canonicalIdentifier, '');
  assert.deepEqual(entity.sources, []);
  assert.equal(entity.identifiers[0].entityPk, 'exact-key');
  assert.equal(entity.identifiers[0].identifier, 'P1');
  assert.throws(() => entityFromWire({ entityPk: '' }), /exact key/);
  const defaults = {
    annotationsTotal: 0,
    childCount: 0,
    hasHierarchy: false,
    identifiersTotal: 0,
    molecularEvidenceTotal: 0,
    parentCount: 0,
    relationCount: 0,
    sourceCount: 0,
  };
  const response = relationSearchFromWire({
    total: 1,
    rows: [
      {
        relation: { relationPk: 'r', isDirected: false, sign: 0 },
        subjectEntity: { entityPk: 's', ...defaults },
        objectEntity: { entityPk: 'o', ...defaults },
      },
    ],
  });
  assert.equal(response.rows[0].subjectEntity.entityPk, 's');
  assert.deepEqual(response.rows[0].relation.sources, []);
});
