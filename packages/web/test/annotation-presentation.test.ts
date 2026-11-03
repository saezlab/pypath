import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import {
  annotationDefaultUnit,
  annotationLabel,
  annotationLabelFor,
  collapseListedAnnotations,
  isNarrativeAnnotation,
  ontologyTermHref,
} from '../src/lib/utils/annotation-presentation.ts';

test('has_topic is labelled by the vocabulary of its value', () => {
  assert.equal(annotationLabelFor('has_topic', 'EC:2.7.7.48'), 'EC number');
  assert.equal(annotationLabelFor('has_topic', 'WP5566_r141133'), 'WikiPathways revision');
  assert.equal(annotationLabelFor('has_topic', 'something else'), 'Classification / reference');
  assert.equal(annotationLabelFor('stoichiometry', '2'), 'Stoichiometry');
});

test('masses default to daltons and transmembrane notes are not descriptions', () => {
  assert.equal(annotationDefaultUnit('chemrof:mass'), 'Da');
  assert.equal(annotationDefaultUnit('chemrof:charge'), undefined);
  assert.equal(isNarrativeAnnotation('up:Function_Annotation'), true);
  assert.equal(isNarrativeAnnotation('up:Transmembrane_Annotation'), false);
});

test('gene mapping candidates collapse into one row', () => {
  const rows = [
    { term: 'omnipath:gene_mapping_status', value: 'ambiguous' },
    ...Array.from({ length: 7 }, (_, i) => ({
      term: 'omnipath:gene_mapping_candidate',
      value: `entrez:${i}`,
    })),
  ];
  const collapsed = collapseListedAnnotations(rows);
  assert.equal(collapsed.length, 2);
  assert.equal(collapsed[1].value, 'entrez:0, entrez:1, entrez:2, entrez:3, entrez:4 and 2 more');
});

test('relation ontology CURIEs read as names', async () => {
  const { annotationLabelFor, annotationValueText } =
    await import('../src/lib/utils/annotation-presentation.ts');
  assert.equal(annotationLabelFor('RO:0018038'), 'Has functional parent');
  assert.equal(
    annotationValueText('original_predicate', 'RO:0002007'),
    'Bounding layer of (RO:0002007)',
  );
  assert.equal(annotationValueText('description', 'RO:0002007'), 'RO:0002007');
});

test('ontology term values link to identifiers.org', () => {
  assert.equal(ontologyTermHref('CHEBI:22587'), 'https://identifiers.org/CHEBI:22587');
  assert.equal(ontologyTermHref('http://example.org/x'), null);
  assert.equal(ontologyTermHref('Food part: leaf'), null);
  assert.equal(annotationLabel('RO:0000087'), 'Has role');
});
