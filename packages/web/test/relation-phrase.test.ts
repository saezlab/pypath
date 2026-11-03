import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import { relationPhrase } from '../src/lib/domain/relation-phrase.ts';

test('effects read as a verb phrase with their aspect and mechanism', () => {
  assert.deepEqual(
    relationPhrase(
      {
        predicate: 'affects',
        sign: 1,
        qualifiers: {
          object_direction_qualifier: 'increased',
          object_aspect_qualifier: 'activity',
          causal_mechanism_qualifier: 'partial_agonism',
        },
      },
      'Affects (increased)',
    ),
    {
      text: 'Increases activity of',
      short: 'activity',
      detail: 'via partial agonism',
      direction: 'up',
    },
  );
  assert.deepEqual(relationPhrase({ predicate: 'affects', sign: -1 }, 'Affects (decreased)'), {
    text: 'Decreases',
    short: 'Decreases',
    detail: undefined,
    direction: 'down',
  });
  assert.equal(
    relationPhrase(
      { predicate: 'affects', qualifiers: { object_aspect_qualifier: 'expression' } },
      'Affects',
    ).text,
    'Affects expression of',
  );
  // Without a direction the short phrase keeps the verb: there is no arrow to say it.
  assert.equal(
    relationPhrase(
      { predicate: 'affects', qualifiers: { object_aspect_qualifier: 'expression' } },
      'Affects',
    ).short,
    'Affects expression',
  );
});

test('qualified labels keep their wording and do not repeat their mechanism', () => {
  assert.deepEqual(
    relationPhrase(
      {
        predicate: 'interacts_with',
        displayLabel: 'Binds',
        qualifiers: { causal_mechanism_qualifier: 'binding' },
      },
      'Interacts with',
    ),
    { text: 'Binds', short: 'Binds', detail: undefined, direction: undefined },
  );
  assert.deepEqual(relationPhrase({ predicate: 'has_participant' }, 'Has participant'), {
    text: 'Has participant',
    short: 'Has participant',
    detail: undefined,
    direction: undefined,
  });
});
