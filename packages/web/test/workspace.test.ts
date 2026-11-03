import { strict as assert } from 'node:assert';
import { test } from 'node:test';
import {
  defaultLayout,
  normalizeLayout,
  dock,
  resizeSplit,
  geometry,
  leaves,
  prune,
  type Tile,
} from '../src/lib/components/workspace/layout.ts';
const ids = ['results', 'types', 'sources', 'taxonomy'];
function check(tree: Tile) {
  assert.deepEqual([...leaves(tree)].sort(), [...ids].sort());
  const branchIds: string[] = [];
  const collect = (n: Tile) => {
    if (n.type === 'split') {
      branchIds.push(n.id);
      collect(n.first);
      collect(n.second);
    }
  };
  collect(tree);
  assert.equal(new Set(branchIds).size, branchIds.length);
  const { panels } = geometry(tree, { x: 0, y: 0, width: 1200, height: 700 });
  const boxes = Object.values(panels);
  for (const b of boxes) {
    assert.ok(b.width > 0);
    assert.ok(b.height > 0);
  }
  for (let i = 0; i < boxes.length; i++)
    for (let j = i + 1; j < boxes.length; j++) {
      const a = boxes[i],
        b = boxes[j];
      assert.ok(
        !(
          a.x < b.x + b.width &&
          a.x + a.width > b.x &&
          a.y < b.y + b.height &&
          a.y + a.height > b.y
        ),
      );
    }
}
test('tiling splits and swaps preserve every panel with no overlap', () => {
  let tree = defaultLayout(ids).tree;
  for (let i = 0; i < 80; i++) {
    tree = dock(
      tree,
      ids[i % 4],
      ids[(i + 1) % 4],
      (['left', 'right', 'top', 'bottom', 'center'] as const)[i % 5],
    );
    check(tree);
  }
});
test('resizing a shared split redistributes space while preserving total bounds', () => {
  const tree = defaultLayout(ids).tree;
  const next = resizeSplit(tree, 'root', 0.5);
  check(next);
  const boxes = geometry(next, { x: 0, y: 0, width: 1210, height: 700 }).panels;
  assert.equal(boxes.results.width, 600);
  assert.equal(boxes.types.x, 610);
  assert.equal((resizeSplit(tree, 'root', 2) as { ratio: number }).ratio, 0.9);
});
test('saved trees validate duplicate leaves and removed panels, add new panels, and prune collapsed tiles', () => {
  const initial = defaultLayout(ids);
  assert.deepEqual(normalizeLayout(JSON.parse(JSON.stringify(initial)), ids), initial);
  assert.deepEqual(normalizeLayout({ version: 2 }, ids), initial);
  const corrupt = {
    version: 3,
    tree: {
      type: 'split',
      id: 'r',
      axis: 'x',
      ratio: Infinity,
      first: { type: 'leaf', id: 'results' },
      second: { type: 'leaf', id: 'results' },
    },
    collapsed: ['removed', 'types'],
  };
  const restored = normalizeLayout(corrupt, ids);
  check(restored.tree);
  assert.deepEqual(restored.collapsed, ['types']);
  const visible = prune(initial.tree, new Set(['results', 'taxonomy']));
  assert.deepEqual(leaves(visible!).sort(), ['results', 'taxonomy']);
});

test('relation workspace defaults to filters and effect over three filter columns', () => {
  const ids = ['results', 'filters', 'effect', 'interaction_types', 'sources', 'ncbi_tax_id'];
  const { panels } = geometry(defaultLayout(ids).tree, { x: 0, y: 0, width: 2000, height: 1400 });
  assert.ok(panels.results.width > 1000 && panels.results.height === 1400);
  assert.equal(panels.filters.y, panels.effect.y);
  assert.ok(panels.filters.width > panels.effect.width);
  const columns = ['sources', 'interaction_types', 'ncbi_tax_id'].map((id) => panels[id]);
  assert.ok(columns.every((rect) => rect.y > panels.filters.y + panels.filters.height));
  assert.ok(columns.every((rect) => Math.abs(rect.width - columns[0].width) < 5));
});
