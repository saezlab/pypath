export type Leaf = { type: 'leaf'; id: string };
export type Split = {
  type: 'split';
  id: string;
  axis: 'x' | 'y';
  ratio: number;
  first: Tile;
  second: Tile;
};
export type Tile = Leaf | Split;
export type Layout = { version: 3; tree: Tile; collapsed: string[] };
export type Edge = 'left' | 'right' | 'top' | 'bottom' | 'center';
export type Rect = { x: number; y: number; width: number; height: number };
export type Divider = { id: string; axis: 'x' | 'y'; ratio: number; rect: Rect };
const leaf = (id: string): Leaf => ({ type: 'leaf', id });
function stack(ids: string[]): Tile {
  if (ids.length === 1) return leaf(ids[0]);
  return {
    type: 'split',
    id: `stack:${ids.join(':')}`,
    axis: 'y',
    ratio: 1 / ids.length,
    first: leaf(ids[0]),
    second: stack(ids.slice(1)),
  };
}
export function defaultLayout(ids: string[]): Layout {
  // Default entity workspace: results, full-height sources, then types over taxonomy.
  if (
    ids.length === 4 &&
    ['results', 'sources', 'entity_types', 'ncbi_tax_id'].every((id) => ids.includes(id))
  ) {
    return {
      version: 3,
      collapsed: [],
      tree: {
        type: 'split',
        id: 'root',
        axis: 'x',
        ratio: 0.54,
        first: leaf('results'),
        second: {
          type: 'split',
          id: 'filter-columns',
          axis: 'x',
          ratio: 0.475,
          first: leaf('sources'),
          second: {
            type: 'split',
            id: 'types-taxonomy',
            axis: 'y',
            ratio: 0.46,
            first: leaf('entity_types'),
            second: leaf('ncbi_tax_id'),
          },
        },
      },
    };
  }
  // Default relation workspace: results; relation filters and effect side by side over
  // sources, participant types and taxonomy in three equal columns.
  if (
    ids.length === 6 &&
    ['results', 'filters', 'effect', 'interaction_types', 'sources', 'ncbi_tax_id'].every((id) =>
      ids.includes(id),
    )
  ) {
    return {
      version: 3,
      collapsed: [],
      tree: {
        type: 'split',
        id: 'root',
        axis: 'x',
        ratio: 0.535,
        first: leaf('results'),
        second: {
          type: 'split',
          id: 'filter-rows',
          axis: 'y',
          ratio: 0.49,
          first: {
            type: 'split',
            id: 'relation-filters',
            axis: 'x',
            ratio: 0.54,
            first: leaf('filters'),
            second: leaf('effect'),
          },
          second: {
            type: 'split',
            id: 'filter-columns',
            axis: 'x',
            ratio: 1 / 3,
            first: leaf('sources'),
            second: {
              type: 'split',
              id: 'types-taxonomy',
              axis: 'x',
              ratio: 0.5,
              first: leaf('interaction_types'),
              second: leaf('ncbi_tax_id'),
            },
          },
        },
      },
    };
  }
  const others = ids.filter((id) => id !== 'results');
  return {
    version: 3,
    collapsed: [],
    tree: others.length
      ? {
          type: 'split',
          id: 'root',
          axis: 'x',
          ratio: others.length > 4 ? 0.6 : 0.68,
          first: leaf('results'),
          second:
            others.length > 4
              ? {
                  type: 'split',
                  id: 'filter-columns',
                  axis: 'x',
                  ratio: 0.5,
                  first: stack(others.slice(0, Math.ceil(others.length / 2))),
                  second: stack(others.slice(Math.ceil(others.length / 2))),
                }
              : stack(others),
        }
      : leaf(ids[0] || 'results'),
  };
}
export function leaves(tree: Tile): string[] {
  return tree.type === 'leaf' ? [tree.id] : [...leaves(tree.first), ...leaves(tree.second)];
}
export function prune(tree: Tile, ids: Set<string>): Tile | null {
  if (tree.type === 'leaf') return ids.has(tree.id) ? tree : null;
  const first = prune(tree.first, ids),
    second = prune(tree.second, ids);
  return first && second ? { ...tree, first, second } : first || second;
}
export function normalizeLayout(raw: unknown, ids: string[]): Layout {
  const fallback = defaultLayout(ids);
  if (!raw || typeof raw !== 'object' || (raw as Layout).version !== 3) return fallback;
  const seen = new Set<string>(),
    branches = new Set<string>();
  function read(v: unknown, depth = 0): Tile | null {
    if (!v || typeof v !== 'object' || depth > 32) return null;
    const n = v as Tile;
    if (n.type === 'leaf') {
      if (!ids.includes(n.id) || seen.has(n.id)) return null;
      seen.add(n.id);
      return leaf(n.id);
    }
    if (n.type !== 'split') return null;
    const first = read(n.first, depth + 1),
      second = read(n.second, depth + 1);
    if (!first || !second) return first || second;
    const id = typeof n.id === 'string' && !branches.has(n.id) ? n.id : `restored:${branches.size}`;
    branches.add(id);
    return {
      type: 'split',
      id,
      axis: n.axis === 'y' ? 'y' : 'x',
      ratio: Number.isFinite(n.ratio) ? Math.max(0.1, Math.min(0.9, n.ratio)) : 0.5,
      first,
      second,
    };
  }
  let tree = read((raw as Layout).tree);
  // The former individual relation filters cannot map to the combined filter panel.
  if (ids.length === 2 && ids.includes('filters') && !seen.has('filters')) return fallback;
  for (const id of ids.filter((id) => !seen.has(id)))
    tree = tree
      ? { type: 'split', id: `added:${id}`, axis: 'y', ratio: 0.75, first: tree, second: leaf(id) }
      : leaf(id);
  return {
    version: 3,
    tree: tree || fallback.tree,
    collapsed: Array.isArray((raw as Layout).collapsed)
      ? (raw as Layout).collapsed.filter((id) => ids.includes(id))
      : [],
  };
}
export function resizeSplit(tree: Tile, id: string, ratio: number): Tile {
  if (tree.type === 'leaf') return tree;
  if (tree.id === id) return { ...tree, ratio: Math.max(0.1, Math.min(0.9, ratio)) };
  return {
    ...tree,
    first: resizeSplit(tree.first, id, ratio),
    second: resizeSplit(tree.second, id, ratio),
  };
}
export function dock(tree: Tile, id: string, target: string, edge: Edge): Tile {
  if (id === target || !leaves(tree).includes(id) || !leaves(tree).includes(target)) return tree;
  if (edge === 'center') {
    const swap = (n: Tile): Tile =>
      n.type === 'leaf'
        ? leaf(n.id === id ? target : n.id === target ? id : n.id)
        : { ...n, first: swap(n.first), second: swap(n.second) };
    return swap(tree);
  }
  const rest = prune(tree, new Set(leaves(tree).filter((key) => key !== id)))!;
  const used = new Set<string>();
  const collect = (n: Tile) => {
    used.add(n.id);
    if (n.type === 'split') {
      collect(n.first);
      collect(n.second);
    }
  };
  collect(rest);
  const base = `dock:${id}:${target}:${edge}`;
  let branch = base;
  let suffix = 0;
  while (used.has(branch)) branch = `${base}:${++suffix}`;
  const insert = (n: Tile): Tile => {
    if (n.type === 'split') return { ...n, first: insert(n.first), second: insert(n.second) };
    if (n.id !== target) return n;
    const before = edge === 'left' || edge === 'top';
    return {
      type: 'split',
      id: branch,
      axis: edge === 'left' || edge === 'right' ? 'x' : 'y',
      ratio: 0.5,
      first: before ? leaf(id) : n,
      second: before ? n : leaf(id),
    };
  };
  return insert(rest);
}
export function geometry(
  tree: Tile | null,
  rect: Rect,
  gap = 10,
): { panels: Record<string, Rect>; dividers: Divider[] } {
  const panels: Record<string, Rect> = {},
    dividers: Divider[] = [];
  function visit(n: Tile, r: Rect) {
    if (n.type === 'leaf') {
      panels[n.id] = r;
      return;
    }
    const minimum = (node: Tile, axis: 'x' | 'y'): number =>
      node.type === 'leaf'
        ? axis === 'x'
          ? node.id === 'results'
            ? 280
            : 180
          : 90
        : node.axis === axis
          ? minimum(node.first, axis) + minimum(node.second, axis) + gap
          : Math.max(minimum(node.first, axis), minimum(node.second, axis));
    const span = n.axis === 'x' ? r.width : r.height,
      space = Math.min(gap, span / 4),
      total = span - space;
    const low = minimum(n.first, n.axis),
      high = minimum(n.second, n.axis);
    const first =
        total >= low + high
          ? Math.max(low, Math.min(total - high, total * n.ratio))
          : (total * low) / (low + high),
      second = total - first;
    dividers.push({ id: n.id, axis: n.axis, ratio: first / total, rect: r });
    if (n.axis === 'x') {
      visit(n.first, { ...r, width: first });
      visit(n.second, { ...r, x: r.x + first + space, width: second });
    } else {
      visit(n.first, { ...r, height: first });
      visit(n.second, { ...r, y: r.y + first + space, height: second });
    }
  }
  if (tree) visit(tree, rect);
  return { panels, dividers };
}
