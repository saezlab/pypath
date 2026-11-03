<script lang="ts">
  import { onMount, setContext, untrack, type Snippet } from 'svelte';
  import {
    defaultLayout,
    normalizeLayout,
    geometry,
    prune,
    dock,
    resizeSplit,
    leaves,
    type Layout,
    type Edge,
    type Divider,
  } from './layout';
  import { WORKSPACE, type WorkspaceContext } from './context';
  import { workspaceWidth } from '$lib/stores/workspace-width.svelte';
  let { name, panels, children }: { name: string; panels: string[]; children: Snippet } = $props();
  let layout = $state<Layout>(untrack(() => defaultLayout(panels)));
  let titles = $state<Record<string, string>>({});
  let active = $state<string[]>([]),
    width = $state(0),
    height = $state(700);
  let canvas: HTMLDivElement;
  let dragging = $state<string | null>(null),
    dragTitle = $state(''),
    pointer = $state({ x: 0, y: 0 });
  let destination = $state<{ id: string; edge: Edge } | null>(null);
  let resizing = $state<string | null>(null);
  const key = untrack(() => `omnipath-workspace-v3-${name}`);
  const visible = $derived(
    prune(layout.tree, new Set(active.filter((id) => !layout.collapsed.includes(id)))),
  );
  const positions = $derived(geometry(visible, { x: 0, y: 0, width, height }));
  const preview = $derived.by(() => {
    if (!destination) return null;
    const r = positions.panels[destination.id];
    if (!r) return null;
    const e = destination.edge;
    return {
      ...r,
      x: r.x + (e === 'right' ? r.width / 2 : 0),
      y: r.y + (e === 'bottom' ? r.height / 2 : 0),
      width: e === 'left' || e === 'right' ? r.width / 2 : r.width,
      height: e === 'top' || e === 'bottom' ? r.height / 2 : r.height,
    };
  });
  function save() {
    try {
      localStorage.setItem(key, JSON.stringify(layout));
    } catch {}
  }
  function cancel() {
    dragging = null;
    destination = null;
  }
  onMount(() => {
    try {
      layout = normalizeLayout(JSON.parse(localStorage.getItem(key) || 'null'), panels);
    } catch {}
    workspaceWidth.load();
  });
  function move(x: number, y: number) {
    if (!dragging || !canvas) return;
    const bounds = canvas.getBoundingClientRect();
    pointer = { x: x - bounds.left, y: y - bounds.top };
    const match = Object.entries(positions.panels).find(
      ([id, r]) =>
        id !== dragging &&
        pointer.x >= r.x &&
        pointer.x <= r.x + r.width &&
        pointer.y >= r.y &&
        pointer.y <= r.y + r.height,
    );
    if (!match) {
      destination = null;
      return;
    }
    const [id, r] = match;
    const dx = (pointer.x - r.x) / r.width,
      dy = (pointer.y - r.y) / r.height;
    const distances = [
      { edge: 'left' as Edge, n: dx },
      { edge: 'right' as Edge, n: 1 - dx },
      { edge: 'top' as Edge, n: dy },
      { edge: 'bottom' as Edge, n: 1 - dy },
    ].sort((a, b) => a.n - b.n);
    destination = { id, edge: distances[0].n < 0.25 ? distances[0].edge : 'center' };
  }
  setContext<WorkspaceContext>(WORKSPACE, {
    get rectangles() {
      return positions.panels;
    },
    get dragging() {
      return dragging;
    },
    get collapsed() {
      return layout.collapsed;
    },
    register(id, title) {
      titles = { ...titles, [id]: title };
      active = [...new Set([...active, id])];
      return () => {
        active = active.filter((key) => key !== id);
      };
    },
    begin(id, x, y, title) {
      dragging = id;
      dragTitle = title;
      move(x, y);
    },
    drag: move,
    drop() {
      if (dragging && destination) {
        layout = { ...layout, tree: dock(layout.tree, dragging, destination.id, destination.edge) };
        save();
      }
      cancel();
    },
    cancel,
    collapse(id) {
      if (layout.collapsed.includes(id))
        layout = { ...layout, collapsed: layout.collapsed.filter((key) => key !== id) };
      else if (active.filter((key) => !layout.collapsed.includes(key)).length > 1)
        layout = { ...layout, collapsed: [...layout.collapsed, id] };
      save();
    },
    keyboard(id, edge) {
      const ids = leaves(visible || layout.tree);
      const index = ids.indexOf(id);
      const target = ids[(index + 1) % ids.length];
      if (target) {
        layout = { ...layout, tree: dock(layout.tree, id, target, edge) };
        save();
      }
    },
  });
  let resizeOrigin: { divider: Divider; x: number; y: number } | null = null;
  function resize(e: PointerEvent) {
    if (!resizeOrigin) return;
    const { divider: d, x, y } = resizeOrigin;
    const delta = d.axis === 'x' ? e.clientX - x : e.clientY - y;
    const span = (d.axis === 'x' ? d.rect.width : d.rect.height) - 10;
    layout = { ...layout, tree: resizeSplit(layout.tree, d.id, d.ratio + delta / span) };
  }
  function finish() {
    resizeOrigin = null;
    resizing = null;
    save();
  }
</script>

{#if layout.collapsed.some((id) => active.includes(id))}
  <div class="workspace-tools">
    {#each layout.collapsed.filter((id) => active.includes(id)) as id}<button
        class="rounded border px-3 py-1 text-xs"
        onclick={() => {
          layout = { ...layout, collapsed: layout.collapsed.filter((key) => key !== id) };
          save();
        }}>+ {titles[id] || id}</button
      >{/each}
  </div>
{/if}
<div
  class="workspace"
  class:resizing
  bind:this={canvas}
  bind:clientWidth={width}
  bind:clientHeight={height}
>
  {@render children()}
  {#each positions.dividers as divider (divider.id)}
    {@const d = divider}
    <!-- svelte-ignore a11y_no_noninteractive_tabindex, a11y_no_noninteractive_element_interactions (A focusable separator implements the keyboard-resizable window-splitter pattern.) -->
    <div
      tabindex="0"
      class="divider"
      class:vertical={d.axis === 'x'}
      aria-label={`Resize ${d.axis === 'x' ? 'columns' : 'rows'}`}
      role="separator"
      aria-orientation={d.axis === 'x' ? 'vertical' : 'horizontal'}
      aria-valuenow={Math.round(d.ratio * 100)}
      aria-valuemin={10}
      aria-valuemax={90}
      style={d.axis === 'x'
        ? `left:${d.rect.x + (d.rect.width - 10) * d.ratio}px;top:${d.rect.y}px;width:10px;height:${d.rect.height}px`
        : `left:${d.rect.x}px;top:${d.rect.y + (d.rect.height - 10) * d.ratio}px;width:${d.rect.width}px;height:10px`}
      onpointerdown={(e) => {
        e.preventDefault();
        resizeOrigin = { divider: d, x: e.clientX, y: e.clientY };
        resizing = d.id;
        e.currentTarget.setPointerCapture(e.pointerId);
      }}
      onpointermove={resize}
      onpointerup={finish}
      onpointercancel={finish}
      onlostpointercapture={() => {
        if (resizeOrigin) finish();
      }}
      onkeydown={(e) => {
        if (['ArrowLeft', 'ArrowUp', 'ArrowRight', 'ArrowDown'].includes(e.key)) {
          e.preventDefault();
          layout = {
            ...layout,
            tree: resizeSplit(
              layout.tree,
              d.id,
              d.ratio + (['ArrowLeft', 'ArrowUp'].includes(e.key) ? -0.03 : 0.03),
            ),
          };
          save();
        }
      }}
    ></div>
  {/each}
  {#if preview}<div
      class="drop-preview"
      style={`left:${preview.x}px;top:${preview.y}px;width:${preview.width}px;height:${preview.height}px`}
    >
      <span>{destination?.edge === 'center' ? 'Swap' : `Dock ${destination?.edge}`}</span>
    </div>{/if}
  {#if dragging}<div class="drag-label" style={`left:${pointer.x + 14}px;top:${pointer.y + 14}px`}>
      {dragTitle}
    </div>{/if}
</div>

<style>
  .workspace-tools {
    flex-shrink: 0;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 12px;
    margin: 0 0 10px;
  }
  .workspace {
    position: relative;
    width: 100%;
    flex: 1;
    min-height: 0;
    overflow: hidden;
  }
  .divider {
    position: absolute;
    z-index: 15;
    cursor: row-resize;
    touch-action: none;
    border-radius: 6px;
  }
  .divider.vertical {
    cursor: col-resize;
  }
  .divider:hover,
  .divider:focus-visible {
    background: color-mix(in srgb, var(--primary) 45%, transparent);
    outline: none;
  }
  .drop-preview {
    position: absolute;
    z-index: 30;
    pointer-events: none;
    background: color-mix(in srgb, var(--primary) 18%, transparent);
    border: 2px solid var(--primary);
    border-radius: 12px;
    display: grid;
    place-items: center;
  }
  .drop-preview span,
  .drag-label {
    background: var(--background);
    border: 1px solid var(--primary);
    border-radius: 8px;
    padding: 8px 12px;
    font-size: 12px;
  }
  .drag-label {
    position: absolute;
    z-index: 40;
    pointer-events: none;
    box-shadow: 0 8px 24px #0004;
  }
</style>
