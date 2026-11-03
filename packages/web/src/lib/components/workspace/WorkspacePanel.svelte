<script lang="ts">
  import { getContext, onMount, type Snippet } from 'svelte';
  import { WORKSPACE, type WorkspaceContext } from './context';
  let {
    id,
    title,
    children,
    enabled = true,
    bodyClass = 'p-3',
    actions,
    busy = false,
  }: {
    id: string;
    title: string;
    children: Snippet;
    enabled?: boolean;
    /** Padding of the panel body; lists that run edge to edge pass an empty class. */
    bodyClass?: string;
    /** Controls shown in the panel header, such as display options for its content. */
    actions?: Snippet;
    /** The content is being refreshed: it stays readable, dimmed, under a progress bar. */
    busy?: boolean;
  } = $props();
  const workspace = getContext<WorkspaceContext | undefined>(WORKSPACE);
  const rect = $derived(workspace?.rectangles[id]);
  let origin: { x: number; y: number } | null = null;
  onMount(() => {
    if (enabled && workspace) return workspace.register(id, title);
  });
  function move(e: PointerEvent) {
    if (!origin || !workspace) return;
    if (!workspace.dragging && Math.hypot(e.clientX - origin.x, e.clientY - origin.y) < 5) return;
    if (!workspace.dragging) workspace.begin(id, e.clientX, e.clientY, title);
    workspace.drag(e.clientX, e.clientY);
  }
  function drop() {
    if (origin) {
      origin = null;
      workspace?.drop();
    }
  }
</script>

{#if enabled && workspace}
  <section
    aria-label={`${title} panel`}
    data-workspace-panel={id}
    class="panel"
    class:busy
    aria-busy={busy}
    class:dragging={workspace.dragging === id}
    hidden={!rect}
    style={rect
      ? `left:${rect.x}px;top:${rect.y}px;width:${rect.width}px;height:${rect.height}px`
      : ''}
  >
    <header class="relative flex shrink-0 items-center gap-2 border-b bg-muted/20 px-3 py-2">
      <span class="busy-bar" aria-hidden="true"></span>
      <button
        class="min-w-0 flex-1 cursor-grab truncate text-left text-sm font-medium"
        style="touch-action:none"
        aria-label={`Move ${title} panel`}
        onpointerdown={(e) => {
          e.preventDefault();
          origin = { x: e.clientX, y: e.clientY };
          e.currentTarget.setPointerCapture(e.pointerId);
        }}
        onpointermove={move}
        onpointerup={drop}
        onpointercancel={() => {
          origin = null;
          workspace.cancel();
        }}
        onlostpointercapture={() => {
          if (origin) {
            origin = null;
            workspace.cancel();
          }
        }}
        onkeydown={(e) => {
          const edges: Record<string, 'left' | 'right' | 'top' | 'bottom' | 'center'> = {
            ArrowLeft: 'left',
            ArrowRight: 'right',
            ArrowUp: 'top',
            ArrowDown: 'bottom',
            Enter: 'center',
          };
          if (edges[e.key]) {
            e.preventDefault();
            workspace.keyboard(id, edges[e.key]);
          }
        }}>{title}</button
      >
      {#if actions}<div class="flex shrink-0 items-center">{@render actions()}</div>{/if}
      <button
        class="panel-action rounded px-2 hover:bg-muted"
        aria-label={`Collapse ${title} panel`}
        onclick={() => workspace.collapse(id)}>−</button
      >
    </header>
    <div class={`panel-body min-h-0 flex-1 overflow-auto ${bodyClass}`}>{@render children()}</div>
  </section>
{:else}
  <div class="relative" class:busy aria-busy={busy}>
    <span class="busy-bar" aria-hidden="true"></span>
    <div class="busy-body">{@render children()}</div>
  </div>
{/if}

<style>
  .panel {
    position: absolute;
    display: flex;
    flex-direction: column;
    min-width: 0;
    min-height: 0;
    border: 1px solid var(--border);
    border-radius: 12px;
    background: var(--background);
    overflow: hidden;
  }
  .panel[hidden] {
    display: none;
  }
  .dragging {
    opacity: 0.45;
  }
  .panel-action {
    opacity: 0;
    transition: opacity 120ms;
  }
  header:hover .panel-action,
  header:focus-within .panel-action {
    opacity: 1;
  }
  /* Busy: shown only after 150 ms, so cached responses do not flicker. */
  .busy-bar {
    position: absolute;
    left: 0;
    right: 0;
    bottom: -1px;
    height: 2px;
    overflow: hidden;
    opacity: 0;
    pointer-events: none;
  }
  .relative > .busy-bar {
    top: 0;
    bottom: auto;
  }
  .busy-bar::before {
    content: '';
    position: absolute;
    inset: 0;
    width: 40%;
    background: var(--primary);
    border-radius: 2px;
    animation: busy-slide 1.1s ease-in-out infinite;
  }
  .busy .busy-bar {
    animation: busy-show 0s linear 150ms forwards;
  }
  .panel-body,
  .busy-body {
    transition: opacity 150ms ease;
  }
  .busy .panel-body,
  .busy .busy-body {
    opacity: 0.6;
    transition-delay: 150ms;
  }
  @keyframes busy-slide {
    from {
      transform: translateX(-100%);
    }
    to {
      transform: translateX(250%);
    }
  }
  @keyframes busy-show {
    to {
      opacity: 1;
    }
  }
  @media (prefers-reduced-motion: reduce) {
    .busy-bar::before {
      width: 100%;
      opacity: 0.6;
      animation: none;
    }
  }
  .panel-body :global(h4) {
    display: none;
  }
  .panel-body :global(.max-h-64) {
    max-height: none;
    overflow: visible;
  }
</style>
