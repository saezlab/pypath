<script lang="ts">
  import { ChevronDown, ChevronRight, Loader2 } from '@lucide/svelte';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import OntologyHierarchyNode from './OntologyHierarchyNode.svelte';
  import type { HierarchyNode } from './ontology-hierarchy-types';

  interface Props {
    node: HierarchyNode;
    seedTermId?: string | null;
    selectedNodeId?: string | null;
    depth?: number;
    onExpand: (node: HierarchyNode) => Promise<void>;
    onSelect: (node: HierarchyNode) => void;
  }

  let {
    node,
    seedTermId = null,
    selectedNodeId = null,
    depth = 0,
    onExpand,
    onSelect,
  }: Props = $props();

  const isOpen = $derived(node.open ?? depth < 3);
  const isSeedTerm = $derived(node.id === seedTermId);
  const isSelectedNode = $derived(node.id === selectedNodeId);
  const hasChildren = $derived((node.children?.length ?? 0) > 0);
  const canExpand = $derived(hasChildren || !node.childrenLoaded);

  async function activateNode(event?: MouseEvent | KeyboardEvent) {
    event?.stopPropagation();
    const wasSelected = isSelectedNode;
    onSelect(node);
    if (isOpen && wasSelected) {
      node.open = false;
      return;
    }
    if (!node.childrenLoaded) await onExpand(node);
    node.open = true;
  }

  async function toggleOpen(event?: MouseEvent) {
    await activateNode(event);
  }

  function handleRowKeydown(event: KeyboardEvent) {
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      void activateNode(event);
    }
  }
</script>

<div class="w-max min-w-full select-none" style={`padding-left: ${depth * 0.35}rem`}>
  <div
    class="group flex w-max cursor-pointer items-center gap-2 rounded-md px-2 py-1.5 hover:bg-muted/60 {isSelectedNode
      ? 'bg-primary/10 ring-1 ring-primary/30'
      : isSeedTerm
        ? 'bg-primary/5 ring-1 ring-primary/15'
        : ''}"
    role="treeitem"
    tabindex="0"
    aria-selected={isSelectedNode}
    onclick={activateNode}
    onkeydown={handleRowKeydown}
  >
    <button
      type="button"
      class="flex h-6 w-6 shrink-0 items-center justify-center rounded hover:bg-background disabled:cursor-default disabled:opacity-40"
      disabled={!canExpand || node.childrenLoading}
      onclick={toggleOpen}
      aria-label={isOpen ? 'Collapse term' : 'Expand term'}
    >
      {#if node.childrenLoading}
        <Loader2 class="h-4 w-4 animate-spin text-muted-foreground" />
      {:else if canExpand && isOpen}
        <ChevronDown class="h-4 w-4 text-muted-foreground" />
      {:else if canExpand}
        <ChevronRight class="h-4 w-4 text-muted-foreground" />
      {/if}
    </button>

    <div class="min-w-0 flex-1 pr-4">
      <div class="flex min-w-0 flex-wrap items-center gap-2">
        <span
          class="whitespace-nowrap text-sm font-medium {isSelectedNode || isSeedTerm
            ? 'text-primary'
            : ''}">{node.name || node.id}</span
        >
        <Badge variant="outline" class="shrink-0 font-mono text-[10px]">{node.id}</Badge>
        {#if isSeedTerm}
          <Badge class="text-[10px]">Selected term</Badge>
        {/if}
        {#if node.childrenLoaded}
          <span class="text-xs text-muted-foreground"
            >{(node.children || []).length} of {node.childrenTotal ?? (node.children || []).length} children</span
          >
        {/if}
      </div>
      {#if node.error}
        <div class="mt-1 text-xs text-destructive">{node.error}</div>
      {/if}
    </div>
  </div>

  {#if isOpen && hasChildren}
    <div class="ml-3 w-max min-w-full border-l border-border/60">
      {#each node.children || [] as child (child.id)}
        <OntologyHierarchyNode
          node={child}
          {seedTermId}
          {selectedNodeId}
          depth={depth + 1}
          {onExpand}
          {onSelect}
        />
      {/each}
      {#if node.childrenNextCursor != null}
        <button
          type="button"
          class="m-2 rounded border px-3 py-1 text-sm"
          disabled={node.childrenLoading}
          onclick={() => onExpand(node)}
        >
          {node.childrenLoading ? 'Loading children…' : 'Load more children'}
        </button>
      {/if}
    </div>
  {/if}
</div>
