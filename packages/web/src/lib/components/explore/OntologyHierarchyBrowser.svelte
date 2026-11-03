<script lang="ts">
  import type { WireOntologyChildren } from '$lib/api/contracts';
  import { releaseFetch } from '$lib/api/release';

  import { AlertCircle, Check, Link, Loader2, Plus } from '@lucide/svelte';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { getSelectionStore } from '$lib/stores/selection.svelte';
  import OntologyHierarchyNode from './OntologyHierarchyNode.svelte';
  import type { ApiTermInfo, HierarchyNode } from './ontology-hierarchy-types';

  export type OntologyBrowserTerm = {
    termId: string;
    ontologyPrefix: string | null;
    label: string | null;
    definition: string | null;
    ontologyId?: string | null;
  };

  interface Props {
    term: OntologyBrowserTerm | null;
  }

  let { term }: Props = $props();

  const selection = getSelectionStore();

  let root = $state<HierarchyNode | null>(null);
  let loading = $state(false);
  let error = $state<string | null>(null);
  let loadedTermId = $state<string | null>(null);
  let selectedNode = $state<HierarchyNode | null>(null);
  let selectedDetails = $state<ApiTermInfo | null>(null);
  let detailsLoading = $state(false);
  let detailsError = $state<string | null>(null);

  const selectedNodeInSelection = $derived(
    selectedNode ? selection.isAnnotationSelected(selectedNode.id) : false,
  );

  async function fetchJson<T>(url: string, init?: RequestInit): Promise<T> {
    const res = await releaseFetch(url, init);
    if (!res.ok) {
      const text = await res.text().catch(() => '');
      throw new Error(text || `Request failed with status ${res.status}`);
    }
    return res.json() as Promise<T>;
  }

  function findNode(node: HierarchyNode | null, termId: string): HierarchyNode | null {
    if (!node) return null;
    if (node.id === termId) return node;
    for (const child of node.children || []) {
      const match = findNode(child, termId);
      if (match) return match;
    }
    return null;
  }

  async function ontologyIdForTermId(termId: string): Promise<string | null> {
    const node = findNode(root, termId);
    if (node?.ontologyId) return node.ontologyId;
    if (term?.termId === termId && term.ontologyId) return term.ontologyId;
    if (selectedDetails?.id === termId && selectedDetails.ontologyId)
      return selectedDetails.ontologyId;
    return null;
  }

  function cloneNode(node: HierarchyNode, selectedTermId: string): HierarchyNode {
    return {
      id: node.id,
      name: node.name || node.id,
      ontologyId: node.ontologyId,
      distance: node.distance,
      children: (node.children || []).map((child) => cloneNode(child, selectedTermId)),
      open: node.id !== selectedTermId,
      childrenLoaded: false,
      childrenLoading: false,
      error: null,
    };
  }

  async function fetchTermTree(termId: string): Promise<HierarchyNode | null> {
    const data = await fetchJson<{ root: HierarchyNode | null }>('/app-api/ontology/tree', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ term_ids: [termId], ontologyId: term?.ontologyId || null }),
    });
    return data.root;
  }

  async function fetchApiTerms(termIds: string[]): Promise<Record<string, ApiTermInfo | null>> {
    if (termIds.length === 0) return {};
    const data = await fetchJson<{ terms: Record<string, ApiTermInfo | null> }>('/app-api/terms', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ term_ids: termIds }),
    });
    return data.terms || {};
  }

  function mergeChildren(existing: HierarchyNode[], fetched: HierarchyNode[]): HierarchyNode[] {
    const byId = new Map<string, HierarchyNode>();
    for (const child of existing) byId.set(child.id, child);
    for (const child of fetched) {
      const previous = byId.get(child.id);
      byId.set(
        child.id,
        previous ? { ...child, ...previous, name: previous.name || child.name } : child,
      );
    }
    return Array.from(byId.values()).sort((a, b) => (a.name || a.id).localeCompare(b.name || b.id));
  }

  async function selectNode(node: HierarchyNode) {
    selectedNode = node;
    selectedDetails = {
      id: node.id,
      name: node.name || node.id,
      definition: null,
      namespace: node.id.includes(':') ? node.id.split(':', 1)[0].toLowerCase() : null,
      relationCount: 0,
    };
    detailsLoading = true;
    detailsError = null;
    try {
      const terms = await fetchApiTerms([node.id]);
      const info = terms[node.id];
      if (selectedNode?.id === node.id && info) selectedDetails = info;
    } catch (err) {
      if (selectedNode?.id === node.id)
        detailsError = err instanceof Error ? err.message : 'Failed to load term details';
    } finally {
      if (selectedNode?.id === node.id) detailsLoading = false;
    }
  }

  async function expandNode(node: HierarchyNode) {
    if (node.childrenLoading || (node.childrenLoaded && node.childrenNextCursor == null)) return;
    const ontologyId = await ontologyIdForTermId(node.id);

    node.childrenLoading = true;
    node.error = null;
    try {
      const params = new URLSearchParams({
        termId: node.id,
        limit: '100',
        offset: String(node.childrenNextCursor ?? 0),
      });
      if (ontologyId) params.set('ontologyId', ontologyId);
      const data = await fetchJson<WireOntologyChildren>(
        `/app-api/ontology/children?${params.toString()}`,
      );
      const fetchedChildren = (data.children || []).map(
        (child) =>
          ({
            id: child.termId || child.id,
            name: child.label || child.name || child.termId || child.id,
            ontologyId: child.ontologyId || ontologyId || node.ontologyId || null,
            children: [],
            open: false,
            childrenLoaded: false,
            childrenLoading: false,
            error: null,
          }) satisfies HierarchyNode,
      );
      node.children = mergeChildren(node.children || [], fetchedChildren);
      node.childrenLoaded = true;
      node.childrenNextCursor = data.nextCursor;
      node.childrenTotal = data.total;
    } catch (err) {
      node.error = err instanceof Error ? err.message : 'Failed to load children';
    } finally {
      node.childrenLoading = false;
    }
  }

  async function loadInitialTree(termId: string) {
    loading = true;
    error = null;
    root = null;
    try {
      const tree = await fetchTermTree(termId);
      root = tree
        ? cloneNode(tree, termId)
        : {
            id: termId,
            name: term?.label || termId,
            ontologyId: term?.ontologyId || null,
            children: [],
            open: true,
            childrenLoaded: false,
          };
      loadedTermId = termId;
      const initialNode = findNode(root, termId) || root;
      await selectNode(initialNode);
      await expandNode(initialNode);
      initialNode.open = true;
    } catch (err) {
      error = err instanceof Error ? err.message : 'Failed to load ontology hierarchy';
      root = {
        id: termId,
        name: term?.label || termId,
        ontologyId: term?.ontologyId || null,
        children: [],
        open: true,
        childrenLoaded: false,
      };
      await selectNode(root);
      await expandNode(root);
      root.open = true;
    } finally {
      loading = false;
    }
  }

  function addSelectedNode() {
    if (!selectedNode) return;
    selection.addAnnotation({
      id: selectedNode.id,
      label: selectedDetails?.name || selectedNode.name || selectedNode.id,
      namespace:
        selectedDetails?.namespace ||
        (selectedNode.id.includes(':')
          ? selectedNode.id.split(':', 1)[0].toLowerCase()
          : undefined),
      definition: selectedDetails?.definition,
    });
  }

  $effect(() => {
    const termId = term?.termId || null;
    if (!termId || loadedTermId === termId) return;
    void loadInitialTree(termId);
  });
</script>

<div
  class="grid min-h-0 grid-cols-1 overflow-hidden rounded-lg border bg-background md:grid-cols-[minmax(0,1fr)_minmax(18rem,32%)]"
>
  <div class="min-h-[24rem] overflow-auto p-3">
    {#if loading}
      <div class="flex items-center justify-center gap-2 py-16 text-sm text-muted-foreground">
        <Loader2 class="h-4 w-4 animate-spin" />
        Loading hierarchy...
      </div>
    {:else if error}
      <div
        class="mb-4 flex items-start gap-2 rounded-lg border border-destructive/30 bg-destructive/10 p-3 text-sm text-destructive"
      >
        <AlertCircle class="mt-0.5 h-4 w-4 shrink-0" />
        <div>
          <div class="font-medium">Could not load initial tree</div>
          <div class="opacity-90">{error}</div>
        </div>
      </div>
    {/if}

    {#if root && term}
      <div class="inline-block min-w-full bg-background p-2">
        <OntologyHierarchyNode
          node={root}
          seedTermId={term.termId}
          selectedNodeId={selectedNode?.id || null}
          onExpand={expandNode}
          onSelect={selectNode}
        />
      </div>
    {:else if !loading}
      <div class="rounded-lg border border-dashed py-12 text-center text-sm text-muted-foreground">
        No hierarchy available.
      </div>
    {/if}
  </div>

  <aside class="min-h-0 min-w-0 overflow-y-auto border-t bg-muted/10 p-4 md:border-l md:border-t-0">
    {#if selectedNode}
      <div class="space-y-4">
        <div class="min-w-0 space-y-2">
          <div class="flex flex-wrap items-center gap-2">
            <Badge variant="outline" class="font-mono">{selectedNode.id}</Badge>
            {#if selectedDetails?.namespace}
              <Badge variant="outline">{selectedDetails.namespace}</Badge>
            {/if}
            {#if selectedNode.id === term?.termId}
              <Badge>Selected term</Badge>
            {/if}
          </div>
          <h3 class="text-lg font-semibold leading-tight">
            {selectedDetails?.name || selectedNode.name || selectedNode.id}
          </h3>
          {#if (selectedDetails?.relationCount || 0) > 0}
            <div
              class="inline-flex items-center gap-1.5 rounded-md border px-2 py-1 text-sm text-muted-foreground"
            >
              <Link class="size-4" />
              <span class="tabular-nums"
                >{(selectedDetails?.relationCount || 0).toLocaleString()}</span
              >
              <span>relations</span>
            </div>
          {/if}
          {#if detailsLoading}
            <div class="flex items-center gap-2 text-sm text-muted-foreground">
              <Loader2 class="h-4 w-4 animate-spin" />
              Loading term details...
            </div>
          {/if}
          {#if detailsError}
            <div
              class="rounded-md border border-destructive/30 bg-destructive/10 p-2 text-sm text-destructive"
            >
              {detailsError}
            </div>
          {/if}
        </div>

        <Button
          size="sm"
          variant={selectedNodeInSelection ? 'default' : 'outline'}
          onclick={addSelectedNode}
          disabled={selectedNodeInSelection}
        >
          {#if selectedNodeInSelection}
            <Check class="mr-1 h-4 w-4" />
            Added to selection
          {:else}
            <Plus class="mr-1 h-4 w-4" />
            Add this term
          {/if}
        </Button>

        <section class="min-w-0 space-y-2">
          <h4 class="text-xs font-semibold uppercase tracking-wide text-muted-foreground">
            Definition
          </h4>
          <p
            class="whitespace-pre-wrap rounded-lg bg-background p-3 text-sm leading-6 text-foreground/90"
          >
            {selectedDetails?.definition || 'No definition available for this ontology term.'}
          </p>
        </section>
      </div>
    {/if}
  </aside>
</div>
