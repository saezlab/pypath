<script lang="ts">
  import EntityResultsList from '$lib/components/entity/EntityResultsList.svelte';
  import { Filter } from '@lucide/svelte';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Card, CardContent } from '$lib/components/ui/card/index.js';
  import { Checkbox } from '$lib/components/ui/checkbox/index.js';
  import { Label } from '$lib/components/ui/label/index.js';
  import {
    Sheet,
    SheetContent,
    SheetHeader,
    SheetTitle,
    SheetTrigger,
  } from '$lib/components/ui/sheet/index.js';
  import EntityResultCard from '$lib/components/entity/EntityResultCard.svelte';
  import ExplorerWorkspace from '$lib/components/workspace/ExplorerWorkspace.svelte';
  import WorkspacePanel from '$lib/components/workspace/WorkspacePanel.svelte';
  import type { EntityWithIdentifiers } from '$lib/types/entities';
  import { IsMobile } from '$lib/hooks/is-mobile.svelte';
  import {
    fetchOntologySearch,
    fetchScopedOntologySearch,
    fetchScopedOntologyIdCounts,
    type SelectionScopeRequest,
  } from '$lib/api/client';

  import { createPagedQuery, type PageState } from '$lib/features/explorer/paged-query';
  import { page as currentPage } from '$app/state';
  import { untrack } from 'svelte';

  import type { SearchFilters } from '$lib/types/search';

  interface Props {
    query: string;
    filters: SearchFilters;
    onFiltersChange: (filters: SearchFilters) => void;
    selectedEntityIds?: string[];
    selectedEntityPks?: Array<string | number>;
    selectedAnnotationIds?: string[];
    selectionScope?: SelectionScopeRequest;
  }

  let {
    query,
    filters,
    onFiltersChange,
    selectedEntityPks,
    selectedAnnotationIds,
    selectionScope,
  }: Props = $props();

  const isMobile = new IsMobile();
  const RESULTS_PER_PAGE = 30;

  type Term = Awaited<ReturnType<typeof fetchOntologySearch>>[number];
  let pageState = $state<PageState<Term>>({
    items: [],
    loading: true,
    loadingMore: false,
    hasMore: true,
    error: null,
  });
  const paging = createPagedQuery<Term>(RESULTS_PER_PAGE, (next) => {
    pageState = next;
  });
  const results = $derived(pageState.items);
  const loading = $derived(pageState.loading);
  const loadingMore = $derived(pageState.loadingMore);
  const hasMore = $derived(pageState.hasMore);
  const error = $derived(pageState.error);
  const loadMore = paging.loadMore;
  let ontologyOptions = $state<Array<{ value: string; count: number }>>([]);
  let loadingOntologies = $state(true);

  const isScoped = $derived(
    !!(selectionScope || selectedEntityPks?.length || selectedAnnotationIds?.length),
  );
  const selectedOntologyIds = $derived(filters.ontology_ids || []);

  // Fetch ontology counts from terms associated with entities in the current scope.
  $effect(() => {
    const q = query;
    const scope = selectionScope;
    const ePks = selectedEntityPks || [];
    const tIds = selectedAnnotationIds || [];

    const controller = new AbortController();
    void currentPage.data.selectedRelease;
    loadingOntologies = true;
    fetchScopedOntologyIdCounts(
      {
        entityPks: scope ? undefined : ePks.length > 0 ? ePks : undefined,
        annotationTermIds: scope ? undefined : tIds.length > 0 ? tIds : undefined,
        selectionScope: scope,
        query: q || undefined,
      },
      controller.signal,
    )
      .then((counts) => {
        if (controller.signal.aborted) return;
        ontologyOptions = counts.map((c) => ({ value: c.ontologyId, count: c.scopedCount }));
      })
      .catch(() => {
        if (!controller.signal.aborted) ontologyOptions = [];
      })
      .finally(() => {
        if (!controller.signal.aborted) loadingOntologies = false;
      });
    return () => controller.abort();
  });

  // Snapshot the complete query, including release, for every initial/page request.
  const requestKey = $derived(
    JSON.stringify([
      currentPage.data.selectedRelease,
      query,
      selectedOntologyIds,
      isScoped,
      selectionScope,
      selectedEntityPks || [],
      selectedAnnotationIds || [],
    ]),
  );
  $effect(() => {
    const [, q, ontologyIds, scoped, scope, entityPks, termIds] = JSON.parse(requestKey);
    untrack(() => {
      void paging.reset((offset, signal) =>
        scoped
          ? fetchScopedOntologySearch(
              {
                entityPks: scope ? undefined : entityPks,
                termIds: scope ? undefined : termIds,
                selectionScope: scope,
                query: q,
                ontologyIds,
                limit: RESULTS_PER_PAGE,
                offset,
              },
              signal,
            )
          : fetchOntologySearch({ query: q, ontologyIds, limit: RESULTS_PER_PAGE, offset }, signal),
      );
    });
    return paging.cancel;
  });

  function toggleOntologyId(ontologyId: string) {
    const next = selectedOntologyIds.includes(ontologyId)
      ? selectedOntologyIds.filter((item) => item !== ontologyId)
      : [...selectedOntologyIds, ontologyId];
    onFiltersChange({
      ...filters,
      ontology_ids: next.length > 0 ? next : undefined,
    });
  }

  function handleClearFilters() {
    onFiltersChange({
      ...filters,
      ontology_ids: undefined,
    });
  }

  function termEntity(term: (typeof results)[number]): EntityWithIdentifiers {
    return {
      entityPk: term.entityPk,
      canonicalIdentifier: term.termId,
      canonicalIdentifierType: term.ontologyId || term.ontologyPrefix || '',
      label: term.label,
      entityType: 'ontology_class',
      taxonomyId: null,
      entityAttributes: {},
      identifiers: [],
      sources: term.sources,
      relationCount: term.annotatedRelationCount,
      ontologyHierarchy: {
        termId: term.termId,
        ontologyPrefix: term.ontologyPrefix,
        ontologyId: term.ontologyId,
        label: term.label,
        definition: term.definition,
        childCount: 0,
        parentCount: 0,
      },
    };
  }
</script>

{#snippet filterSidebarContent()}
  <div class="space-y-6">
    <div class="mb-3 flex items-center justify-between gap-2">
      <h4 class="text-[11px] font-semibold uppercase tracking-[0.16em] text-muted-foreground">
        Ontologies
      </h4>
      {#if !isMobile.current && selectedOntologyIds.length > 0}
        <Button
          variant="ghost"
          size="xs"
          onclick={handleClearFilters}
          class="h-auto px-0 text-[11px] text-muted-foreground"
        >
          Clear all
        </Button>
      {/if}
    </div>
    <div class="max-h-[calc(100vh-14rem)] space-y-1 overflow-y-auto pr-2">
      {#each ontologyOptions as ontology}
        {@const selected = selectedOntologyIds.includes(ontology.value)}
        <div class="flex items-center justify-between gap-2 py-0.5">
          <Label
            class="flex min-w-0 flex-1 cursor-pointer items-center gap-2 text-sm font-normal leading-5 text-foreground {selected
              ? 'font-medium'
              : ''}"
          >
            <Checkbox
              checked={selected}
              onCheckedChange={() => toggleOntologyId(ontology.value)}
              class={selected ? 'h-4 w-4 flex-shrink-0 border-primary' : 'h-4 w-4 flex-shrink-0'}
            />
            <span class="truncate">{ontology.value}</span>
          </Label>
          <span class="text-xs text-muted-foreground tabular-nums flex-shrink-0"
            >{ontology.count.toLocaleString()}</span
          >
        </div>
      {:else}
        {#if loadingOntologies}
          <p class="text-sm text-muted-foreground">Loading filters...</p>
        {:else}
          <p class="text-sm text-muted-foreground">No ontology filters available</p>
        {/if}
      {/each}
    </div>
  </div>
{/snippet}

{#snippet resultsPane()}
  <div>
    {#if error}
      <div
        class="m-4 rounded-lg border border-destructive/30 bg-destructive/10 p-3 text-sm text-destructive"
      >
        {error}
      </div>
    {/if}

    {#if loading && results.length === 0}
      <EntityResultsList>
        {#each Array.from({ length: 6 }) as _, _i}
          <div class="h-12 animate-pulse bg-muted/30"></div>
        {/each}
      </EntityResultsList>
    {:else if results.length > 0}
      <div class="space-y-4">
        <EntityResultsList>
          {#each results as term}
            <EntityResultCard result={termEntity(term)} />
          {/each}
        </EntityResultsList>

        {#if hasMore}
          <div class="flex justify-center py-4" style="min-height: 40px;">
            <Button variant="outline" onclick={loadMore} disabled={loadingMore}>
              {#if loadingMore}
                <div
                  class="h-4 w-4 animate-spin rounded-full border-2 border-primary border-t-transparent mr-2"
                ></div>
                <span>Loading...</span>
              {:else}
                <span>Load more</span>
              {/if}
            </Button>
          </div>
        {/if}
      </div>
    {:else}
      <Card class="m-4 border-dashed">
        <CardContent class="flex flex-col items-center justify-center gap-2 py-16 text-center">
          <div class="text-lg font-semibold">No ontology terms found</div>
          <p class="max-w-2xl text-sm text-muted-foreground">
            {#if query.trim().length > 0}
              {#if isScoped}
                Try a different ontology term, synonym, or ID within the current entity scope.
              {:else}
                Try a different ontology term, synonym, or ID such as GO:0005634, KW-0001, or
                MI:0217.
              {/if}
            {:else}
              {#if isScoped}
                No associated ontology terms are available for the current entity scope.
              {:else}
                No ontology terms are available to browse right now.
              {/if}
            {/if}
          </p>
        </CardContent>
      </Card>
    {/if}
  </div>
{/snippet}

{#if isMobile.current}
  <div class="flex h-full min-h-0 flex-col overflow-hidden">
    <div class="border-b p-4">
      <Sheet>
        <SheetTrigger>
          {#snippet child({ props })}
            <Button {...props} variant="outline" class="w-full">
              <Filter class="mr-2 size-4" />
              Filters
              {#if selectedOntologyIds.length > 0}
                <Badge variant="secondary" class="ml-2">{selectedOntologyIds.length}</Badge>
              {/if}
            </Button>
          {/snippet}
        </SheetTrigger>
        <SheetContent side="left" class="w-[85%] overflow-y-auto sm:w-[400px]">
          <SheetHeader>
            <SheetTitle>Ontology filters</SheetTitle>
          </SheetHeader>
          <div class="pt-4">{@render filterSidebarContent()}</div>
        </SheetContent>
      </Sheet>
    </div>
    <div class="min-h-0 flex-1">
      {@render resultsPane()}
    </div>
  </div>
{:else}
  <ExplorerWorkspace name="ontology" panels={['results', 'filters']}>
    <WorkspacePanel id="results" title="Results" bodyClass=""
      >{@render resultsPane()}</WorkspacePanel
    >
    <WorkspacePanel id="filters" title="Filters">{@render filterSidebarContent()}</WorkspacePanel>
  </ExplorerWorkspace>
{/if}
