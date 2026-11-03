<script lang="ts">
  import ExportMenu from '$lib/components/explore/ExportMenu.svelte';
  import ExplorerWorkspace from '$lib/components/workspace/ExplorerWorkspace.svelte';
  import WorkspacePanel from '$lib/components/workspace/WorkspacePanel.svelte';
  import { untrack } from 'svelte';
  import { page as currentPage } from '$app/state';
  import { Filter, X } from '@lucide/svelte';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Alert, AlertDescription } from '$lib/components/ui/alert/index.js';
  import {
    Sheet,
    SheetContent,
    SheetHeader,
    SheetTitle,
    SheetTrigger,
  } from '$lib/components/ui/sheet/index.js';
  import EntityDetailsDialog from '$lib/components/entity/EntityDetailsDialog.svelte';
  import InteractionFilterSidebar from '$lib/components/interactions/InteractionFilterSidebar.svelte';
  import InteractionDetailsSheet from '$lib/components/interactions/InteractionDetailsSheet.svelte';
  import RelationsTable from '$lib/components/interactions/RelationsTable.svelte';
  import { IsMobile } from '$lib/hooks/is-mobile.svelte';
  import { getUiPreferences } from '$lib/stores/ui-preferences.svelte';
  import { fetchRelationsSearch, fetchTopHit, type TopHit } from '$lib/api/client';
  import { formatNumber } from '$lib/utils/format';
  import type { EntityLike } from '$lib/domain/display';
  import type { SearchFilters } from '$lib/types/search';
  import type { InteractionListRow } from '$lib/types/interactions';

  interface Props {
    filters: SearchFilters;
    onFilterChange: (filters: SearchFilters) => void;
    query?: string;
    entitySearchFilters?: SearchFilters;
    scopedEntityIds?: string[];
    scopedAnnotationIds?: string[];
    scopeEndpointMode?: 'any' | 'both';
    scopeMode?: 'union' | 'intersection';
  }

  let {
    filters,
    onFilterChange,
    query = '',
    entitySearchFilters = {},
    scopedEntityIds,
    scopedAnnotationIds,
    scopeEndpointMode = 'any',
    scopeMode = 'union',
  }: Props = $props();

  const isMobile = new IsMobile();
  const RESULTS_PER_PAGE = 20;

  let results = $state<InteractionListRow[]>([]);
  let loading = $state(true);
  let loadingMore = $state(false);
  let hasMore = $state(true);
  // The relations matching the current scope and filters (the filter panels show totals).
  let total = $state<number | null>(null);
  let offset = $state(0);
  let error = $state<string | null>(null);
  let selectedInteraction = $state<InteractionListRow | null>(null);
  let detailsOpen = $state(false);
  let detailsEntity = $state<EntityLike | null>(null);
  let entityDetailsOpen = $state(false);
  let queryEntityIds = $state<string[]>([]);
  // The entity (or group) a text search resolved to; its relations are shown.
  let topHit = $state<TopHit | null>(null);
  let queryEntityIdsLoading = $state(false);
  let queryEntityIdsError = $state<string | null>(null);

  function rowsFromSearch(
    data: Awaited<ReturnType<typeof fetchRelationsSearch>>,
  ): InteractionListRow[] {
    if (data.rows?.length) {
      return data.rows.filter((row) => row.relation && row.subjectEntity && row.objectEntity);
    }
    return [];
  }

  const activeFilterCount = $derived(
    Object.entries(filters).reduce((count, [, value]) => {
      if (Array.isArray(value)) return count + value.length;
      if (value !== null && value !== undefined) return count + 1;
      return count;
    }, 0),
  );

  const effectiveFilters = $derived({
    ...filters,
    ...(query.trim() && queryEntityIds.length > 0 ? { entity_ids: queryEntityIds } : {}),
    ...(scopedEntityIds && scopedEntityIds.length > 0 ? { scope_entity_ids: scopedEntityIds } : {}),
    ...(scopedEntityIds && scopedEntityIds.length > 0
      ? { scope_endpoint_mode: scopeEndpointMode }
      : {}),
    ...(scopedAnnotationIds && scopedAnnotationIds.length > 0
      ? { scope_annotation_ids: scopedAnnotationIds }
      : {}),
    ...(scopedEntityIds?.length || scopedAnnotationIds?.length
      ? { selection_scope_mode: scopeMode }
      : {}),
  });
  // An export names what it contains: filters, or entity and group keys (a text search
  // contributes its top hit's key), never the whole collection.
  // The filters that select something, so a copied API command lists only those.
  const exportFilters = $derived(
    Object.fromEntries(
      Object.entries(effectiveFilters).filter(([, value]) =>
        Array.isArray(value) ? value.length > 0 : value != null && value !== '',
      ),
    ),
  );
  const exportable = $derived(!queryEntityIdsLoading && Object.keys(exportFilters).length > 0);

  const preferences = getUiPreferences();
  const queryKey = $derived(
    JSON.stringify([
      currentPage.data.selectedRelease,
      query.trim(),
      entitySearchFilters,
      preferences.groupResults,
    ]),
  );
  // A text search shows the relations of its top hit: the entity (or group) the entity
  // search finds first, not of everything that matches the text.
  $effect(() => {
    const [, q, searchFilters, grouped] = JSON.parse(queryKey);
    if (!q) {
      queryEntityIds = [];
      topHit = null;
      queryEntityIdsLoading = false;
      return;
    }

    let cancelled = false;
    const controller = new AbortController();
    queryEntityIds = [];
    topHit = null;
    queryEntityIdsError = null;
    queryEntityIdsLoading = true;
    untrack(() => fetchTopHit(q, searchFilters, grouped, controller.signal))
      .then((hit) => {
        if (cancelled) return;
        topHit = hit;
        queryEntityIds = hit ? [hit.key] : [];
      })
      .catch((err) => {
        if (!cancelled) {
          queryEntityIds = [];
          queryEntityIdsError =
            err instanceof Error ? err.message : 'Failed to resolve entity search';
        }
      })
      .finally(() => {
        if (!cancelled) queryEntityIdsLoading = false;
      });
    return () => {
      cancelled = true;
      controller.abort();
    };
  });

  let searchGeneration = 0;
  const searchKey = $derived(
    JSON.stringify([
      currentPage.data.selectedRelease,
      effectiveFilters,
      query.trim(),
      queryEntityIdsLoading,
      queryEntityIds.length,
      queryEntityIdsError,
    ]),
  );
  $effect(() => {
    const current = ++searchGeneration;
    const [, f, q, resolvingQueryEntities, matchCount, resolutionError] = JSON.parse(searchKey);
    const queryHasNoEntityMatches = !!q && !resolvingQueryEntities && matchCount === 0;
    const controller = new AbortController();

    loadingMore = false;
    offset = 0;
    hasMore = true;
    loading = true;
    error = null;
    total = null;

    if (resolvingQueryEntities) {
      return;
    }

    if (queryHasNoEntityMatches) {
      error = resolutionError;
      results = [];
      loading = false;
      hasMore = false;
      return;
    }

    (async () => {
      try {
        const data = await untrack(() =>
          fetchRelationsSearch(
            { filters: f, limit: RESULTS_PER_PAGE, offset: 0 },
            controller.signal,
          ),
        );
        if (current !== searchGeneration) return;
        const hits = rowsFromSearch(data);

        results = hits;
        total = typeof data.total === 'number' ? data.total : hits.length;
        hasMore = hits.length === RESULTS_PER_PAGE;
        offset = RESULTS_PER_PAGE;
      } catch (err) {
        if (current !== searchGeneration) return;
        error = err instanceof Error ? err.message : 'Failed to load relations';
        hasMore = false;
      } finally {
        if (current === searchGeneration) loading = false;
      }
    })();
    return () => {
      searchGeneration++;
      controller.abort();
    };
  });

  async function loadMore() {
    if (loadingMore || !hasMore) return;
    const current = searchGeneration;
    loadingMore = true;
    try {
      const data = await fetchRelationsSearch({
        filters: effectiveFilters,
        limit: RESULTS_PER_PAGE,
        offset,
      });
      if (current !== searchGeneration) return;
      const hits = rowsFromSearch(data);

      results = [...results, ...hits];
      hasMore = hits.length === RESULTS_PER_PAGE;
      offset += RESULTS_PER_PAGE;
    } catch (err) {
      if (current !== searchGeneration) return;
      error = err instanceof Error ? err.message : 'Failed to load more relations';
    } finally {
      if (current === searchGeneration) loadingMore = false;
    }
  }

  function handleClearFilters() {
    onFilterChange({});
  }

  function handleRowClick(row: InteractionListRow) {
    selectedInteraction = row;
    detailsOpen = true;
  }

  function openEntityDetails(entity: EntityLike) {
    detailsEntity = entity;
    entityDetailsOpen = true;
  }
</script>

{#snippet searchPanel()}
  <div class="relative">
    {#if isMobile.current}
      <div class="lg:hidden flex items-center gap-2 p-4 border-b">
        {@render resultsActions()}
        <Sheet>
          <SheetTrigger>
            {#snippet child({ props })}
              <Button {...props} variant="outline" class="flex-1">
                <Filter class="h-4 w-4 mr-2" />
                Filters
                {#if activeFilterCount > 0}
                  <Badge variant="secondary" class="ml-2">
                    {activeFilterCount}
                  </Badge>
                {/if}
              </Button>
            {/snippet}
          </SheetTrigger>
          <SheetContent side="left" class="w-[85%] sm:w-[400px] p-0">
            <SheetHeader class="px-6 py-4 border-b">
              <div class="flex items-center justify-between">
                <SheetTitle class="flex items-center gap-2">
                  <Filter class="h-5 w-5 text-primary" />
                  Filters
                </SheetTitle>
                {#if activeFilterCount > 0}
                  <Button
                    variant="ghost"
                    size="sm"
                    onclick={handleClearFilters}
                    class="flex items-center gap-1 text-muted-foreground hover:text-foreground"
                  >
                    <X class="h-4 w-4" />
                    Clear all
                  </Button>
                {/if}
              </div>
            </SheetHeader>
            <div class="h-[calc(100%-4rem)] overflow-y-auto">
              <InteractionFilterSidebar
                filters={effectiveFilters}
                {onFilterChange}
                onClearFilters={handleClearFilters}
                isMobile
                {scopedEntityIds}
                {scopedAnnotationIds}
                {scopeEndpointMode}
                {scopeMode}
                {queryEntityIds}
              />
            </div>
          </SheetContent>
        </Sheet>
      </div>
    {/if}

    {#if topHit}
      <p class="flex flex-wrap items-baseline gap-x-2 border-b px-3 py-2 text-sm">
        <span class="text-muted-foreground">Relations of</span>
        <span class="font-medium">{topHit.label}</span>
        {#if topHit.detail}<span class="text-xs text-muted-foreground">{topHit.detail}</span>{/if}
        <span class="ml-auto font-mono text-xs text-muted-foreground">{topHit.key}</span>
      </p>
    {/if}
    {#if error}
      <div class="p-6">
        <Alert variant="destructive">
          <AlertDescription>{error}</AlertDescription>
        </Alert>
      </div>
    {:else if loading && results.length === 0}
      <!-- the panel's progress bar shows the loading, as for the filters -->
      <p role="status" class="p-4 text-sm text-muted-foreground">Loading relations…</p>
    {:else if results.length > 0}
      <RelationsTable
        rows={results}
        {hasMore}
        {loadingMore}
        onLoadMore={loadMore}
        onRowClick={handleRowClick}
        onEntityClick={openEntityDetails}
      />
    {:else if !loading && results.length === 0}
      <div class="p-6 flex-1 flex items-center justify-center">
        <p class="text-muted-foreground text-center">
          {Object.keys(effectiveFilters).length > 0
            ? 'No relations found matching your criteria.'
            : 'Select filters to explore relations.'}
        </p>
      </div>
    {/if}
  </div>
{/snippet}

{#snippet resultsActions()}
  {#if total !== null}
    <span class="text-xs text-muted-foreground tabular-nums" aria-live="polite"
      >{formatNumber(total)} relations</span
    >
  {/if}
  <ExportMenu
    endpoint="/export"
    body={{ filters: exportFilters }}
    name="relations"
    limit={100_000}
    disabledReason={exportable
      ? undefined
      : 'Choose an entity or filters first: exports cover specific relations, not the whole collection.'}
  />
{/snippet}

{#snippet desktopSidebar()}
  <ExplorerWorkspace
    name="relations"
    panels={['results', 'filters', 'effect', 'interaction_types', 'sources', 'ncbi_tax_id']}
  >
    <WorkspacePanel id="results" title="Results" busy={loading} actions={resultsActions}
      >{@render searchPanel()}</WorkspacePanel
    >
    <InteractionFilterSidebar
      filters={effectiveFilters}
      {onFilterChange}
      onClearFilters={handleClearFilters}
      {scopedEntityIds}
      {scopedAnnotationIds}
      {scopeEndpointMode}
      {scopeMode}
      {queryEntityIds}
    />
  </ExplorerWorkspace>
{/snippet}

<div class="relative flex min-h-0 flex-1 flex-col overflow-hidden">
  {#if isMobile.current}
    {@render searchPanel()}
  {:else}
    {@render desktopSidebar()}
  {/if}

  <!-- Interaction Details Sheet -->
  <InteractionDetailsSheet
    open={detailsOpen}
    onOpenChange={(open) => {
      detailsOpen = open;
    }}
    interaction={selectedInteraction}
  />

  <EntityDetailsDialog bind:open={entityDetailsOpen} entity={detailsEntity} />
</div>
