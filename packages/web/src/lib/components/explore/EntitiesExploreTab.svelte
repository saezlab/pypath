<script lang="ts">
  import { getUiPreferences } from '$lib/stores/ui-preferences.svelte';
  import EntityResultsList from '$lib/components/entity/EntityResultsList.svelte';
  import EntityResultCard from '$lib/components/entity/EntityResultCard.svelte';
  import ExplorerWorkspace from '$lib/components/workspace/ExplorerWorkspace.svelte';
  import WorkspacePanel from '$lib/components/workspace/WorkspacePanel.svelte';
  import { untrack } from 'svelte';
  import { page as currentPage } from '$app/state';
  import { formatTaxonomy } from '$lib/utils/taxonomy';
  import { Filter, X } from '@lucide/svelte';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Checkbox } from '$lib/components/ui/checkbox/index.js';
  import { Label } from '$lib/components/ui/label/index.js';
  import {
    Sheet,
    SheetContent,
    SheetHeader,
    SheetTitle,
    SheetTrigger,
  } from '$lib/components/ui/sheet/index.js';
  import EntityGroups from '$lib/components/explore/EntityGroups.svelte';
  import GroupToggle from '$lib/components/explore/GroupToggle.svelte';
  import TopicSourcePicker from './TopicSourcePicker.svelte';
  import EntityDetailsDialog from '$lib/components/entity/EntityDetailsDialog.svelte';
  import {
    fetchEntityExamples,
    fetchEntitiesSearch,
    fetchScopedEntityFacetCounts,
    type EntitySearchCursor,
  } from '$lib/api/client';
  import { getIdentifierTypeLabel } from '$lib/domain/display';
  import { IsMobile } from '$lib/hooks/is-mobile.svelte';
  import type { SearchFilters } from '$lib/types/search';
  import type { EntityWithIdentifiers } from '$lib/types/entities';
  import { getEntityTypeEmoji } from '$lib/utils/entity-types';
  import { formatNumber } from '$lib/utils/format';
  import { splitEmptyOptions } from '$lib/utils/facets';
  import { truncateTitle } from '$lib/actions/truncate-title';

  interface Props {
    query: string;
    filters: SearchFilters;
    onFiltersChange: (filters: SearchFilters) => void;
    selectedEntityIds?: string[];
    selectedEntityPks?: Array<string | number>;
    selectedAnnotationIds?: string[];
  }

  let { query, filters, onFiltersChange, selectedEntityPks, selectedAnnotationIds }: Props =
    $props();

  const isMobile = new IsMobile();
  const RESULTS_PER_PAGE = 20;
  const FACET_PAGE_SIZE = 15;

  interface FilterOption {
    value: string;
    displayName: string;
    icon?: string;
    id?: string | null;
  }

  type EntityResult = EntityWithIdentifiers;
  type ExamplesPage = Awaited<ReturnType<typeof fetchEntityExamples>>;

  const GROUP_EXPLANATION =
    'Group genes and their products by gene reference, and chemicals by connectivity. Source types and molecular evidence stay distinct. Entities without an unambiguous grouping stay separate.';

  let results = $state<EntityResult[]>([]);
  let searchGeneration = 0;
  let loading = $state(true);
  let loadingMore = $state(false);
  let hasMore = $state(true);
  let cursor = $state<EntitySearchCursor | null>(null);
  let facetCountsLoading = $state(true);
  // The option whose click is being applied, until the new counts arrive.
  let pendingOption = $state<string | null>(null);
  let groupsLoading = $state(false);
  let entityTypeOptions = $state<FilterOption[]>([]);
  let sourceOptions = $state<FilterOption[]>([]);
  let taxonomyOptions = $state<FilterOption[]>([]);
  const entityTypeFacetLimit = Infinity;
  const sourceFacetLimit = Infinity;
  let taxonomyFacetLimit = $state(FACET_PAGE_SIZE);
  let scopedFacetCounts = $state<Map<string, number>>(new Map());
  let detailsEntity = $state<EntityResult | null>(null);
  let detailsOpen = $state(false);

  const effectiveFilters = $derived({
    ...filters,
    ...(selectedEntityPks && selectedEntityPks.length > 0
      ? { entity_pks: selectedEntityPks.map(String) }
      : {}),
    ...(selectedAnnotationIds && selectedAnnotationIds.length > 0
      ? { annotation_term_ids: selectedAnnotationIds }
      : {}),
  });

  const facetQueryLimit = $derived(Math.max(taxonomyFacetLimit + 1, 100));

  const activeFilterCount = $derived(
    Object.entries(filters).reduce((count, [, value]) => {
      if (Array.isArray(value)) return count + value.length;
      if (value !== null && value !== undefined) return count + 1;
      return count;
    }, 0),
  );

  // Fetch facet counts (scoped when selection present, global otherwise) and build filter options.
  // Counts reflect the current scope AND query AND all OTHER active filters (cross-facet filtering).
  const facetKey = $derived(
    JSON.stringify([
      currentPage.data.selectedRelease,
      {
        entityIds: selectedEntityPks?.length ? selectedEntityPks : undefined,
        annotationTermIds: selectedAnnotationIds?.length ? selectedAnnotationIds : undefined,
        entityTypes: filters.entity_types,
        sources: filters.sources,
        ncbi_tax_id: filters.taxonomy_ids,
        query: query || undefined,
        facetLimit: facetQueryLimit,
      },
    ]),
  );
  $effect(() => {
    const [, scope] = JSON.parse(facetKey);
    const controller = new AbortController();
    let cancelled = false;
    facetCountsLoading = true;
    untrack(() => fetchScopedEntityFacetCounts(scope, controller.signal))
      .then((counts) => {
        if (cancelled) return;
        const map = new Map<string, number>();
        const types: FilterOption[] = [];
        const sources: FilterOption[] = [];
        const taxonomies: FilterOption[] = [];
        for (const c of counts) {
          map.set(`${c.facetName}:${c.facetValue}`, c.scopedCount);
          if (c.facetName === 'entity_type') {
            types.push(mapEntityTypeOption(c.facetValue));
          } else if (c.facetName === 'source') {
            sources.push({ value: c.facetValue, displayName: c.facetValue });
          } else if (c.facetName === 'taxonomy_id') {
            taxonomies.push(mapTaxonomyOption(c.facetValue, c.facetLabel));
          }
        }
        scopedFacetCounts = map;
        entityTypeOptions = types;
        sourceOptions = sources;
        taxonomyOptions = taxonomies;
      })
      .catch(() => {
        if (!cancelled) {
          scopedFacetCounts = new Map();
          entityTypeOptions = [];
          sourceOptions = [];
          taxonomyOptions = [];
        }
      })
      .finally(() => {
        if (!cancelled) {
          facetCountsLoading = false;
          pendingOption = null;
        }
      });
    return () => {
      cancelled = true;
      controller.abort();
    };
  });

  // Grouping is a remembered preference, shared by the explorer and the selection.
  const preferences = getUiPreferences();
  const groupResults = $derived(preferences.groupResults);
  const showExamples = $derived(
    !query.trim() &&
      !selectedEntityPks &&
      !selectedAnnotationIds &&
      !Object.values(effectiveFilters).some((value) =>
        Array.isArray(value) ? value.length > 0 : value != null && value !== '',
      ),
  );
  const searchKey = $derived(
    JSON.stringify([
      currentPage.data.selectedRelease,
      query,
      effectiveFilters,
      groupResults,
      showExamples,
    ]),
  );
  $effect(() => {
    const [, q, f, grouped, examples] = JSON.parse(searchKey);
    if (grouped && !examples) return; // EntityGroups searches on its own
    const controller = new AbortController();
    const current = ++searchGeneration;

    // Reset internal state
    cursor = null;
    hasMore = true;
    loading = true;
    loadingMore = true;

    untrack(() =>
      examples
        ? fetchEntityExamples(controller.signal)
        : fetchEntitiesSearch(
            { query: q || '', limit: RESULTS_PER_PAGE, cursor: null, filters: f },
            controller.signal,
          ),
    )
      .then((data) => {
        if (current !== searchGeneration) return;
        const groups = grouped && examples ? (data as ExamplesPage).groups : null;
        results = (groups
          ? groups.map((group) => group.entity)
          : data.entities) as unknown as EntityResult[];
        cursor = data.nextCursor;
        hasMore = data.nextCursor !== null;
      })
      .catch(() => {
        if (current !== searchGeneration) return;
        results = [];
        hasMore = false;
      })
      .finally(() => {
        if (current !== searchGeneration) return;
        loading = false;
        loadingMore = false;
      });
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
      const data = await fetchEntitiesSearch({
        query: query || '',
        limit: RESULTS_PER_PAGE,
        cursor,
        filters: effectiveFilters,
      });
      if (current !== searchGeneration) return;
      results = [...results, ...(data.entities as unknown as EntityResult[])];
      cursor = data.nextCursor;
      hasMore = data.nextCursor !== null;
    } finally {
      if (current === searchGeneration) loadingMore = false;
    }
  }

  function handleClearFilters() {
    onFiltersChange({});
  }

  function handleFilterChange(next: {
    entity_types?: string[];
    sources?: string[];
    taxonomy_ids?: string[];
  }) {
    onFiltersChange({
      ...filters,
      ...next,
    });
  }

  function handleFilterToggle(
    filterKey: 'entity_types' | 'sources' | 'taxonomy_ids',
    value: string,
  ) {
    pendingOption = `${filterKey}:${value}`;
    const currentValues = filters[filterKey] || [];
    const nextValues = currentValues.includes(value)
      ? currentValues.filter((entry) => entry !== value)
      : [...currentValues, value];

    handleFilterChange({
      [filterKey]: nextValues.length > 0 ? nextValues : undefined,
    });
  }

  function mapTaxonomyOption(value: string, name?: string | null): FilterOption {
    return {
      value,
      displayName: formatTaxonomy(value, name) || value,
      id: value,
    };
  }

  function facetNameForFilter(filterKey: 'entity_types' | 'sources' | 'taxonomy_ids'): string {
    if (filterKey === 'entity_types') return 'entity_type';
    if (filterKey === 'taxonomy_ids') return 'taxonomy_id';
    return 'source';
  }

  function mapEntityTypeOption(value: string): FilterOption {
    // Handle raw format from DB: "MI:0326:Protein"
    const rawMatch = value.match(/^([A-Z]+):(\d+):(.+)$/i);
    if (rawMatch) {
      const displayName = rawMatch[3];
      return {
        value,
        displayName,
        icon: getEntityTypeEmoji(displayName),
        id: `${rawMatch[1]}:${rawMatch[2]}`,
      };
    }

    // Handle old formatted value: "protein:MI:0326"
    const formattedMatch = value.match(/^(.+):([A-Z]+:\d+)$/i);
    if (formattedMatch) {
      const displayName = formattedMatch[1];
      return {
        value,
        displayName,
        icon: getEntityTypeEmoji(displayName),
        id: formattedMatch[2],
      };
    }

    const displayName = getIdentifierTypeLabel(value) || value;
    return {
      value,
      displayName,
      icon: getEntityTypeEmoji(displayName),
      id: value,
    };
  }

  function openDetails(entity: EntityResult) {
    detailsEntity = entity;
    detailsOpen = true;
  }
</script>

{#if isMobile.current}
  <div class="flex h-full min-h-0 flex-col overflow-hidden">
    <div class="flex items-center gap-2 border-b p-4">
      {@render groupToggle()}
      <Sheet>
        <SheetTrigger>
          {#snippet child({ props })}
            <Button {...props} variant="outline" class="flex-1">
              <Filter class="mr-2 size-4" />
              Filters
              {#if activeFilterCount > 0}
                <Badge variant="secondary" class="ml-2">{activeFilterCount}</Badge>
              {/if}
            </Button>
          {/snippet}
        </SheetTrigger>
        <SheetContent side="left" class="w-[85%] sm:w-[400px] overflow-y-auto p-0">
          <SheetHeader class="border-b px-6 py-4">
            <div class="flex items-center justify-between">
              <SheetTitle class="flex items-center gap-2">
                <Filter class="size-5 text-primary" />
                Filters
              </SheetTitle>
              {#if activeFilterCount > 0}
                <Button
                  variant="ghost"
                  size="sm"
                  onclick={handleClearFilters}
                  class="flex items-center gap-1 text-muted-foreground"
                >
                  <X class="size-4" />
                  Clear all
                </Button>
              {/if}
            </div>
          </SheetHeader>
          <div class="p-4">{@render filterSnippet(true)}</div>
        </SheetContent>
      </Sheet>
    </div>
    <div class="min-h-0 flex-1">
      {@render resultsPane()}
    </div>
  </div>
{:else}
  <ExplorerWorkspace name="entities" panels={['results', 'entity_types', 'sources', 'ncbi_tax_id']}>
    <WorkspacePanel
      id="results"
      title="Results"
      bodyClass=""
      actions={groupToggle}
      busy={groupResults && !showExamples ? groupsLoading : loading}
      >{@render resultsPane()}</WorkspacePanel
    >
    <WorkspacePanel id="entity_types" title="Entity types" busy={facetCountsLoading}
      >{@render filterSection(
        'Entity Types',
        'entity_types',
        entityTypeOptions,
        filters.entity_types || [],
        entityTypeFacetLimit,
        () => {},
        handleClearFilters,
      )}</WorkspacePanel
    >
    <WorkspacePanel id="sources" title="Sources" busy={facetCountsLoading}
      >{@render topicCards()}{@render filterSection(
        'Data Sources',
        'sources',
        sourceOptions,
        filters.sources || [],
        sourceFacetLimit,
        () => {},
      )}</WorkspacePanel
    >
    <WorkspacePanel id="ncbi_tax_id" title="Taxonomy" busy={facetCountsLoading}
      >{@render filterSection(
        'Taxonomy',
        'taxonomy_ids',
        taxonomyOptions,
        filters.taxonomy_ids || [],
        taxonomyFacetLimit,
        () => (taxonomyFacetLimit += FACET_PAGE_SIZE),
      )}</WorkspacePanel
    >
  </ExplorerWorkspace>
{/if}

<EntityDetailsDialog bind:open={detailsOpen} entity={detailsEntity} />

{#snippet filterSnippet(mobile = false)}
  <div
    class={mobile ? 'space-y-6' : 'space-y-6 px-3 pt-8 pb-4'}
    class:opacity-70={facetCountsLoading}
  >
    {@render filterSection(
      'Entity Types',
      'entity_types',
      entityTypeOptions,
      filters.entity_types || [],
      entityTypeFacetLimit,
      () => {},
      mobile ? undefined : handleClearFilters,
    )}
    {@render topicCards()}
    {@render filterSection(
      'Data Sources',
      'sources',
      sourceOptions,
      filters.sources || [],
      sourceFacetLimit,
      () => {},
    )}
    {@render filterSection(
      'Taxonomy',
      'taxonomy_ids',
      taxonomyOptions,
      filters.taxonomy_ids || [],
      taxonomyFacetLimit,
      () => (taxonomyFacetLimit += FACET_PAGE_SIZE),
    )}
  </div>
{/snippet}

{#snippet filterSection(
  title: string,
  filterKey: 'entity_types' | 'sources' | 'taxonomy_ids',
  options: FilterOption[],
  selectedValues: string[],
  visibleLimit: number,
  onLoadMore: () => void,
  onClear?: () => void,
)}
  {#if options.length > 0}
    {@const split = splitEmptyOptions(
      options,
      (option) => scopedFacetCounts.get(`${facetNameForFilter(filterKey)}:${option.value}`),
      (option) => selectedValues.includes(option.value),
    )}
    {@const listed = split.shown}
    <div>
      <div class="mb-3 flex items-center justify-between gap-2">
        <h4 class="text-[11px] font-semibold uppercase tracking-[0.16em] text-muted-foreground">
          {title}
        </h4>
        {#if onClear && activeFilterCount > 0}
          <Button
            variant="ghost"
            size="xs"
            onclick={onClear}
            class="h-auto px-0 text-[11px] text-muted-foreground"
          >
            Clear all
          </Button>
        {/if}
      </div>
      <div class="max-h-64 space-y-1 overflow-y-auto pr-2">
        {#each listed.slice(0, visibleLimit) as option (option.value)}
          {@const selected = selectedValues.includes(option.value)}
          {@const count = scopedFacetCounts.get(`${facetNameForFilter(filterKey)}:${option.value}`)}
          <div class="flex items-center justify-between gap-2 py-0.5">
            <Label
              for={`${filterKey}-${option.value}`}
              class="flex min-w-0 flex-1 cursor-pointer items-center gap-2 text-sm font-normal leading-5 text-foreground {selected
                ? 'font-medium'
                : ''}"
            >
              <Checkbox
                id={`${filterKey}-${option.value}`}
                checked={selected}
                onCheckedChange={() => handleFilterToggle(filterKey, option.value)}
                class={selected ? 'h-4 w-4 flex-shrink-0 border-primary' : 'h-4 w-4 flex-shrink-0'}
              />
              <span class="truncate" use:truncateTitle={option.displayName}>
                {#if option.icon}<span class="mr-1.5">{option.icon}</span>{/if}
                {option.displayName}
              </span>
            </Label>
            {#if facetCountsLoading && pendingOption === `${filterKey}:${option.value}`}
              <span class="pending-spinner text-muted-foreground" aria-hidden="true"></span>
            {:else if count != null}
              <span class="text-xs text-muted-foreground tabular-nums flex-shrink-0">
                {formatNumber(count)}
              </span>
            {/if}
          </div>
        {/each}
        {#if listed.length > visibleLimit}
          <Button variant="ghost" size="sm" class="mt-2 w-full text-xs" onclick={onLoadMore}>
            Load more
          </Button>
        {/if}
      </div>
    </div>
  {:else if facetCountsLoading}
    <div class="space-y-2">
      <h4 class="mb-3 text-[11px] font-semibold uppercase tracking-[0.16em] text-muted-foreground">
        {title}
      </h4>
      <p class="text-sm text-muted-foreground">Loading filters…</p>
    </div>
  {/if}
{/snippet}

<!-- Topic cards that select the sources of a topic. -->
{#snippet topicCards()}
  <TopicSourcePicker
    selected={filters.sources || []}
    available={sourceOptions.map((option) => option.value)}
    onChange={(sources) => handleFilterChange({ sources: sources.length ? sources : undefined })}
  />
{/snippet}

{#snippet groupToggle()}
  <GroupToggle
    bind:pressed={() => preferences.groupResults, (pressed) => preferences.setGroupResults(pressed)}
    explanation={GROUP_EXPLANATION}
    class={isMobile.current ? 'h-9' : 'h-5'}
  />
{/snippet}

{#snippet resultsPane()}
  <div>
    {#if groupResults && !showExamples}
      <EntityGroups
        {query}
        filters={effectiveFilters}
        renderMember={resultCard}
        bind:loading={groupsLoading}
      />
    {:else if loading && results.length === 0}
      <!-- the panel's progress bar shows the loading, as for the filters -->
      <p role="status" class="p-4 text-sm text-muted-foreground">Loading results…</p>
    {:else if results.length > 0}
      <EntityResultsList>
        {#each results as result}
          {@render resultCard(result)}
        {/each}
      </EntityResultsList>
      {#if hasMore && !showExamples}
        <div class="flex justify-center py-8">
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
    {:else}
      <div class="flex flex-col items-center justify-center gap-2 px-4 py-16 text-center">
        <div class="text-lg font-semibold">
          {showExamples ? 'Search this collection' : 'No entities found'}
        </div>
        <p class="max-w-2xl text-sm text-muted-foreground">
          Try a gene symbol, UniProt identifier, chemical name, or broader text query.
        </p>
      </div>
    {/if}
  </div>
{/snippet}

{#snippet resultCard(result: (typeof results)[0])}
  <EntityResultCard {result} onopen={openDetails} />
{/snippet}
