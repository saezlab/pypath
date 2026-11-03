<script lang="ts">
  import { page } from '$app/stores';
  import { goto } from '$app/navigation';
  import { browser } from '$app/environment';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Card, CardContent } from '$lib/components/ui/card/index.js';
  import ExploreBrowserShell from '$lib/components/explore/ExploreBrowserShell.svelte';
  import AnnotationBrowserTab from '$lib/components/explore/AnnotationBrowserTab.svelte';
  import EntitiesExploreTab from '$lib/components/explore/EntitiesExploreTab.svelte';
  import RelationsExploreTab from '$lib/components/explore/RelationsExploreTab.svelte';
  import SelectionSheet from '$lib/components/selection/SelectionSheet.svelte';
  import { getSelectionStore } from '$lib/stores/selection.svelte';
  import { getSelectionScope, getSelectionScopeSettings } from '$lib/stores/selection-scope.svelte';
  import type { SearchFilters } from '$lib/types/search';

  const selection = getSelectionStore();
  const scopeSettings = getSelectionScopeSettings();
  const rawTab = $derived($page.url.searchParams.get('tab') || 'entities');
  const tab = $derived(
    rawTab === 'relations' || rawTab === 'ontology'
      ? rawTab
      : rawTab === 'summary'
        ? 'ontology'
        : 'entities',
  );
  const query = $derived($page.url.searchParams.get('q') || '');

  let inputRef = $state<HTMLInputElement | null>(null);
  let draftQuery = $state('');
  let filters = $state<SearchFilters>({});

  const scope = $derived(
    getSelectionScope(selection.entityIds, selection.selectedEntityPks, selection.annotationIds, {
      includeAssociatedEntities: scopeSettings.expandSelection,
      includeMembersParticipants: scopeSettings.expandSelection,
      mode: scopeSettings.mode,
      resolveAnnotationEntities: true,
    }),
  );
  const scopeEndpointMode = $derived(scopeSettings.mode === 'intersection' ? 'both' : 'any');
  const entityTabScopePks = $derived(
    Array.from(new Set([...scope.termEntityPks, ...scope.scopedEntityPks])),
  );
  const scopeAnnotationIds = $derived(
    Array.from(new Set([...selection.annotationIds, ...scope.ontologyTermIds])),
  );
  const ontologySelectionScope = $derived({
    entityPks: selection.selectedEntityPks,
    annotationTermIds: selection.annotationIds,
    includeAssociatedEntities: scopeSettings.expandSelection,
    includeMembersParticipants: scopeSettings.expandSelection,
    mode: scopeSettings.mode,
  });
  const hasEntityTabScope = $derived(entityTabScopePks.length > 0);
  const hasRelationScope = $derived(
    scope.scopedEntityPks.length > 0 || scopeAnnotationIds.length > 0,
  );

  $effect(() => {
    draftQuery = query;
  });

  function setTab(next: string) {
    const url = new URL($page.url);
    url.searchParams.set('tab', next);
    goto(url, { replaceState: true, keepFocus: true, noScroll: true });
  }

  function setQuery(next: string) {
    const url = new URL($page.url);
    if (next.trim()) {
      url.searchParams.set('q', next.trim());
    } else {
      url.searchParams.delete('q');
    }
    goto(url, { replaceState: true, keepFocus: true, noScroll: true });
  }

  function submitSearch() {
    setQuery(draftQuery);
  }

  $effect(() => {
    if (!browser) return;
    const onKeyDown = (event: KeyboardEvent) => {
      const target = event.target as HTMLElement | null;
      const tagName = target?.tagName;
      const isTypingTarget =
        tagName === 'INPUT' ||
        tagName === 'TEXTAREA' ||
        tagName === 'SELECT' ||
        target?.isContentEditable;

      if (event.key === '/' && !isTypingTarget) {
        event.preventDefault();
        inputRef?.focus();
        inputRef?.select();
      }
    };

    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  });

  const hasSelection = $derived(
    selection.selectedEntities.length > 0 || selection.selectedAnnotations.length > 0,
  );
  const isSelectionEmpty = $derived(!hasSelection);
  const hasResolvedScopeForTab = $derived(
    tab === 'relations'
      ? hasRelationScope
      : tab === 'ontology'
        ? hasEntityTabScope || scopeAnnotationIds.length > 0
        : hasEntityTabScope,
  );
  const explanationTitle = $derived(
    tab === 'relations'
      ? 'Selected relations'
      : tab === 'ontology'
        ? 'Selected ontology terms'
        : 'Selected entities',
  );
  const explanationText = $derived(
    tab === 'relations'
      ? scopeSettings.mode === 'intersection'
        ? 'Relations connect two entities through an interaction or association. This view shows relations shared by the selected items.'
        : 'Relations connect two entities through an interaction or association. This view shows relations linked to any selected item.'
      : tab === 'ontology'
        ? 'Ontology terms are standardized concepts used to describe and group biological data. This view shows terms linked to your selection.'
        : 'Entities are biological objects such as genes, proteins, chemicals, complexes, and pathways. This view shows entities matching your selection.',
  );
</script>

{#if !(tab === 'relations' && scope.isLoading) && isSelectionEmpty}
  <div class="flex min-h-[60vh] items-center justify-center p-8">
    <Card class="w-full max-w-2xl border-dashed">
      <CardContent class="space-y-4 py-12 text-center">
        <h1 class="text-2xl font-semibold">Selection is empty</h1>
        <p class="text-muted-foreground">Add entities or terms from Explore.</p>
        {#if selection.selectedAnnotations.length > 0}
          <div class="flex flex-wrap justify-center gap-2">
            {#each selection.selectedAnnotations.slice(0, 8) as annotation}
              <Badge variant="secondary">{annotation.label}</Badge>
            {/each}
          </div>
        {/if}
      </CardContent>
    </Card>
  </div>
{:else}
  <ExploreBrowserShell
    fitViewport
    {query}
    {draftQuery}
    onDraftQueryChange={(value) => (draftQuery = value)}
    onSubmitSearch={submitSearch}
    {tab}
    onTabChange={setTab}
    tabs={[
      { value: 'entities', label: 'Entities' },
      { value: 'relations', label: 'Relations' },
      { value: 'ontology', label: 'Ontology' },
    ]}
    searchPlaceholder={tab === 'relations'
      ? 'Search scoped relations…'
      : tab === 'ontology'
        ? 'Search associated ontology terms…'
        : 'Search scoped entities…'}
    {explanationTitle}
    {explanationText}
    bind:searchInputRef={inputRef}
  >
    {#snippet content()}
      {#if scope.isLoading}
        <div class="flex flex-1 items-center justify-center p-8">
          <div class="flex items-center gap-2 text-sm text-muted-foreground">
            <div
              class="h-4 w-4 animate-spin rounded-full border-2 border-primary border-t-transparent"
            ></div>
            <span>Resolving selection scope...</span>
          </div>
        </div>
      {:else if !hasResolvedScopeForTab}
        <div class="flex flex-1 items-center justify-center p-8">
          <Card class="w-full max-w-2xl border-dashed">
            <CardContent class="space-y-3 py-10 text-center">
              <h2 class="text-xl font-semibold">
                {tab === 'relations'
                  ? 'No relation scope'
                  : tab === 'ontology'
                    ? 'No ontology scope'
                    : 'No entities in scope'}
              </h2>
              <p class="text-sm text-muted-foreground">Adjust the scope or add an entity.</p>
            </CardContent>
          </Card>
        </div>
      {:else if tab === 'entities'}
        <EntitiesExploreTab
          {query}
          {filters}
          onFiltersChange={(f) => (filters = f)}
          selectedEntityPks={entityTabScopePks}
        />
      {:else if tab === 'relations'}
        <RelationsExploreTab
          {query}
          {filters}
          onFilterChange={(f) => (filters = f)}
          scopedEntityIds={scope.scopedEntityPks}
          scopedAnnotationIds={scopeAnnotationIds}
          {scopeEndpointMode}
          scopeMode={scopeSettings.mode}
        />
      {:else if tab === 'ontology'}
        <div class="flex min-h-0 flex-1 flex-col overflow-hidden">
          <AnnotationBrowserTab
            {query}
            {filters}
            onFiltersChange={(f) => (filters = f)}
            selectedEntityPks={entityTabScopePks}
            selectedAnnotationIds={scopeAnnotationIds}
            selectionScope={ontologySelectionScope}
          />
        </div>
      {/if}
    {/snippet}

    {#snippet searchActions()}
      {#if hasSelection}<SelectionSheet triggerClass="h-9 shrink-0" />{/if}
    {/snippet}
  </ExploreBrowserShell>
{/if}
