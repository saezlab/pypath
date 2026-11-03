<script lang="ts">
  import { goto } from '$app/navigation';
  import { page } from '$app/stores';
  import { browser } from '$app/environment';
  import ExploreBrowserShell from '$lib/components/explore/ExploreBrowserShell.svelte';
  import EntitiesExploreTab from '$lib/components/explore/EntitiesExploreTab.svelte';
  import RelationsExploreTab from '$lib/components/explore/RelationsExploreTab.svelte';
  import AppliedFilters from '$lib/components/explore/AppliedFilters.svelte';
  import SelectionSheet from '$lib/components/selection/SelectionSheet.svelte';
  import { getSelectionStore } from '$lib/stores/selection.svelte';
  import type { SearchFilters } from '$lib/types/search';

  const selection = getSelectionStore();

  let inputRef = $state<HTMLInputElement | null>(null);
  let draftQuery = $state('');
  let entityFilters = $state<SearchFilters>({});
  let interactionFilters = $state<SearchFilters>({});
  let selectionSheetOpen = $state(false);

  const rawTab = $derived($page.url.searchParams.get('tab') || 'entity');
  const tab = $derived(rawTab === 'relations' ? 'relations' : 'entity');
  const query = $derived($page.url.searchParams.get('q') || '');
  const sourceKey = $derived(JSON.stringify($page.url.searchParams.getAll('source')));
  const sourceScope: string[] = $derived(JSON.parse(sourceKey));
  const entityMatchFilters = $derived({ ...entityFilters, sources: sourceScope });
  const interactionMatchFilters = $derived({ ...interactionFilters, sources: sourceScope });
  function updateFilters(next: SearchFilters, relation = false) {
    if (relation) interactionFilters = next;
    else entityFilters = next;
    const sources = next.sources || [];
    entityFilters = { ...entityFilters, sources };
    interactionFilters = { ...interactionFilters, sources };
    const url = new URL($page.url);
    url.searchParams.delete('source');
    for (const source of sources) url.searchParams.append('source', source);
    if (url.href !== $page.url.href)
      goto(url, { replaceState: true, keepFocus: true, noScroll: true });
  }

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
    setQuery(inputRef?.value ?? draftQuery);
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
        return;
      }

      if (
        (event.key === 's' || event.key === 'S') &&
        !isTypingTarget &&
        selection.totalSelectionCount > 0
      ) {
        event.preventDefault();
        selectionSheetOpen = true;
      }
    };

    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  });

  const searchPlaceholder = $derived(
    tab === 'relations' ? 'Search relations…' : 'Search entities…',
  );
  const explanationTitle = $derived(tab === 'relations' ? 'Search relations' : 'Search entities');
  const explanationText = $derived(
    tab === 'relations'
      ? 'Relations connect two entities through an interaction or association. Search and filter them by participant, type, or source.'
      : 'Entities are biological objects such as genes, proteins, chemicals, complexes, and pathways. Search and filter the full index.',
  );
</script>

<svelte:window />

<ExploreBrowserShell
  fitViewport
  {query}
  {draftQuery}
  onDraftQueryChange={(value) => (draftQuery = value)}
  onSubmitSearch={submitSearch}
  {tab}
  onTabChange={setTab}
  tabs={[
    { value: 'entity', label: 'Entities' },
    { value: 'relations', label: 'relations' },
  ]}
  {searchPlaceholder}
  {explanationTitle}
  {explanationText}
  bind:searchInputRef={inputRef}
>
  {#snippet content()}
    <AppliedFilters
      filters={tab === 'entity' ? entityMatchFilters : interactionMatchFilters}
      onClear={(keys) => {
        const next: SearchFilters = {
          ...(tab === 'entity' ? entityMatchFilters : interactionMatchFilters),
        };
        for (const key of keys) delete next[key];
        updateFilters(next, tab === 'relations');
      }}
    />
    {#if tab === 'entity'}
      <EntitiesExploreTab
        {query}
        filters={entityMatchFilters}
        onFiltersChange={(f) => updateFilters(f)}
      />
    {:else}
      <RelationsExploreTab
        {query}
        entitySearchFilters={entityMatchFilters}
        filters={interactionMatchFilters}
        onFilterChange={(f) => updateFilters(f, true)}
      />
    {/if}
  {/snippet}

  {#snippet footerCta()}
    <!-- Opened with the "s" shortcut; the header's Selection link shows the count. -->
    <SelectionSheet bind:open={selectionSheetOpen} showTrigger={false} />
  {/snippet}
</ExploreBrowserShell>
