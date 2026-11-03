<script lang="ts">
  import { goto } from '$app/navigation';
  import { page } from '$app/state';
  import { ArrowDown, ArrowUp, Database, Filter, Search } from '@lucide/svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Checkbox } from '$lib/components/ui/checkbox/index.js';
  import { Label } from '$lib/components/ui/label/index.js';
  import * as Table from '$lib/components/ui/table/index.js';
  import ResourceRow from '$lib/components/resources/ResourceRow.svelte';
  import ResourceDetailsDialog from '$lib/components/resources/ResourceDetailsDialog.svelte';
  // Topic cards above the list; they replace the topic group of the filters.
  import TopicBlocks from '$lib/components/resources/TopicBlocks.svelte';
  import TopicDownloadDialog from '$lib/components/resources/TopicDownloadDialog.svelte';
  import type { ResourceTag } from '$lib/resources/types';
  import type { ResourceRecord } from '$lib/resources/types';
  import {
    defaultDirection,
    dimensions,
    resourceTags,
    filterResources,
    sortResources,
    formatBytes,
    type SortDirection,
  } from '$lib/resources/catalog';
  import type { PageData } from './$types';

  let { data }: { data: PageData } = $props();
  let inputRef = $state<HTMLInputElement | null>(null);
  let draftQuery = $state('');
  let showFilters = $state(false);
  let detailsResource = $state<ResourceRecord | null>(null);
  let detailsOpen = $state(false);
  let downloadTopics = $state<ResourceTag[]>([]);
  let downloadOpen = $state(false);
  const query = $derived(page.url.searchParams.get('q') || '');
  const selectedTags = $derived(page.url.searchParams.getAll('tag'));
  const licenses = $derived(page.url.searchParams.getAll('use'));
  const sort = $derived(page.url.searchParams.get('sort') || 'name');
  const order = $derived(
    (page.url.searchParams.get('order') as SortDirection | null) || defaultDirection(sort),
  );
  // Sample builds cap every dataset at the same size, so the note is shown once for the page.
  const sampleLimit = $derived(
    data.resources.find((resource) => resource.sample_build && resource.max_records)?.max_records,
  );
  const tags = $derived(resourceTags(data.resources));
  const filtered = $derived(
    sortResources(filterResources(data.resources, query, selectedTags, licenses), sort, order),
  );
  const totals = $derived(
    filtered.reduce(
      (sum, r) => ({
        entities: sum.entities + r.entity_count,
        relations: sum.relations + r.interaction_count,
        bytes: sum.bytes + r.total_size_bytes,
      }),
      { entities: 0, relations: 0, bytes: 0 },
    ),
  );
  const activeCount = $derived(selectedTags.length + licenses.length);
  const number = new Intl.NumberFormat('en-US', { notation: 'compact', maximumFractionDigits: 1 });
  const columns = [
    { key: 'name', label: 'Resource', align: 'left' },
    { key: 'entities', label: 'Entities', align: 'right' },
    { key: 'relations', label: 'Relations', align: 'right' },
    { key: 'size', label: 'Size', align: 'right' },
  ] as const;
  $effect(() => {
    draftQuery = query;
  });

  function tagCount(id: string, dimension: string) {
    const otherGroups = selectedTags.filter(
      (selected) => tags.find((tag) => tag.id === selected)?.dimension !== dimension,
    );
    return filterResources(data.resources, query, [...otherGroups, id], licenses).length;
  }
  function licenseCount(audience: string) {
    return filterResources(data.resources, query, selectedTags, [
      ...new Set([...licenses, audience]),
    ]).length;
  }
  function navigate(url: URL) {
    url.searchParams.delete('tab');
    void goto(url, { keepFocus: true, noScroll: true });
  }
  function setParam(key: string, value: string) {
    const url = new URL(page.url);
    if (value.trim()) url.searchParams.set(key, value.trim());
    else url.searchParams.delete(key);
    navigate(url);
  }
  function toggle(key: string, id: string) {
    const url = new URL(page.url);
    const values = url.searchParams.getAll(key);
    url.searchParams.delete(key);
    for (const value of values.filter((value) => value !== id)) url.searchParams.append(key, value);
    if (!values.includes(id)) url.searchParams.append(key, id);
    navigate(url);
  }
  function setSort(key: string) {
    const url = new URL(page.url);
    const next: SortDirection =
      key === sort ? (order === 'asc' ? 'desc' : 'asc') : defaultDirection(key);
    url.searchParams.set('sort', key);
    // The default direction stays implicit so shared links stay short.
    if (next === defaultDirection(key)) url.searchParams.delete('order');
    else url.searchParams.set('order', next);
    navigate(url);
  }
  function clear() {
    const url = new URL(page.url);
    for (const key of ['q', 'tag', 'use']) url.searchParams.delete(key);
    navigate(url);
  }
  function keydown(event: KeyboardEvent) {
    const target = event.target as HTMLElement;
    if (
      event.key === '/' &&
      !event.metaKey &&
      !event.ctrlKey &&
      !event.altKey &&
      !target.closest('input, textarea, select, [contenteditable="true"]')
    ) {
      event.preventDefault();
      inputRef?.focus();
      inputRef?.select();
    }
  }
</script>

<svelte:head><title>Resources · OmniPath</title></svelte:head>
<svelte:window onkeydown={keydown} />

<h1 class="sr-only">Resources</h1>
<!-- Fits the viewport on large screens (the table and the filters scroll on their own); a
     single scrolling page on small ones. -->
<div
  class="flex h-full min-h-0 flex-col gap-3 overflow-y-auto pt-2 pb-6 lg:overflow-hidden lg:pb-0"
>
  <div class="flex items-center gap-2">
    <form
      class="flex min-w-0 flex-1 items-center gap-2 rounded-md border border-border/60 bg-background px-3 focus-within:border-ring"
      onsubmit={(event) => {
        event.preventDefault();
        setParam('q', draftQuery);
      }}
    >
      <Search class="size-4 shrink-0 text-muted-foreground" />
      <input
        bind:this={inputRef}
        bind:value={draftQuery}
        type="search"
        aria-label="Search resources"
        placeholder="Search resources, topics, molecules…"
        class="h-9 min-w-0 flex-1 bg-transparent text-sm text-foreground placeholder:text-muted-foreground focus:outline-none"
      />
      <button
        type="submit"
        aria-label="Search"
        title="Search (press / to focus)"
        class="inline-flex h-6 shrink-0 items-center rounded border border-border/70 px-1.5 font-mono text-[11px] text-muted-foreground transition-colors hover:border-ring hover:text-foreground"
        >↵</button
      >
    </form>
    <Button
      variant="outline"
      class="lg:hidden"
      aria-expanded={showFilters}
      aria-controls="resource-filters"
      onclick={() => (showFilters = !showFilters)}
    >
      <Filter class="size-4" />
      Filters
      {#if activeCount}<Badge variant="secondary">{activeCount}</Badge>{/if}
    </Button>
  </div>

  <TopicBlocks
    resources={data.resources}
    {selectedTags}
    onToggle={(id) => toggle('tag', id)}
    onClear={(ids) => {
      const url = new URL(page.url);
      const keep = url.searchParams.getAll('tag').filter((tag) => !ids.includes(tag));
      url.searchParams.delete('tag');
      for (const tag of keep) url.searchParams.append('tag', tag);
      navigate(url);
    }}
    onDownload={(topics) => {
      downloadTopics = topics;
      downloadOpen = true;
    }}
  />

  {#if sampleLimit}
    <p class="rounded-md border px-3 py-2 text-xs text-muted-foreground" role="note">
      This release is a sample build: up to {sampleLimit.toLocaleString('en')} input records per dataset.
      Supporting records can make the counts larger.
    </p>
  {/if}

  <div
    class="grid items-start gap-3 lg:min-h-0 lg:flex-1 lg:grid-cols-[minmax(0,1fr)_260px] lg:grid-rows-[minmax(0,1fr)] lg:items-stretch"
  >
    <section
      class="min-w-0 overflow-hidden rounded-xl border bg-background lg:flex lg:min-h-0 lg:flex-col"
    >
      <header
        class="flex shrink-0 items-center justify-between gap-3 border-b bg-muted/20 px-3 py-2 text-sm font-medium"
      >
        <span>Resources</span>
        <span class="flex items-center gap-3 text-xs font-normal text-muted-foreground">
          <span class="hidden items-center gap-1.5 sm:inline-flex"
            ><span class="legend-dot entities"></span>Entities<span class="legend-dot relations"
            ></span>Relations</span
          >
          <span class="tabular-nums" aria-live="polite"
            >{filtered.length} of {data.resources.length}</span
          >
        </span>
      </header>
      {#if filtered.length}
        <div class="resource-table lg:min-h-0 lg:flex-1">
          <Table.Root class="table-fixed">
            <Table.Header class="sticky top-0 z-10 bg-background">
              <Table.Row class="hover:bg-transparent">
                {#each columns as column (column.key)}
                  <Table.Head
                    class={column.key === 'name'
                      ? 'pl-3'
                      : column.key === 'entities'
                        ? 'w-24 text-right'
                        : 'hidden w-24 text-right sm:table-cell'}
                    aria-sort={sort === column.key
                      ? order === 'asc'
                        ? 'ascending'
                        : 'descending'
                      : undefined}
                  >
                    <button
                      type="button"
                      class={`inline-flex items-center gap-1 transition-colors hover:text-foreground ${
                        sort === column.key ? 'font-semibold text-foreground' : ''
                      }`}
                      onclick={() => setSort(column.key)}
                    >
                      {column.label}
                      {#if sort === column.key}{#if order === 'asc'}<ArrowUp
                            class="size-3"
                          />{:else}<ArrowDown class="size-3" />{/if}{/if}
                    </button>
                  </Table.Head>
                {/each}
                <Table.Head class="hidden w-44 xl:table-cell">Topic</Table.Head>
                <Table.Head class="w-14"><span class="sr-only">Download</span></Table.Head>
              </Table.Row>
            </Table.Header>
            <Table.Body>
              {#each filtered as resource (resource.resource_id)}<ResourceRow
                  {resource}
                  {selectedTags}
                  onTagClick={(id) => toggle('tag', id)}
                  onOpen={(next) => {
                    detailsResource = next;
                    detailsOpen = true;
                  }}
                />{/each}
            </Table.Body>
          </Table.Root>
        </div>
        <footer
          class="flex shrink-0 flex-wrap gap-x-4 gap-y-1 border-t px-3 py-2 text-xs text-muted-foreground tabular-nums"
        >
          <span class="inline-flex items-center gap-1.5 font-medium text-foreground"
            ><Database class="size-3.5" />{filtered.length} resources</span
          >
          <span>{number.format(totals.entities)} entities</span>
          <span>{number.format(totals.relations)} relations</span>
          <span>{formatBytes(totals.bytes)} total size</span>
        </footer>
      {:else}
        <div class="flex flex-col items-center gap-2 px-4 py-16 text-center">
          <Database class="size-7 text-muted-foreground" />
          <h2 class="text-base font-semibold">
            {data.resourcesUnavailable ? 'Resources are unavailable' : 'No matching resources'}
          </h2>
          <p class="text-sm text-muted-foreground">
            {data.resourcesUnavailable
              ? 'The API has no resources available for this release.'
              : 'Try another search or remove a filter to broaden your results.'}
          </p>
          {#if !data.resourcesUnavailable}<Button variant="outline" size="sm" onclick={clear}
              >Clear all filters</Button
            >{/if}
        </div>
      {/if}
    </section>

    <aside
      id="resource-filters"
      aria-label="Resource filters"
      class={`order-first overflow-hidden rounded-xl border bg-background lg:order-none lg:flex lg:min-h-0 lg:flex-col ${
        showFilters ? '' : 'hidden'
      }`}
    >
      <header
        class="flex shrink-0 items-center justify-between gap-2 border-b bg-muted/20 px-3 py-2 text-sm font-medium"
      >
        <span>Filters</span>
        {#if activeCount || query}<Button
            variant="ghost"
            size="xs"
            class="h-auto px-1 text-xs text-muted-foreground"
            onclick={clear}>Clear all</Button
          >{/if}
      </header>
      <div class="lg:min-h-0 lg:flex-1 lg:overflow-y-auto">
        <!-- The topics are the cards above the list. -->
        {#each dimensions.filter((dimension) => dimension.id !== 'topic') as dimension (dimension.id)}
          {@const options = tags
            .filter((tag) => tag.dimension === dimension.id)
            .map((tag) => ({ tag, count: tagCount(tag.id, tag.dimension) }))
            .filter(({ tag, count }) => count > 0 || selectedTags.includes(tag.id))
            .sort((a, b) => b.count - a.count || a.tag.label.localeCompare(b.tag.label))}
          {#if options.length}
            <div
              role="group"
              aria-label={dimension.label}
              class="space-y-1 border-b px-3 py-2.5 last:border-b-0"
            >
              <h3
                class="mb-1 text-[11px] font-semibold tracking-wide text-muted-foreground uppercase"
              >
                {dimension.label}
              </h3>
              {#each options as { tag, count } (tag.id)}
                {@render option(
                  `tag-${tag.id}`,
                  tag.label,
                  count,
                  selectedTags.includes(tag.id),
                  () => toggle('tag', tag.id),
                )}
              {/each}
            </div>
          {/if}
        {/each}
        <div role="group" aria-label="License use" class="space-y-1 px-3 py-2.5">
          <h3 class="mb-1 text-[11px] font-semibold tracking-wide text-muted-foreground uppercase">
            License use
          </h3>
          {#each ['academic', 'commercial'] as audience}
            {@render option(
              `use-${audience}`,
              audience === 'academic' ? 'Academic use' : 'Commercial use',
              licenseCount(audience),
              licenses.includes(audience),
              () => toggle('use', audience),
            )}
          {/each}
        </div>
      </div>
    </aside>
  </div>
</div>

<TopicDownloadDialog bind:open={downloadOpen} topics={downloadTopics} resources={data.resources} />

<ResourceDetailsDialog
  bind:open={detailsOpen}
  resource={detailsResource}
  {selectedTags}
  onTagClick={(id) => toggle('tag', id)}
/>

{#snippet option(id: string, label: string, count: number, checked: boolean, onToggle: () => void)}
  <div class="flex items-center justify-between gap-2">
    <Label
      for={id}
      class="flex min-w-0 flex-1 cursor-pointer items-center gap-2 text-sm leading-5 font-normal {checked
        ? 'font-medium'
        : ''}"
    >
      <Checkbox {id} {checked} onCheckedChange={onToggle} class="size-4 shrink-0" />
      <span class="truncate">{label}</span>
    </Label>
    <span class="shrink-0 text-xs text-muted-foreground tabular-nums">{count}</span>
  </div>
{/snippet}

<style>
  /* Entity and relation colours, shared by the legend and each row's ring. */
  :global(:root) {
    --resource-entities: #0891b2;
    --resource-relations: #16a34a;
  }
  :global(.dark) {
    --resource-entities: #22d3ee;
    --resource-relations: #4ade80;
  }
  /* On large screens the table scrolls inside its panel, under a sticky header. */
  @media (min-width: 1024px) {
    .resource-table :global([data-slot='table-container']) {
      height: 100%;
      overflow: auto;
    }
  }
  .legend-dot {
    display: inline-block;
    width: 8px;
    height: 8px;
    border-radius: 2px;
  }
  .legend-dot.entities {
    background: var(--resource-entities);
  }
  .legend-dot.relations {
    margin-left: 6px;
    background: var(--resource-relations);
  }
</style>
