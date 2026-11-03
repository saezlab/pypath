<script lang="ts">
  // The filters in use, one chip per category with the number of chosen values: a chip
  // clears its category, "Clear all" clears every filter. The values themselves are in the
  // filter panels.
  import { X } from '@lucide/svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import type { SearchFilters } from '$lib/types/search';

  type Key = keyof SearchFilters;
  let {
    filters,
    onClear,
  }: {
    filters: SearchFilters;
    /** Clears the filter keys of one category, or of all categories. */
    onClear: (keys: Key[]) => void;
  } = $props();

  // The categories as the filter panels name them, with the filter keys of each.
  const CATEGORIES: { label: string; keys: Key[] }[] = [
    { label: 'Sources', keys: ['sources'] },
    { label: 'Entity types', keys: ['entity_types'] },
    { label: 'Relation types', keys: ['predicates', 'relation_categories'] },
    {
      label: 'Effect',
      keys: [
        'object_aspect_qualifier',
        'object_direction_qualifier',
        'causal_mechanism_qualifier',
        'signs',
      ],
    },
    { label: 'Participant types', keys: ['interaction_types', 'partner_types'] },
    { label: 'Taxonomy', keys: ['taxonomy_ids'] },
    { label: 'Annotations', keys: ['annotation_term_ids'] },
  ];

  function chosen(key: Key): number {
    const value = filters[key];
    if (Array.isArray(value)) return value.length;
    return value === undefined || value === null || value === '' ? 0 : 1;
  }
  const applied = $derived(
    CATEGORIES.map((category) => ({
      ...category,
      count: category.keys.reduce((sum, key) => sum + chosen(key), 0),
    })).filter((category) => category.count > 0),
  );
</script>

{#if applied.length}
  <div
    class="mb-3 flex min-w-0 flex-wrap items-center gap-1.5 text-sm"
    role="group"
    aria-label="Applied filters"
  >
    {#each applied as category (category.label)}
      <button
        type="button"
        class="inline-flex items-center gap-1.5 rounded-md border bg-background py-0.5 pr-1 pl-2 text-xs hover:border-ring"
        aria-label={`Clear ${category.label} (${category.count})`}
        title={`Clear ${category.label}`}
        onclick={() => onClear(category.keys)}
        >{category.label}<span class="rounded bg-muted px-1.5 text-[11px] font-medium tabular-nums"
          >{category.count}</span
        ><X class="size-3 text-muted-foreground" /></button
      >
    {/each}
    <Button
      variant="ghost"
      size="sm"
      class="h-6 px-2 text-xs text-muted-foreground"
      onclick={() => onClear(applied.flatMap((category) => category.keys))}>Clear all</Button
    >
  </div>
{/if}
