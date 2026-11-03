<script lang="ts">
  import { Minus } from '@lucide/svelte';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import EntityBadge from '$lib/components/entity/EntityBadge.svelte';
  import RelationLabel from './RelationLabel.svelte';
  import {
    getEntityBadgeIdentifier,
    getEntityDisplayName,
    getEntityTypeLabel,
    type EntityLike,
  } from '$lib/domain/display';
  import type { InteractionListRow } from '$lib/types/interactions';
  import { formatNumber } from '$lib/utils/format';

  interface Props {
    rows: InteractionListRow[];
    hasMore?: boolean;
    loadingMore?: boolean;
    onLoadMore?: () => void;
    onRowClick: (row: InteractionListRow) => void;
    onEntityClick: (entity: EntityLike) => void;
  }

  let {
    rows,
    hasMore = false,
    loadingMore = false,
    onLoadMore,
    onRowClick,
    onEntityClick,
  }: Props = $props();

  // The layout follows the width of the panel, not of the window (panels are resized and
  // moved): from 32rem a table of four columns, the relation column never narrower than
  // 7.5rem; below it each relation is two lines, "From → To", then the relation and its
  // evidence across the whole width, so that nothing overlaps.
  const NARROW =
    "grid-cols-[minmax(0,1fr)_auto_minmax(0,1fr)] [grid-template-areas:'from_arrow_to'_'rel_rel_ev']";
  const WIDE =
    "@lg:grid-cols-[minmax(0,1fr)_minmax(7.5rem,0.75fr)_minmax(0,1fr)_3.5rem] @lg:[grid-template-areas:'from_rel_to_ev']";
</script>

{#snippet endpoint(entity: EntityLike)}
  <button
    type="button"
    onclick={(event) => {
      event.stopPropagation();
      onEntityClick(entity);
    }}
    class="block w-full cursor-pointer rounded-lg p-0.5 text-left transition-all hover:bg-primary/10 hover:ring-2 hover:ring-primary/35 focus:outline-none focus-visible:ring-2 focus-visible:ring-ring"
  >
    <EntityBadge
      displayName={getEntityDisplayName(entity)}
      canonicalIdentifier={getEntityBadgeIdentifier(entity)}
      entityType={getEntityTypeLabel(entity)}
      resolutionStatus={entity.resolutionStatus}
    />
  </button>
{/snippet}

<div class="@container" role="table" aria-label="Relations">
  <div role="rowgroup">
    <div
      role="row"
      class="hidden items-center gap-x-2 border-b px-3 py-2 text-sm font-medium text-muted-foreground @lg:grid {WIDE}"
    >
      <span role="columnheader" class="[grid-area:from]">From</span>
      <span role="columnheader" class="text-center [grid-area:rel]">Relation</span>
      <span role="columnheader" class="[grid-area:to]">To</span>
      <span role="columnheader" class="text-center [grid-area:ev]">Evidence</span>
    </div>
  </div>
  <div role="rowgroup">
    {#each rows as row}
      <div
        role="row"
        tabindex="0"
        onclick={() => onRowClick(row)}
        onkeydown={(event) => {
          if (event.key === 'Enter' && event.target === event.currentTarget) onRowClick(row);
        }}
        class="grid cursor-pointer items-center gap-x-2 gap-y-1 border-b px-3 py-2 transition-colors hover:bg-muted/50 focus-visible:bg-muted/50 focus-visible:outline-none {NARROW} {WIDE}"
      >
        <div role="cell" class="min-w-0 [grid-area:from]">
          {@render endpoint(row.subjectEntity)}
        </div>
        <span aria-hidden="true" class="text-muted-foreground [grid-area:arrow] @lg:hidden">→</span>
        <div role="cell" class="min-w-0 [grid-area:to]">
          {@render endpoint(row.objectEntity)}
        </div>
        <div role="cell" class="min-w-0 justify-self-start [grid-area:rel] @lg:justify-self-center">
          {#if row.relation.predicate.trim()}
            <RelationLabel relation={row.relation} responsive />
          {:else}
            <Minus class="h-4 w-4 text-muted-foreground" />
          {/if}
        </div>
        <div role="cell" class="justify-self-end [grid-area:ev] @lg:justify-self-center">
          <Badge variant="outline" title="Evidence">
            {formatNumber(row.relation.evidenceCount || 0)}
          </Badge>
        </div>
      </div>
    {/each}
  </div>
  {#if hasMore}
    <div class="flex justify-center py-4" style="min-height: 40px;">
      <Button variant="outline" onclick={onLoadMore} disabled={loadingMore}>
        {#if loadingMore}
          <div
            class="mr-2 h-4 w-4 animate-spin rounded-full border-2 border-primary border-t-transparent"
          ></div>
          <span>Loading...</span>
        {:else}
          <span>Load more</span>
        {/if}
      </Button>
    </div>
  {/if}
</div>
