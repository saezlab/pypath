<script lang="ts">
  import { untrack } from 'svelte';
  import { page } from '$app/state';
  import { Alert, AlertDescription } from '$lib/components/ui/alert/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import RelationsTable from '$lib/components/interactions/RelationsTable.svelte';
  import InteractionDetailsSheet from '$lib/components/interactions/InteractionDetailsSheet.svelte';
  import { fetchRelationsSearch } from '$lib/api/client';
  import type { EntityLike } from '$lib/domain/display';
  import type { InteractionListRow } from '$lib/types/interactions';
  import { ALL_SCOPE, type RelationScope } from '$lib/utils/entity-overview';
  import { formatNumber } from '$lib/utils/format';

  interface Props {
    scopes: RelationScope[];
    scopeId: string;
    totals?: Record<string, number>;
    onEntitySelect: (entity: EntityLike) => void;
  }

  let { scopes, scopeId = $bindable(ALL_SCOPE), totals = {}, onEntitySelect }: Props = $props();

  const PAGE_SIZE = 20;
  let rows = $state<InteractionListRow[]>([]);
  let loading = $state(false);
  let loadingMore = $state(false);
  let hasMore = $state(false);
  let error = $state<string | null>(null);
  let selected = $state<InteractionListRow | null>(null);
  let detailsOpen = $state(false);

  const products = $derived(
    scopes
      .filter((item) => item.kind === 'product')
      .map((item) => ({
        entityPk: item.id,
        label: item.title,
        canonicalIdentifier: item.identifier,
      })),
  );
  const scope = $derived(scopes.find((item) => item.id === scopeId) ?? scopes[0]);
  const searchKey = $derived(JSON.stringify([page.data.selectedRelease, scope?.filters ?? null]));

  $effect(() => {
    const [, filters] = JSON.parse(searchKey);
    rows = [];
    hasMore = false;
    error = null;
    if (!filters) return;
    const controller = new AbortController();
    loading = true;
    untrack(() => fetchRelationsSearch({ filters, limit: PAGE_SIZE, offset: 0 }, controller.signal))
      .then((data) => {
        rows = data.rows;
        hasMore = !!data.nextCursor;
      })
      .catch((err) => {
        if (!controller.signal.aborted)
          error = err instanceof Error ? err.message : 'Relations could not be loaded.';
      })
      .finally(() => {
        if (!controller.signal.aborted) loading = false;
      });
    return () => controller.abort();
  });

  async function loadMore() {
    if (!scope || loadingMore || !hasMore) return;
    const key = searchKey;
    loadingMore = true;
    try {
      const data = await fetchRelationsSearch({
        filters: scope.filters,
        limit: PAGE_SIZE,
        offset: rows.length,
      });
      if (key !== searchKey) return;
      rows = [...rows, ...data.rows];
      hasMore = !!data.nextCursor;
    } catch (err) {
      if (key === searchKey)
        error = err instanceof Error ? err.message : 'More relations could not be loaded.';
    } finally {
      if (key === searchKey) loadingMore = false;
    }
  }
</script>

<div class="space-y-3">
  {#if scopes.length > 1}
    <div class="space-y-2">
      <div class="flex flex-wrap items-center gap-1.5" role="group" aria-label="Relation scope">
        {#each scopes as item (item.id)}
          <Button
            size="sm"
            variant={item.id === scope?.id ? 'secondary' : 'outline'}
            aria-pressed={item.id === scope?.id}
            title={item.title}
            onclick={() => (scopeId = item.id)}
          >
            {item.label}
            {#if totals[item.id] !== undefined}<span class="tabular-nums text-muted-foreground"
                >{formatNumber(totals[item.id])}</span
              >{/if}
          </Button>
        {/each}
      </div>
      <p class="text-xs text-muted-foreground">
        {#if !scope || scope.kind === 'all'}
          All relations of this entity and the entities grouped with it.
        {:else if scope.kind === 'product'}
          Relations whose evidence names {scope.title}{scope.identifier &&
          scope.identifier !== scope.title
            ? ` (${scope.identifier})`
            : ''}.
        {:else}
          Relations of {scope.title} only.
        {/if}
      </p>
    </div>
  {/if}

  {#if error}
    <Alert variant="destructive"><AlertDescription>{error}</AlertDescription></Alert>
  {:else if loading && !rows.length}
    <p role="status" class="text-sm text-muted-foreground">Loading relations…</p>
  {:else if rows.length}
    <div class="overflow-hidden rounded-lg border">
      <RelationsTable
        {rows}
        {hasMore}
        {loadingMore}
        onLoadMore={loadMore}
        onRowClick={(row) => {
          selected = row;
          detailsOpen = true;
        }}
        onEntityClick={onEntitySelect}
      />
    </div>
  {:else}
    <p class="text-sm text-muted-foreground">No relations are recorded for this entity.</p>
  {/if}
</div>

<InteractionDetailsSheet
  open={detailsOpen}
  onOpenChange={(open) => (detailsOpen = open)}
  interaction={selected}
  filters={scope?.kind === 'product' ? scope.filters : undefined}
  {products}
/>
