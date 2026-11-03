<script lang="ts">
  import EntityResultsList from '$lib/components/entity/EntityResultsList.svelte';
  import { entityFromWire } from '$lib/api/adapters';
  import type { WireEntity } from '$lib/api/contracts';
  import { untrack, type Snippet } from 'svelte';
  import { page } from '$app/state';
  import { releaseFetch } from '$lib/api/release';
  import { Button } from '$lib/components/ui/button/index.js';
  import type { SearchFilters } from '$lib/types/search';
  import type { EntityWithIdentifiers } from '$lib/types/entities';
  type Group = {
    group_key: string;
    group_label: string | null;
    is_group: boolean;
    member_count: number;
    entity: EntityWithIdentifiers;
    members: EntityWithIdentifiers[];
    nextMemberCursor: string | null;
  };
  let {
    query,
    filters,
    renderMember,
    strategy = 'auto',
    loading = $bindable(false),
  }: {
    strategy?: string;
    /** A page of groups is being loaded. */
    loading?: boolean;
    query: string;
    filters: SearchFilters;
    renderMember: Snippet<[EntityWithIdentifiers]>;
  } = $props();
  let groups = $state<Group[]>([]);
  let nextCursor = $state<string | null>(null);
  let error = $state('');
  let generation = 0;
  const requestKey = $derived(
    JSON.stringify([query, filters, strategy, page.data.selectedRelease]),
  );
  let controller: AbortController | undefined;

  async function request(
    extra: Record<string, unknown>,
    q = query,
    f = filters,
    selectedStrategy = strategy,
  ) {
    const response = await releaseFetch('/app-api/entities/groups', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      signal: controller?.signal,
      body: JSON.stringify({
        strategy: selectedStrategy,
        query: q,
        filters: f,
        member_limit: 1,
        ...extra,
      }),
    });
    if (!response.ok) throw new Error('Unable to load groups. Please retry.');
    const data = (await response.json()) as {
      groups: Array<
        Omit<Group, 'entity' | 'members'> & { entity: WireEntity; members: WireEntity[] }
      >;
      nextCursor: string | null;
    };
    return {
      ...data,
      groups: data.groups.map((group) => ({
        ...group,
        entity: entityFromWire(group.entity),
        members: group.members.map(entityFromWire),
      })),
    };
  }

  $effect(() => {
    const [q, f, selectedStrategy] = JSON.parse(requestKey);
    controller = new AbortController();
    const current = ++generation;
    // The current groups stay, dimmed by the panel, until the new page replaces them.
    loading = true;
    error = '';
    untrack(() => request({}, q, f, selectedStrategy))
      .then((data) => {
        if (current !== generation) return;
        groups = data.groups;
        nextCursor = data.nextCursor;
      })
      .catch((e) => {
        if (current !== generation) return;
        groups = [];
        nextCursor = null;
        error = e.message;
      })
      .finally(() => {
        if (current === generation) loading = false;
      });
    return () => {
      generation++;
      controller?.abort();
    };
  });

  async function moreGroups() {
    if (loading || !nextCursor) return;
    const current = generation;
    loading = true;
    error = '';
    try {
      const data = await request({ cursor: nextCursor });
      if (current !== generation) return;
      groups = [...groups, ...data.groups];
      nextCursor = data.nextCursor;
    } catch (e) {
      if (current === generation) error = (e as Error).message;
    } finally {
      if (current === generation) loading = false;
    }
  }
</script>

{#if error}<p role="alert" class="p-4 text-destructive">{error}</p>{/if}
{#if loading && !groups.length}<p role="status" class="p-4 text-sm text-muted-foreground">
    Loading results…
  </p>{/if}
{#if groups.length}
  <EntityResultsList>
    {#each groups as group (group.group_key)}
      {@render renderMember(group.entity)}
    {/each}
  </EntityResultsList>
{/if}
{#if !loading && !error && !groups.length}<p class="p-4 text-sm text-muted-foreground">
    No matching entities found.
  </p>{/if}
{#if nextCursor}<div class="flex justify-center py-6">
    <Button variant="outline" disabled={loading} onclick={moreGroups}>Load more groups</Button>
  </div>{/if}
