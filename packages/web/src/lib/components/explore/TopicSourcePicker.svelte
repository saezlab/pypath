<script module lang="ts">
  import type { components } from '$lib/api/types.generated';

  type ResourceTopic = components['schemas']['ResourceTopic'];

  // The topic map of a release (each topic with its resources, a few KB), loaded once per
  // page session: every Sources panel (Entities, Relations, the mobile filter sheet) shares
  // it. A failed load is retried next time.
  const topicMaps = new Map<string, Promise<ResourceTopic[]>>();
  export function releaseTopics(
    release: string,
    load: () => Promise<ResourceTopic[]>,
  ): Promise<ResourceTopic[]> {
    let topics = topicMaps.get(release);
    if (!topics) {
      topics = load().catch((error) => {
        topicMaps.delete(release);
        throw error;
      });
      topicMaps.set(release, topics);
    }
    return topics;
  }
</script>

<script lang="ts">
  // Topic cards in the Sources filter of Explore. A card selects every source of
  // its topic that the current results have; clicking it again removes them. A card whose
  // sources are only partly selected shows how many.
  import { page } from '$app/state';
  import { releaseFetch } from '$lib/api/release';
  import { topicIcon } from '$lib/resources/topics';

  let {
    selected,
    available,
    onChange,
  }: {
    selected: string[];
    /** The sources the current view offers; a topic selects only these. */
    available: string[];
    onChange: (sources: string[]) => void;
  } = $props();

  let topicMap = $state<ResourceTopic[]>([]);
  $effect(() => {
    let cancelled = false;
    releaseTopics(page.data.selectedRelease || 'latest', async () => {
      const res = await releaseFetch('/app-api/resources/topics');
      if (!res.ok) throw new Error(`Resource topics returned HTTP ${res.status}`);
      return ((await res.json()) as components['schemas']['ResourceTopicsResponse']).topics ?? [];
    })
      .then((topics) => {
        if (!cancelled) topicMap = topics;
      })
      .catch(() => {});
    return () => {
      cancelled = true;
    };
  });

  const offered = $derived(new Set(available));
  const topics = $derived(
    topicMap
      .map((topic) => ({ topic, ids: (topic.resources ?? []).filter((id) => offered.has(id)) }))
      .filter(({ ids }) => ids.length > 0),
  );

  function toggle(ids: string[], on: boolean) {
    const next = on
      ? selected.filter((id) => !ids.includes(id))
      : [...new Set([...selected, ...ids])];
    onChange(next);
  }
</script>

{#if topics.length}
  <div class="mb-3 space-y-1.5">
    <h4 class="text-[11px] font-semibold tracking-[0.16em] text-muted-foreground uppercase">
      Topics
    </h4>
    <div class="flex flex-wrap gap-1.5">
      {#each topics as { topic, ids } (topic.id)}
        {@const Icon = topicIcon(topic.id)}
        {@const picked = ids.filter((id) => selected.includes(id)).length}
        {@const on = picked === ids.length}
        <button
          type="button"
          aria-pressed={on ? 'true' : picked ? 'mixed' : 'false'}
          aria-label={`${topic.label}: ${picked ? `${picked} of ` : ''}${ids.length} sources`}
          title={topic.description}
          onclick={() => toggle(ids, on)}
          class="group inline-flex max-w-full items-center gap-1.5 rounded-md border bg-background py-1 pr-2 pl-1 text-xs transition-colors hover:border-ring {on
            ? 'border-primary bg-primary/5 ring-1 ring-primary'
            : picked
              ? 'border-primary/50'
              : ''}"
        >
          <span
            class="grid size-5 shrink-0 place-items-center rounded {on
              ? 'bg-primary text-primary-foreground'
              : 'bg-muted text-muted-foreground group-hover:text-foreground'}"
            ><Icon class="size-3" /></span
          >
          <span class="truncate">{topic.label}</span>
          <span class="text-muted-foreground tabular-nums"
            >{picked && !on ? `${picked}/` : ''}{ids.length}</span
          >
        </button>
      {/each}
    </div>
  </div>
{/if}
