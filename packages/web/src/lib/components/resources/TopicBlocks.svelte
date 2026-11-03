<script lang="ts">
  // Topic cards above the resource list. Selecting a topic sets the same `tag`
  // filter the sidebar used to offer, so the list, the counts and the footer follow.
  import { Download } from '@lucide/svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import type { ResourceRecord, ResourceTag } from '$lib/resources/types';
  import { formatBytes } from '$lib/resources/catalog';
  import { catalogTopics, topicIcon } from '$lib/resources/topics';

  let {
    resources,
    selectedTags,
    onToggle,
    onClear,
    onDownload,
  }: {
    resources: ResourceRecord[];
    selectedTags: string[];
    onToggle: (id: string) => void;
    onClear: (ids: string[]) => void;
    onDownload: (topics: ResourceTag[]) => void;
  } = $props();

  const topics = $derived(catalogTopics(resources));
  const chosen = $derived(topics.filter((topic) => selectedTags.includes(topic.tag.id)));
  // The resources of the chosen topics, each once (a resource can have several topics).
  const chosenResources = $derived(
    new Map(chosen.flatMap((topic) => topic.members.map((r) => [r.resource_id, r]))),
  );
  const chosenBytes = $derived(
    [...chosenResources.values()].reduce((sum, r) => sum + (r.total_size_bytes || 0), 0),
  );
</script>

<section aria-label="Topics" class="shrink-0 space-y-2">
  <div class="grid grid-cols-2 gap-2 sm:grid-cols-3 md:grid-cols-5">
    {#each topics as topic (topic.tag.id)}
      {@const Icon = topicIcon(topic.tag.id)}
      {@const on = selectedTags.includes(topic.tag.id)}
      <button
        type="button"
        aria-pressed={on}
        title={topic.tag.description}
        onclick={() => onToggle(topic.tag.id)}
        class="group flex min-w-0 items-center gap-2.5 rounded-lg border bg-background px-3 py-2 text-left transition-colors hover:border-ring {on
          ? 'border-primary bg-primary/5 ring-1 ring-primary'
          : ''}"
      >
        <span
          class="grid size-8 shrink-0 place-items-center rounded-md {on
            ? 'bg-primary text-primary-foreground'
            : 'bg-muted text-muted-foreground group-hover:text-foreground'}"
          ><Icon class="size-4" /></span
        >
        <span class="min-w-0">
          <span class="block truncate text-sm font-medium">{topic.tag.label}</span>
          <span class="block text-xs text-muted-foreground tabular-nums"
            >{topic.members.length} resources</span
          >
        </span>
      </button>
    {/each}
  </div>
  {#if chosen.length}
    <div
      class="flex flex-wrap items-center justify-between gap-2 rounded-lg border bg-muted/20 px-3 py-2 text-sm"
    >
      <span>
        <span class="font-medium">{chosen.map((topic) => topic.tag.label).join(' or ')}</span>
        <span class="text-muted-foreground tabular-nums">
          · {chosenResources.size} resources · {formatBytes(chosenBytes)}</span
        >
      </span>
      <span class="flex items-center gap-1.5">
        <Button
          variant="ghost"
          size="sm"
          class="text-muted-foreground"
          onclick={() => onClear(chosen.map((topic) => topic.tag.id))}>Clear</Button
        >
        <Button size="sm" onclick={() => onDownload(chosen.map((topic) => topic.tag))}
          ><Download class="size-4" />Download selection</Button
        >
      </span>
    </div>
  {/if}
</section>
