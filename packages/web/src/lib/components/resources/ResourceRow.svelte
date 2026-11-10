<script lang="ts">
  import { Download, Telescope } from '@lucide/svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import * as Table from '$lib/components/ui/table/index.js';
  import { truncateTitle } from '$lib/actions/truncate-title';
  import { formatBytes } from '$lib/resources/catalog';
  import ResourceDonut from './ResourceDonut.svelte';
  import type { ResourceRecord } from '$lib/resources/types';

  let {
    resource,
    selectedTags = [],
    onTagClick,
    onOpen,
  }: {
    resource: ResourceRecord;
    selectedTags?: string[];
    onTagClick?: (id: string) => void;
    onOpen: (resource: ResourceRecord) => void;
  } = $props();
  const countFormatter = new Intl.NumberFormat('en-US');
  const files = $derived(
    (resource.files || []).filter((file) => file.name !== 'evidence_payloads.parquet'),
  );
  const totalBytes = $derived(
    resource.total_size_bytes || files.reduce((sum, file) => sum + (file.size_bytes || 0), 0),
  );
  const totalRecords = $derived(resource.entity_count + resource.interaction_count);
  const entityPercent = $derived(
    totalRecords > 0 ? (resource.entity_count / totalRecords) * 100 : 0,
  );
  const description = $derived(
    resource.description ||
      (resource.interaction_count
        ? 'Molecular entities and relationships available to explore and download.'
        : 'A reference catalog of molecular entities and annotations.'),
  );
  const topics = $derived((resource.tags ?? []).filter((tag) => tag.dimension === 'topic'));

  const stop = (event: Event) => event.stopPropagation();
</script>

<Table.Row class="cursor-pointer" onclick={() => onOpen(resource)}>
  <Table.Cell class="max-w-0 py-2 pl-3">
    <div class="flex min-w-0 items-center gap-3">
      <ResourceDonut {resource} />
      <div class="min-w-0">
        <!-- The row opens the details; the name is its button, for the keyboard. Explore
             is its own labelled button, so no part of the row leaves the page unexpectedly. -->
        <button
          type="button"
          class="block max-w-full truncate text-left font-semibold text-foreground hover:underline"
          onclick={(event) => {
            stop(event);
            onOpen(resource);
          }}
          use:truncateTitle={resource.resource_name}>{resource.resource_name}</button
        >
        <p class="truncate text-xs text-muted-foreground" use:truncateTitle={description}>
          {description}
        </p>
      </div>
    </div>
  </Table.Cell>
  <Table.Cell class="text-right tabular-nums"
    >{countFormatter.format(resource.entity_count)}<span
      class="block text-[11px] text-muted-foreground">{entityPercent.toFixed(0)}%</span
    ></Table.Cell
  >
  <!-- Relations and size need the width a phone does not have; the dialog shows them. -->
  <Table.Cell class="hidden text-right tabular-nums sm:table-cell"
    >{countFormatter.format(resource.interaction_count)}<span
      class="block text-[11px] text-muted-foreground"
      >{totalRecords ? (100 - entityPercent).toFixed(0) : 0}%</span
    ></Table.Cell
  >
  <Table.Cell class="hidden text-right tabular-nums sm:table-cell"
    >{formatBytes(totalBytes)}</Table.Cell
  >
  <Table.Cell class="hidden xl:table-cell">
    <div class="flex flex-wrap gap-1">
      {#each topics as tag (tag.id)}<button
          type="button"
          class={`rounded-md border px-1.5 py-0.5 text-[11px] transition-colors hover:border-ring ${
            selectedTags.includes(tag.id) ? 'border-primary bg-primary/10' : 'bg-muted/40'
          }`}
          aria-pressed={selectedTags.includes(tag.id)}
          title={tag.description}
          onclick={(event) => {
            stop(event);
            onTagClick?.(tag.id);
          }}>{tag.label}</button
        >{/each}
    </div>
  </Table.Cell>
  <Table.Cell class="pr-3">
    <div class="flex items-center justify-end gap-1">
      <Button
        variant="outline"
        size="sm"
        class="h-8 gap-1.5 px-2"
        href={`/explore?source=${encodeURIComponent(resource.resource_id)}`}
        aria-label={`Explore the entities of ${resource.resource_name}`}
        title="Explore its entities"
        onclick={stop}
        ><Telescope class="size-4" /><span class="hidden sm:inline">Explore</span></Button
      >
      <Button
        variant="outline"
        size="icon-sm"
        href={`/api/resources/${encodeURIComponent(resource.resource_id)}/download`}
        download={`${resource.resource_id}.zip`}
        aria-label={`Download ${resource.resource_name} ZIP`}
        title="Download resource ZIP"
        onclick={stop}><Download class="size-4" /></Button
      >
    </div>
  </Table.Cell>
</Table.Row>
