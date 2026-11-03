<script lang="ts">
  import { releaseUrl } from '$lib/api/release';
  import { Download, ExternalLink } from '@lucide/svelte';
  import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogHeader,
    DialogTitle,
  } from '$lib/components/ui/dialog/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import ResourceDonut from './ResourceDonut.svelte';
  import { formatBytes } from '$lib/resources/catalog';
  import type { ResourceFile, ResourceRecord, LicenseUseStatus } from '$lib/resources/types';

  let {
    open = $bindable(false),
    resource,
    selectedTags = [],
    onTagClick,
  }: {
    open: boolean;
    resource: ResourceRecord | null;
    selectedTags?: string[];
    onTagClick?: (id: string) => void;
  } = $props();

  const countFormatter = new Intl.NumberFormat('en-US');
  const files = $derived(
    (resource?.files || []).filter((file) => file.name !== 'evidence_payloads.parquet'),
  );
  const totalBytes = $derived(
    resource?.total_size_bytes || files.reduce((sum, file) => sum + (file.size_bytes || 0), 0),
  );
  const totalRecords = $derived((resource?.entity_count ?? 0) + (resource?.interaction_count ?? 0));
  const entityPercent = $derived(
    totalRecords > 0 ? ((resource?.entity_count ?? 0) / totalRecords) * 100 : 0,
  );

  function safeUrl(value: string | null | undefined): string | undefined {
    if (!value) return undefined;
    try {
      const url = new URL(value);
      return ['http:', 'https:'].includes(url.protocol) ? url.href : undefined;
    } catch {
      return undefined;
    }
  }
  function fileHref(id: string, file: ResourceFile): string {
    return (
      safeUrl(file.url) ??
      releaseUrl(`/api/resources/${encodeURIComponent(id)}/files/${encodeURIComponent(file.name)}`)
    );
  }
  function fileLabel(file: ResourceFile): string {
    return file.name.replace(/\.parquet$/, '').replaceAll('_', ' ');
  }
  function licenseLabel(status?: LicenseUseStatus): string {
    return {
      allowed: 'Allowed',
      requires_permission: 'Requires permission',
      not_allowed: 'Not allowed',
      unknown: 'Unknown',
    }[status ?? 'unknown'];
  }
</script>

<Dialog bind:open>
  <DialogContent class="max-h-[88dvh] gap-0 overflow-y-auto p-0 sm:max-w-2xl">
    {#if resource}
      <DialogHeader class="gap-1.5 px-6 pt-5 pb-4 pr-14 text-left">
        <DialogTitle class="text-2xl font-semibold tracking-tight"
          >{resource.resource_name}</DialogTitle
        >
        <DialogDescription>{resource.description}</DialogDescription>
        <div class="flex flex-wrap gap-2 pt-2">
          <Button
            size="sm"
            href={releaseUrl(`/api/resources/${encodeURIComponent(resource.resource_id)}/download`)}
            download={`${resource.resource_id}.zip`}><Download class="size-4" />Download ZIP</Button
          >
          <Button
            size="sm"
            variant="outline"
            href={releaseUrl(`/explore?source=${encodeURIComponent(resource.resource_id)}`)}
            >Explore its entities</Button
          >
        </div>
      </DialogHeader>

      <div class="grid gap-6 border-t px-6 py-5 sm:grid-cols-[auto_minmax(0,1fr)]">
        <div class="flex items-center gap-4">
          <ResourceDonut {resource} size={100} ring={11}>
            <span class="text-center leading-tight"
              ><strong class="block text-sm tabular-nums">{formatBytes(totalBytes)}</strong><span
                class="text-[10px] text-muted-foreground">download</span
              ></span
            >
          </ResourceDonut>
          <dl class="grid grid-cols-[auto_auto_auto] items-baseline gap-x-3 gap-y-1.5 text-sm">
            <dt class="flex items-center gap-1.5 text-muted-foreground">
              <span class="size-2 rounded-sm bg-[var(--resource-entities)]"></span>Entities
            </dt>
            <dd class="text-right font-semibold tabular-nums">
              {countFormatter.format(resource.entity_count)}
            </dd>
            <dd class="text-xs text-muted-foreground tabular-nums">
              {entityPercent.toFixed(0)}%
            </dd>
            <dt class="flex items-center gap-1.5 text-muted-foreground">
              <span class="size-2 rounded-sm bg-[var(--resource-relations)]"></span>Relations
            </dt>
            <dd class="text-right font-semibold tabular-nums">
              {countFormatter.format(resource.interaction_count)}
            </dd>
            <dd class="text-xs text-muted-foreground tabular-nums">
              {totalRecords ? (100 - entityPercent).toFixed(0) : 0}%
            </dd>
          </dl>
        </div>
        <div class="min-w-0 space-y-1.5">
          <h3 class="text-xs font-semibold tracking-wide text-muted-foreground uppercase">Files</h3>
          <ul class="space-y-1 text-sm">
            {#each files as file (file.name)}
              <li class="flex items-center justify-between gap-3">
                <a
                  class="inline-flex min-w-0 items-center gap-1.5 text-foreground underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                  href={fileHref(resource.resource_id, file)}
                  download={`${resource.resource_id}-${file.name}`}
                  ><Download class="size-3.5 shrink-0" /><span class="truncate"
                    >{fileLabel(file)}</span
                  ></a
                >
                <span class="shrink-0 text-xs text-muted-foreground tabular-nums"
                  >{formatBytes(file.size_bytes)}</span
                >
              </li>
            {/each}
          </ul>
          {#if resource.dataset_subset?.length}<p class="text-xs text-muted-foreground">
              Selected datasets: {resource.dataset_subset.join(', ')}
            </p>{/if}
        </div>
      </div>

      <div class="space-y-4 border-t px-6 py-5 text-sm">
        <dl class="grid grid-cols-[max-content_minmax(0,1fr)] gap-x-5 gap-y-2">
          <dt class="text-muted-foreground">Version</dt>
          <dd class="tabular-nums">{resource.version ?? 'Not recorded'}</dd>
          <dt class="text-muted-foreground">Website</dt>
          <dd class="break-words">
            {#if safeUrl(resource.website)}<a
                class="inline-flex items-center gap-1 text-foreground underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                href={safeUrl(resource.website)}
                target="_blank"
                rel="noopener noreferrer">{resource.website}<ExternalLink class="size-3" /></a
              >{:else}Not specified{/if}
          </dd>
          <dt class="text-muted-foreground">License</dt>
          <dd class="break-words">
            {#if safeUrl(resource.license_url)}<a
                class="inline-flex items-center gap-1 text-foreground underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                href={safeUrl(resource.license_url)}
                target="_blank"
                rel="noopener noreferrer"
                >{resource.license || 'License terms'}<ExternalLink class="size-3" /></a
              >{:else}{resource.license || 'Not specified'}{/if}
          </dd>
          <dt class="text-muted-foreground">License use</dt>
          <dd>
            Academic: {licenseLabel(resource.license_use?.academic)} · Commercial: {licenseLabel(
              resource.license_use?.commercial,
            )}
          </dd>
          <dt class="text-muted-foreground">pypath input module</dt>
          <dd>
            {#if safeUrl(resource.input_module_url)}<a
                class="inline-flex items-center gap-1 text-foreground underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                href={safeUrl(resource.input_module_url)}
                target="_blank"
                rel="noopener noreferrer">View on GitHub<ExternalLink class="size-3" /></a
              >{:else}Not recorded{/if}
          </dd>
        </dl>
        {#if resource.tags?.length}
          <div class="space-y-1.5">
            <h3 class="text-xs font-semibold tracking-wide text-muted-foreground uppercase">
              Tags
            </h3>
            <div class="flex flex-wrap gap-1">
              {#each resource.tags as tag (tag.id)}<button
                  type="button"
                  class={`rounded-md border px-1.5 py-0.5 text-xs transition-colors hover:border-ring ${
                    selectedTags.includes(tag.id) ? 'border-primary bg-primary/10' : 'bg-muted/40'
                  }`}
                  aria-pressed={selectedTags.includes(tag.id)}
                  title={tag.description
                    ? `${tag.description} Click to filter the list.`
                    : 'Filter the list by this tag'}
                  onclick={() => onTagClick?.(tag.id)}>{tag.label}</button
                >{/each}
            </div>
          </div>
        {/if}
      </div>
    {/if}
  </DialogContent>
</Dialog>
