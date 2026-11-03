<script lang="ts">
  import { Download, FileDown, Terminal } from '@lucide/svelte';
  import { toast } from 'svelte-sonner';
  import { page } from '$app/state';
  import { releaseFetch } from '$lib/api/release';
  import {
    DropdownMenu,
    DropdownMenuContent,
    DropdownMenuItem,
    DropdownMenuLabel,
    DropdownMenuSeparator,
    DropdownMenuTrigger,
  } from '$lib/components/ui/dropdown-menu/index.js';

  /** Exports what a results panel shows as Parquet: a download, or the API request for it. */
  let {
    endpoint,
    body,
    name,
    limit,
    disabledReason,
  }: {
    /** API path of the export, e.g. '/export'. */
    endpoint: string;
    /** The request body: the filters and entity or group keys of the current results. */
    body: Record<string, unknown>;
    /** What is exported, e.g. 'relations'; names the file. */
    name: string;
    /** Most rows one export returns. */
    limit: number;
    /** Why the current results cannot be exported; the menu explains it instead. */
    disabledReason?: string;
  } = $props();

  let downloading = $state(false);
  const release = $derived(page.data.selectedRelease || 'latest');
  const filename = $derived(`omnipath_${name}.parquet`);

  async function download() {
    if (downloading) return;
    downloading = true;
    try {
      const response = await releaseFetch(`/app-api${endpoint}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(body),
      });
      if (!response.ok) throw new Error(`Export failed (${response.status})`);
      const url = URL.createObjectURL(await response.blob());
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      link.click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    } catch (error) {
      toast.error((error as Error).message || 'Export failed');
    } finally {
      downloading = false;
    }
  }

  async function copyCommand() {
    const json = JSON.stringify(body).replaceAll("'", "'\\''");
    const url = `${page.url.origin}/api${endpoint}?release=${encodeURIComponent(release)}`;
    const command = `curl -fsSL -X POST '${url}' -H 'Content-Type: application/json' -d '${json}' -o ${filename}`;
    try {
      await navigator.clipboard.writeText(command);
      toast.success('API command copied');
    } catch {
      toast.error('Could not copy the command');
    }
  }
</script>

<DropdownMenu>
  <DropdownMenuTrigger
    class="inline-flex h-6 items-center gap-1.5 rounded-md px-2 text-xs text-muted-foreground transition-colors hover:bg-muted hover:text-foreground focus-visible:ring-2 focus-visible:ring-ring focus-visible:outline-none"
    aria-label={`Export ${name}`}
  >
    {#if downloading}<span class="export-spinner" aria-hidden="true"></span>{:else}<Download
        class="size-3.5"
        aria-hidden="true"
      />{/if}
    Export
  </DropdownMenuTrigger>
  <DropdownMenuContent align="end" class="w-64">
    <DropdownMenuLabel class="text-xs font-normal text-muted-foreground">
      {disabledReason ?? `Current results as Parquet, up to ${limit.toLocaleString('en')} ${name}`}
    </DropdownMenuLabel>
    <DropdownMenuSeparator />
    <DropdownMenuItem onSelect={download} disabled={downloading || !!disabledReason}>
      <FileDown class="size-4" />Download Parquet
    </DropdownMenuItem>
    <DropdownMenuItem onSelect={copyCommand} disabled={!!disabledReason}>
      <Terminal class="size-4" />Copy API command
    </DropdownMenuItem>
  </DropdownMenuContent>
</DropdownMenu>

<style>
  .export-spinner {
    width: 0.8rem;
    height: 0.8rem;
    border-radius: 9999px;
    border: 2px solid currentColor;
    border-right-color: transparent;
    animation: export-spin 0.7s linear infinite;
  }
  @keyframes export-spin {
    to {
      transform: rotate(360deg);
    }
  }
</style>
