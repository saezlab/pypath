<script lang="ts">
  // Download every resource of one or more topics, as a ZIP of Parquet files or
  // through the Python client (which installs them into a local DuckDB database).
  import { page } from '$app/state';
  import { Check, Copy, Download } from '@lucide/svelte';
  import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogHeader,
    DialogTitle,
  } from '$lib/components/ui/dialog/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { formatBytes } from '$lib/resources/catalog';
  import type { ResourceRecord, ResourceTag } from '$lib/resources/types';

  let {
    open = $bindable(false),
    topics,
    resources,
  }: { open: boolean; topics: ResourceTag[]; resources: ResourceRecord[] } = $props();

  const ids = $derived(new Set(topics.map((topic) => topic.id)));
  const members = $derived(
    resources
      .filter((r) => (r.tags ?? []).some((tag) => ids.has(tag.id)))
      .sort((a, b) => a.resource_name.localeCompare(b.resource_name)),
  );
  const bytes = $derived(members.reduce((sum, r) => sum + (r.total_size_bytes || 0), 0));
  const label = $derived(topics.map((topic) => topic.label).join(' or '));
  const release = $derived(page.data.selectedRelease || 'latest');
  const python = $derived(
    `import omnipath_client\n\nomnipath_client.init(topic=${
      topics.length === 1 ? `"${topics[0]?.id}"` : `[${topics.map((t) => `"${t.id}"`).join(', ')}]`
    })\ndb = omnipath_client.connect()`,
  );
  const shell = $derived(`omnipath-client init --topic ${topics.map((t) => t.id).join(' ')}`);
  let copied = $state('');
  let busy = $state(false);
  let error = $state('');

  async function copy(text: string, key: string) {
    await navigator.clipboard.writeText(text);
    copied = key;
    setTimeout(() => (copied = ''), 1500);
  }
  // The bulk endpoint takes the resource ids in a POST body, so the ZIP is fetched and handed
  // to the browser; a production version would stream a GET link instead (see the note).
  async function downloadZip() {
    busy = true;
    error = '';
    try {
      const res = await fetch(`/api/resources/download?release=${encodeURIComponent(release)}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ resource_ids: members.map((r) => r.resource_id) }),
      });
      if (!res.ok) throw new Error(`The download failed (${res.status})`);
      const url = URL.createObjectURL(await res.blob());
      const link = document.createElement('a');
      link.href = url;
      link.download = `omnipath-${topics.map((t) => t.id).join('-')}.zip`;
      link.click();
      URL.revokeObjectURL(url);
    } catch (exc) {
      error = exc instanceof Error ? exc.message : String(exc);
    } finally {
      busy = false;
    }
  }
</script>

<Dialog bind:open>
  <DialogContent class="sm:max-w-xl">
    <DialogHeader>
      <DialogTitle>Download {label}</DialogTitle>
      <DialogDescription>
        {members.length} resources · {formatBytes(bytes)} of Parquet tables · release {release}
      </DialogDescription>
    </DialogHeader>

    <div class="max-h-32 overflow-y-auto rounded-md border px-3 py-2 text-xs leading-5">
      {#each members as resource, index (resource.resource_id)}<span
          >{resource.resource_name}{index < members.length - 1 ? ', ' : ''}</span
        >{/each}
    </div>

    <div class="space-y-3">
      <section class="space-y-1.5">
        <h3 class="text-sm font-medium">Python client</h3>
        <p class="text-xs text-muted-foreground">
          Installs the resources into a local DuckDB database you can query with SQL.
        </p>
        {@render snippet(python, 'python')}
        {@render snippet(shell, 'shell')}
      </section>
      <section class="flex items-center justify-between gap-3 border-t pt-3">
        <div>
          <h3 class="text-sm font-medium">ZIP of Parquet files</h3>
          <p class="text-xs text-muted-foreground">
            One folder per resource, without raw evidence.
          </p>
        </div>
        <Button onclick={downloadZip} disabled={busy || !members.length}
          ><Download class="size-4" />{busy
            ? 'Preparing…'
            : `Download ${formatBytes(bytes)}`}</Button
        >
      </section>
      {#if error}<p class="text-xs text-destructive">{error}</p>{/if}
    </div>
  </DialogContent>
</Dialog>

{#snippet snippet(text: string, key: string)}
  <div class="relative rounded-md border bg-muted/40">
    <pre class="overflow-x-auto px-3 py-2 pr-10 font-mono text-xs leading-5">{text}</pre>
    <button
      type="button"
      aria-label="Copy"
      class="absolute top-1.5 right-1.5 rounded p-1 text-muted-foreground hover:bg-background hover:text-foreground"
      onclick={() => copy(text, key)}
      >{#if copied === key}<Check class="size-3.5" />{:else}<Copy class="size-3.5" />{/if}</button
    >
  </div>
{/snippet}
