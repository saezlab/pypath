<script lang="ts">
  import { compareResourceVersions } from '$lib/utils/versions';
  import {
    AlertTriangle,
    Calendar,
    CheckCircle2,
    Database,
    FileText,
    HardDrive,
    Loader2,
    Trash2,
  } from '@lucide/svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogHeader,
    DialogTitle,
  } from '$lib/components/ui/dialog/index.js';

  export interface VersionItem {
    version: string;
    resource: string;
    entities: { rows: number; size_bytes: number; mtime: number | null; exists: boolean };
    relations: { rows: number; size_bytes: number; mtime: number | null; exists: boolean };
    payloads: { rows: number; size_bytes: number; mtime: number | null; exists: boolean };
    log?: { exists: boolean; size_bytes?: number; mtime?: number | null };
  }

  interface Props {
    open?: boolean;
    resourceName: string;
    versions: VersionItem[];
    activeVersion: string | null;
    onDeleteVersion: (resource: string, version: string) => Promise<boolean | void>;
    onCleanupOlder?: (resource: string) => Promise<boolean | void>;
    onViewLogs: (resource: string, version: string) => void;
  }

  let {
    open = $bindable(false),
    resourceName,
    versions = [],
    activeVersion,
    onDeleteVersion,
    onCleanupOlder,
    onViewLogs,
  }: Props = $props();

  let confirmDeleteVersion = $state<string | null>(null);
  let deletingVersion = $state<string | null>(null);
  let confirmCleanup = $state(false);
  let cleaningUp = $state(false);
  let error = $state<string | null>(null);

  function formatBytes(bytes: number | null | undefined): string {
    if (bytes == null || bytes === 0) return '0 B';
    const units = ['B', 'KB', 'MB', 'GB', 'TB'];
    const i = Math.floor(Math.log(bytes) / Math.log(1024));
    return `${(bytes / Math.pow(1024, i)).toFixed(i === 0 ? 0 : 1)} ${units[i]}`;
  }

  function formatNumber(num: number | null | undefined): string {
    if (num == null) return '0';
    return num.toLocaleString();
  }

  function versionTimestamp(item: VersionItem): number {
    return Math.max(item.entities.mtime || 0, item.relations.mtime || 0, item.payloads.mtime || 0);
  }

  function formatDate(seconds: number | null | undefined): string {
    if (!seconds) return 'Unknown date';
    const date = new Date(seconds * 1000);
    return date.toLocaleString(undefined, {
      month: 'short',
      day: 'numeric',
      year: 'numeric',
      hour: '2-digit',
      minute: '2-digit',
    });
  }

  function versionTotalSize(item: VersionItem): number {
    return (
      (item.entities.size_bytes || 0) +
      (item.relations.size_bytes || 0) +
      (item.payloads.size_bytes || 0)
    );
  }

  const sortedVersions = $derived.by(() => {
    return [...versions].sort((a, b) =>
      compareResourceVersions(b.version, a.version, versionTimestamp(b), versionTimestamp(a)),
    );
  });

  const totalResourceBytes = $derived.by(() => {
    return versions.reduce((sum, v) => sum + versionTotalSize(v), 0);
  });

  const olderVersionsCount = $derived.by(() => {
    return activeVersion
      ? versions.filter((v) => v.version !== activeVersion).length
      : Math.max(0, versions.length - 1);
  });

  async function handleDelete(ver: string) {
    deletingVersion = ver;
    error = null;
    try {
      await onDeleteVersion(resourceName, ver);
      confirmDeleteVersion = null;
    } catch (err) {
      error = err instanceof Error ? err.message : String(err);
    } finally {
      deletingVersion = null;
    }
  }

  async function handleCleanupOlder() {
    if (!onCleanupOlder) return;
    cleaningUp = true;
    error = null;
    try {
      await onCleanupOlder(resourceName);
      confirmCleanup = false;
    } catch (err) {
      error = err instanceof Error ? err.message : String(err);
    } finally {
      cleaningUp = false;
    }
  }
</script>

<Dialog bind:open>
  <DialogContent
    class="grid max-h-[88vh] grid-rows-[auto_minmax(0,1fr)] gap-4 overflow-hidden p-6 sm:max-w-3xl"
  >
    <DialogHeader class="pr-6">
      <div class="flex flex-wrap items-start justify-between gap-3">
        <div>
          <div class="flex items-center gap-2">
            <Database class="size-5 text-primary" />
            <DialogTitle class="text-lg font-bold capitalize">{resourceName} Versions</DialogTitle>
          </div>
          <DialogDescription class="mt-1 text-xs text-muted-foreground">
            {versions.length} version{versions.length === 1 ? '' : 's'} stored on disk · Total: {formatBytes(
              totalResourceBytes,
            )}
          </DialogDescription>
        </div>

        {#if olderVersionsCount > 0 && onCleanupOlder}
          <Button
            variant="outline"
            size="sm"
            class="h-8 gap-1.5 text-xs text-muted-foreground hover:text-foreground"
            onclick={() => (confirmCleanup = !confirmCleanup)}
            disabled={cleaningUp}
          >
            <Trash2 class="size-3.5" />
            Clean older ({olderVersionsCount})
          </Button>
        {/if}
      </div>

      {#if error}
        <div
          class="mt-2 rounded-lg border border-destructive/40 bg-destructive/10 p-2.5 text-xs text-destructive"
        >
          <strong>Error:</strong>
          {error}
        </div>
      {/if}

      {#if confirmCleanup}
        <div
          class="mt-3 flex flex-wrap items-center justify-between gap-3 rounded-xl border border-destructive/40 bg-destructive/5 p-3 text-xs"
        >
          <div class="flex items-center gap-2 text-destructive">
            <AlertTriangle class="size-4 shrink-0" />
            <span>
              Delete <strong>{olderVersionsCount}</strong> older version{olderVersionsCount === 1
                ? ''
                : 's'} of {resourceName}? The active version will be preserved.
            </span>
          </div>
          <div class="flex items-center gap-2">
            <Button
              variant="destructive"
              size="sm"
              class="h-7 px-3 text-xs"
              disabled={cleaningUp}
              onclick={handleCleanupOlder}
            >
              {#if cleaningUp}
                <Loader2 class="mr-1 size-3 animate-spin" />
                Deleting…
              {:else}
                Yes, Clean Older
              {/if}
            </Button>
            <Button
              variant="ghost"
              size="sm"
              class="h-7 px-2.5 text-xs"
              onclick={() => (confirmCleanup = false)}
            >
              Cancel
            </Button>
          </div>
        </div>
      {/if}
    </DialogHeader>

    <div class="min-h-0 overflow-y-auto pr-1">
      {#if sortedVersions.length === 0}
        <div
          class="flex flex-col items-center justify-center py-12 text-center text-muted-foreground"
        >
          <Database class="mb-2 size-8 stroke-1 text-muted-foreground/50" />
          <p class="text-sm font-medium">No built versions</p>
          <p class="mt-1 text-xs text-muted-foreground">
            Build this resource to create parquet tables and view them here.
          </p>
        </div>
      {:else}
        <div class="flex flex-col gap-3">
          {#each sortedVersions as ver (ver.version)}
            {@const isActive = ver.version === activeVersion}
            {@const stamp = versionTimestamp(ver)}
            {@const totalSize = versionTotalSize(ver)}
            <div
              class="flex flex-col gap-3 rounded-xl border p-4 transition-colors {isActive
                ? 'border-primary/50 bg-primary/5'
                : 'bg-card/70'}"
            >
              <div class="flex flex-wrap items-center justify-between gap-2">
                <div class="flex flex-wrap items-center gap-2">
                  <span class="font-mono text-sm font-semibold">{ver.version}</span>
                  {#if isActive}
                    <span
                      class="flex items-center gap-1 rounded-full border border-emerald-500/30 bg-emerald-500/10 px-2 py-0.5 text-[11px] font-medium text-emerald-600"
                    >
                      <CheckCircle2 class="size-3" />
                      Active Version
                    </span>
                  {/if}
                  <span class="flex items-center gap-1 text-xs text-muted-foreground">
                    <Calendar class="size-3" />
                    {formatDate(stamp)}
                  </span>
                  <span class="text-xs text-muted-foreground">·</span>
                  <span class="flex items-center gap-1 font-mono text-xs text-muted-foreground">
                    <HardDrive class="size-3" />
                    {formatBytes(totalSize)}
                  </span>
                </div>

                <div class="flex items-center gap-2">
                  {#if ver.log?.exists}
                    <Button
                      variant="outline"
                      size="sm"
                      class="h-7 gap-1 px-2 text-xs"
                      onclick={() => onViewLogs(resourceName, ver.version)}
                    >
                      <FileText class="size-3" />
                      Log
                    </Button>
                  {/if}

                  {#if confirmDeleteVersion === ver.version}
                    <div
                      class="flex items-center gap-1.5 rounded-lg border border-destructive/40 bg-destructive/10 px-2 py-1 text-xs text-destructive"
                    >
                      <AlertTriangle class="size-3.5 shrink-0" />
                      <span class="text-[11px] font-medium">Delete?</span>
                      <Button
                        size="sm"
                        variant="destructive"
                        class="h-6 px-2 text-[11px]"
                        disabled={deletingVersion === ver.version}
                        onclick={() => handleDelete(ver.version)}
                      >
                        {#if deletingVersion === ver.version}
                          <Loader2 class="size-3 animate-spin" />
                        {:else}
                          Yes
                        {/if}
                      </Button>
                      <Button
                        size="sm"
                        variant="ghost"
                        class="h-6 px-1.5 text-[11px]"
                        onclick={() => (confirmDeleteVersion = null)}
                      >
                        Cancel
                      </Button>
                    </div>
                  {:else}
                    <Button
                      variant="ghost"
                      size="sm"
                      class="h-7 gap-1 px-2 text-xs text-destructive hover:bg-destructive/10 hover:text-destructive"
                      onclick={() => (confirmDeleteVersion = ver.version)}
                    >
                      <Trash2 class="size-3" />
                      Delete
                    </Button>
                  {/if}
                </div>
              </div>

              <!-- File rows breakdown -->
              <div class="grid grid-cols-1 gap-2 sm:grid-cols-3">
                <div class="flex flex-col rounded-lg border bg-background/60 px-3 py-2 text-xs">
                  <div class="flex items-center justify-between text-muted-foreground">
                    <span>entity.parquet</span>
                    <span class="font-mono text-[11px]">{formatBytes(ver.entities.size_bytes)}</span
                    >
                  </div>
                  <div class="mt-1 font-mono font-semibold text-foreground">
                    {ver.entities.exists ? `${formatNumber(ver.entities.rows)} rows` : 'missing'}
                  </div>
                </div>

                <div class="flex flex-col rounded-lg border bg-background/60 px-3 py-2 text-xs">
                  <div class="flex items-center justify-between text-muted-foreground">
                    <span>relation.parquet</span>
                    <span class="font-mono text-[11px]"
                      >{formatBytes(ver.relations.size_bytes)}</span
                    >
                  </div>
                  <div class="mt-1 font-mono font-semibold text-foreground">
                    {ver.relations.exists ? `${formatNumber(ver.relations.rows)} rows` : 'missing'}
                  </div>
                </div>

                <div class="flex flex-col rounded-lg border bg-background/60 px-3 py-2 text-xs">
                  <div class="flex items-center justify-between text-muted-foreground">
                    <span>evidence_payloads.parquet</span>
                    <span class="font-mono text-[11px]">{formatBytes(ver.payloads.size_bytes)}</span
                    >
                  </div>
                  <div class="mt-1 font-mono font-semibold text-foreground">
                    {ver.payloads.exists ? `${formatNumber(ver.payloads.rows)} rows` : 'missing'}
                  </div>
                </div>
              </div>
            </div>
          {/each}
        </div>
      {/if}
    </div>
  </DialogContent>
</Dialog>
