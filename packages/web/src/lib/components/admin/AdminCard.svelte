<script lang="ts">
  import { Button } from '$lib/components/ui/button/index.js';
  import { Loader2 } from '@lucide/svelte';

  type FileChip = {
    name: string;
    size_bytes: number;
    rows?: number;
  };

  interface Props {
    title: string;
    subtitle: string;
    actionLabel: string;
    onAction: () => void;
    disabled?: boolean;
    busy?: boolean;
    statusMessage?: string;
    hasLogs?: boolean;
    onOpenLogs?: () => void;
    versionsCount?: number;
    onManageVersions?: () => void;
    primaryLabel: string;
    primaryValue: number;
    primaryMax: number;
    secondaryLabel?: string;
    secondaryValue?: number;
    secondaryMax?: number;
    files?: FileChip[];
  }

  let {
    title,
    subtitle,
    actionLabel,
    onAction,
    disabled = false,
    busy = false,
    statusMessage = '',
    hasLogs = false,
    onOpenLogs,
    versionsCount,
    onManageVersions,
    primaryLabel,
    primaryValue,
    primaryMax,
    secondaryLabel,
    secondaryValue = 0,
    secondaryMax = 1,
    files = [],
  }: Props = $props();

  const countFormatter = new Intl.NumberFormat('en-US');
  const radius = 28;
  const circumference = 2 * Math.PI * radius;

  function formatCount(value: number | null | undefined): string {
    return countFormatter.format(value || 0);
  }

  function formatBytes(value: number | null | undefined): string {
    const bytes = value || 0;
    if (bytes < 1024) return `${bytes} B`;
    if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
    if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
    return `${(bytes / (1024 * 1024 * 1024)).toFixed(1)} GB`;
  }

  function fileLabel(file: FileChip): string {
    return file.name.replace(/\.parquet$/, '').replaceAll('_', ' ');
  }

  function fileColor(name: string): string {
    if (
      name.startsWith('entities') ||
      name === 'gene_to_entrez' ||
      name === 'gene_to_uniprot' ||
      name === 'entrez_aliases'
    )
      return 'var(--chart-1)';
    if (
      name.startsWith('relations') ||
      name === 'chemical_to_inchikey' ||
      name === 'inchikey_aliases'
    )
      return 'var(--chart-2)';
    if (name.startsWith('evidence')) return 'var(--chart-3)';
    return 'var(--chart-4, var(--chart-3))';
  }

  function logWidth(value: number, max: number): number {
    if (value <= 0 || max <= 0) return 0;
    return Math.max(6, (Math.log10(value + 1) / Math.log10(max + 1)) * 100);
  }

  const totalBytes = $derived(files.reduce((sum, file) => sum + (file.size_bytes || 0), 0));
  const segments = $derived.by(() => {
    let offset = 0;
    const total = totalBytes || 1;
    return files.map((file) => {
      const fraction = (file.size_bytes || 0) / total;
      const length = fraction * circumference;
      const segment = {
        file,
        color: fileColor(file.name),
        dash: length,
        offset: -offset,
      };
      offset += length;
      return segment;
    });
  });
</script>

<article
  class={`flex flex-col gap-4 rounded-[1.25rem] border bg-card/70 p-4 ${
    busy ? 'border-primary/60' : 'border-border/50'
  }`}
>
  <div class="flex items-start justify-between gap-3">
    <div class="min-w-0">
      <div class="truncate text-lg font-semibold tracking-tight">{title}</div>
      <div class="mt-0.5 text-xs text-muted-foreground">{subtitle}</div>
    </div>
    <Button size="sm" onclick={onAction} {disabled}>
      {#if busy}
        <Loader2 class="size-3.5 animate-spin" />
      {/if}
      {actionLabel}
    </Button>
  </div>

  <div class="flex items-center gap-4">
    <svg viewBox="0 0 80 80" class="size-20 shrink-0" aria-hidden="true">
      <circle
        cx="40"
        cy="40"
        r={radius}
        fill="none"
        stroke="currentColor"
        class="text-muted/60"
        stroke-width="12"
      />
      {#each segments as segment (segment.file.name)}
        <circle
          cx="40"
          cy="40"
          r={radius}
          fill="none"
          stroke={segment.color}
          stroke-width="12"
          stroke-dasharray={`${segment.dash} ${circumference}`}
          stroke-dashoffset={segment.offset}
          transform="rotate(-90 40 40)"
        >
          <title>{fileLabel(segment.file)} · {formatBytes(segment.file.size_bytes)}</title>
        </circle>
      {/each}
    </svg>

    <div class="min-w-0 flex-1 space-y-3">
      <div>
        <div class="mb-1 flex items-baseline justify-between gap-2 text-xs">
          <span class="text-muted-foreground">{primaryLabel}</span>
          <span class="font-mono text-foreground">{formatCount(primaryValue)}</span>
        </div>
        <div class="h-2 overflow-hidden rounded-full bg-muted">
          <div
            class="h-full rounded-full bg-chart-1"
            style={`width: ${logWidth(primaryValue, primaryMax)}%`}
          ></div>
        </div>
      </div>
      {#if secondaryLabel}
        <div>
          <div class="mb-1 flex items-baseline justify-between gap-2 text-xs">
            <span class="text-muted-foreground">{secondaryLabel}</span>
            <span class="font-mono text-foreground">{formatCount(secondaryValue)}</span>
          </div>
          <div class="h-2 overflow-hidden rounded-full bg-muted">
            <div
              class="h-full rounded-full bg-chart-2"
              style={`width: ${logWidth(secondaryValue, secondaryMax)}%`}
            ></div>
          </div>
        </div>
      {/if}
    </div>
  </div>

  {#if files.length || hasLogs}
    <div class="flex flex-wrap items-center gap-2">
      {#each files as file (file.name)}
        <span
          class="inline-flex items-center gap-1.5 rounded-lg border border-border/60 bg-muted/25 px-2.5 py-1.5 text-xs text-muted-foreground"
        >
          <span class="size-2 rounded-full" style={`background: ${fileColor(file.name)}`}></span>
          <span class="font-mono">{fileLabel(file)}</span>
          <span class="text-[11px]">{formatBytes(file.size_bytes)}</span>
        </span>
      {/each}
      {#if versionsCount && onManageVersions}
        <button
          type="button"
          class="inline-flex items-center gap-1 rounded-md border border-border/80 bg-background/80 px-2 py-1 text-xs font-medium text-foreground transition-colors hover:bg-accent"
          onclick={onManageVersions}
        >
          <span>{versionsCount} version{versionsCount === 1 ? '' : 's'}</span>
        </button>
      {/if}
      {#if hasLogs && onOpenLogs}
        <button
          type="button"
          class="text-xs font-medium text-primary underline-offset-4 hover:underline"
          onclick={onOpenLogs}
        >
          Logs
        </button>
      {/if}
    </div>
  {/if}

  {#if busy && statusMessage}
    <p class="truncate text-xs text-muted-foreground">{statusMessage}</p>
  {/if}
</article>
