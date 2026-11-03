<script lang="ts">
  import { Check, Circle, Loader2, X } from '@lucide/svelte';

  import type { Stage } from '$lib/features/admin/types';

  interface Props {
    stages: Stage[];
    nowSeconds?: number;
  }

  let { stages, nowSeconds = 0 }: Props = $props();

  const activeIndex = $derived.by(() => {
    let index = -1;
    stages.forEach((stage, i) => {
      if (stage.status === 'running') index = i;
    });
    return index;
  });

  const verbs: Record<string, string> = {
    emit: 'Export',
    ingest: 'Ingest',
    build: 'Build',
    resolve: 'Resolve',
  };

  function stageLabel(id: string): string {
    const split = id.indexOf(':');
    if (split < 0) {
      return id.replaceAll('_', ' ');
    }
    const kind = id.slice(0, split);
    const rest = id
      .slice(split + 1)
      .replaceAll('_', ' ')
      .replaceAll('.', ' · ');
    const verb = verbs[kind];
    return verb ? `${verb} ${rest}` : `${kind} ${rest}`;
  }

  function formatDuration(seconds: number | null | undefined): string {
    if (seconds == null || Number.isNaN(seconds)) return '';
    if (seconds < 10) return `${seconds.toFixed(1)}s`;
    if (seconds < 60) return `${Math.round(seconds)}s`;
    const minutes = Math.floor(seconds / 60);
    const remainder = Math.round(seconds % 60);
    return `${minutes}m ${String(remainder).padStart(2, '0')}s`;
  }

  function stageSeconds(stage: Stage): number | null {
    if (stage.status === 'running' && stage.started_at) {
      const clock = nowSeconds || Date.now() / 1000;
      return Math.max(0, clock - stage.started_at);
    }
    return stage.elapsed_seconds ?? null;
  }

  function countsLabel(stage: Stage): string {
    if (stage.current == null) return '';
    if (stage.total == null) return String(stage.current);
    return `${stage.current} / ${stage.total}`;
  }
</script>

<ol class="max-h-64 space-y-1 overflow-y-auto pr-1">
  {#each stages as stage, index (stage.id)}
    {@const active = index === activeIndex}
    {@const elapsed = formatDuration(stageSeconds(stage))}
    <li
      class={`flex items-start gap-2 rounded-md px-1.5 py-1 text-xs ${
        active ? 'bg-primary/10 text-foreground' : 'text-muted-foreground'
      }`}
    >
      <span class="mt-0.5 flex size-4 shrink-0 items-center justify-center">
        {#if active}
          <Loader2 class="size-3.5 animate-spin text-primary" />
        {:else if stage.status === 'done'}
          <Check class="size-3.5 text-primary" />
        {:else if stage.status === 'failed'}
          <X class="size-3.5 text-destructive" />
        {:else if stage.status === 'cancelled'}
          <X class="size-3.5 text-muted-foreground" />
        {:else}
          <Circle class="size-2.5 text-muted-foreground/70" />
        {/if}
      </span>
      <div class="min-w-0 flex-1">
        <div class="flex items-baseline justify-between gap-2">
          <span class={`truncate ${active || stage.status === 'done' ? 'text-foreground' : ''}`}>
            {stageLabel(stage.id)}
          </span>
          {#if elapsed}
            <span class="shrink-0 font-mono tabular-nums text-[11px]">{elapsed}</span>
          {/if}
        </div>
        {#if active && (stage.message || countsLabel(stage))}
          <p class="truncate text-[11px] text-muted-foreground">
            {stage.message}{countsLabel(stage) ? ` · ${countsLabel(stage)}` : ''}
          </p>
        {/if}
      </div>
    </li>
  {/each}
</ol>
