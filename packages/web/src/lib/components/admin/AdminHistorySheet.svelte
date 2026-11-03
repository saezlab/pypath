<script lang="ts">
  import {
    ArrowLeft,
    Check,
    CheckCircle2,
    ChevronRight,
    Clock,
    Copy,
    History,
    Layers,
    Loader2,
    RefreshCw,
    Terminal,
    XCircle,
  } from '@lucide/svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import {
    Sheet,
    SheetContent,
    SheetDescription,
    SheetHeader,
    SheetTitle,
  } from '$lib/components/ui/sheet/index.js';
  import JobPhases from '$lib/components/admin/JobPhases.svelte';
  import { createPoller } from '$lib/features/admin/polling';
  import { fetchAdminLog } from '$lib/features/admin/operations';
  import type { Job as AdminJob } from '$lib/features/admin/types';

  interface Props {
    open?: boolean;
    jobs?: AdminJob[];
    selectedJobId?: string | null;
    nowSeconds?: number;
    onRefresh?: () => void;
    onCancelJob?: (id: string) => void;
  }

  let {
    open = $bindable(false),
    jobs = [],
    selectedJobId = $bindable(null),
    nowSeconds = 0,
    onRefresh,
    onCancelJob,
  }: Props = $props();

  let logText = $state('');
  let logLoading = $state(false);
  let copied = $state(false);
  let copyTimer: ReturnType<typeof setTimeout> | null = null;

  const selectedJob = $derived(jobs.find((j) => j.id === selectedJobId) || null);

  async function loadJobLogs(jobId: string, signal: AbortSignal) {
    logLoading = true;
    try {
      const text = await fetchAdminLog('job', jobId, '', signal);
      if (!signal.aborted) logText = text;
    } catch {
      if (!signal.aborted) logText = (jobs.find((job) => job.id === jobId)?.logs || []).join('\n');
    } finally {
      if (!signal.aborted) logLoading = false;
    }
  }

  $effect(() => {
    const id = selectedJobId;
    const running = selectedJob?.status === 'running' || selectedJob?.status === 'queued';
    const polling = createPoller(
      (signal) => loadJobLogs(id!, signal),
      () => 1500,
    );
    if (id && open) {
      if (running) polling.start();
      else void polling.refresh();
    } else logText = '';
    return polling.stop;
  });

  async function copyLogs() {
    try {
      await navigator.clipboard.writeText(logText || '');
      copied = true;
      if (copyTimer) clearTimeout(copyTimer);
      copyTimer = setTimeout(() => {
        copied = false;
      }, 1600);
    } catch {
      copied = false;
    }
  }

  function formatTime(timestamp: number | null | undefined): string {
    if (!timestamp) return '';
    const date = new Date(timestamp * 1000);
    return (
      date.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' }) +
      ' · ' +
      date.toLocaleDateString([], { month: 'short', day: 'numeric' })
    );
  }

  function formatJobTitle(job: AdminJob): string {
    if (job.action === 'build_resources' || job.action === 'build_resource') {
      const sources =
        (job.params.sources as string[]) || (job.params.source ? [String(job.params.source)] : []);
      if (sources.length) return `Build ${sources.join(', ')}`;
      return 'Build resources';
    }
    if (job.action === 'export_hubs') {
      const hubs = job.params.hubs as string[] | undefined;
      if (hubs?.length) return `Export hubs (${hubs.join(', ')})`;
      return 'Export all hubs';
    }
    return job.action.replaceAll('_', ' ');
  }

  function formatDuration(seconds: number | null | undefined): string {
    if (seconds == null || Number.isNaN(seconds)) return '0s';
    if (seconds < 10) return `${seconds.toFixed(1)}s`;
    if (seconds < 60) return `${Math.round(seconds)}s`;
    const mins = Math.floor(seconds / 60);
    const secs = Math.round(seconds % 60);
    return `${mins}m ${String(secs).padStart(2, '0')}s`;
  }

  function statusColorClasses(status: string): string {
    if (status === 'done') return 'bg-emerald-500/10 text-emerald-600';
    if (status === 'failed') return 'bg-destructive/10 text-destructive';
    if (status === 'running') return 'bg-primary/10 text-primary';
    return 'bg-muted text-muted-foreground';
  }

  function statusBadgeClasses(status: string): string {
    if (status === 'done') return 'border-emerald-500/30 text-emerald-600';
    if (status === 'failed') return 'border-destructive/30 text-destructive';
    if (status === 'running') return 'border-primary/30 text-primary';
    return 'border-border text-muted-foreground';
  }
</script>

<Sheet bind:open>
  <SheetContent
    side="right"
    class="flex w-[95vw] flex-col gap-0 overflow-hidden p-0 sm:max-w-2xl lg:max-w-3xl"
  >
    <SheetHeader class="flex flex-row items-center justify-between border-b px-6 py-4 pr-12">
      {#if selectedJob}
        <div class="flex items-center gap-2">
          <Button variant="ghost" size="sm" class="h-8 px-2" onclick={() => (selectedJobId = null)}>
            <ArrowLeft class="mr-1 size-4" />
            Back
          </Button>
          <div class="h-4 w-px bg-border"></div>
          <div class="min-w-0">
            <SheetTitle class="truncate text-base font-semibold leading-none">
              {formatJobTitle(selectedJob)}
            </SheetTitle>
            <SheetDescription class="mt-0.5 truncate font-mono text-xs text-muted-foreground">
              #{selectedJob.id} · {formatTime(selectedJob.created_at)}
            </SheetDescription>
          </div>
        </div>
      {:else}
        <div>
          <div class="flex items-center gap-2">
            <History class="size-5 text-primary" />
            <SheetTitle class="text-base font-semibold leading-none">Job History</SheetTitle>
            {#if jobs.length}
              <Badge variant="secondary" class="rounded-full text-xs">{jobs.length}</Badge>
            {/if}
          </div>
          <SheetDescription class="mt-0.5 text-xs text-muted-foreground">
            Previous build & export runs
          </SheetDescription>
        </div>
      {/if}

      {#if !selectedJob && onRefresh}
        <Button variant="ghost" size="icon" class="size-8" onclick={onRefresh}>
          <RefreshCw class="size-3.5" />
        </Button>
      {/if}
    </SheetHeader>

    <div class="min-h-0 flex-1 overflow-y-auto">
      {#if selectedJob}
        <div class="flex flex-col gap-6 p-6">
          <div
            class="flex flex-wrap items-center justify-between gap-3 rounded-xl border bg-card/60 p-4"
          >
            <div class="flex items-center gap-3">
              <span
                class={`flex size-8 shrink-0 items-center justify-center rounded-full ${statusColorClasses(selectedJob.status)}`}
              >
                {#if selectedJob.status === 'running'}
                  <Loader2 class="size-4 animate-spin" />
                {:else if selectedJob.status === 'done'}
                  <CheckCircle2 class="size-4" />
                {:else if selectedJob.status === 'failed'}
                  <XCircle class="size-4" />
                {:else}
                  <Clock class="size-4" />
                {/if}
              </span>
              <div>
                <div class="flex items-center gap-2">
                  <span class="text-sm font-semibold capitalize">{selectedJob.status}</span>
                  <span class="text-xs text-muted-foreground">·</span>
                  <span class="font-mono text-xs text-muted-foreground">
                    {formatDuration(selectedJob.elapsed_seconds)}
                  </span>
                </div>
                <div class="mt-0.5 text-xs text-muted-foreground">
                  Action: <span class="font-mono">{selectedJob.action}</span>
                  {#if selectedJob.params.max_records}
                    · cap: <span class="font-mono">{selectedJob.params.max_records} rows</span>
                  {/if}
                  {#if typeof selectedJob.params.parallel === 'number' && selectedJob.params.parallel > 1}
                    · <span class="font-mono">{selectedJob.params.parallel} workers</span>
                  {/if}
                </div>
              </div>
            </div>

            {#if (selectedJob.status === 'running' || selectedJob.status === 'queued') && onCancelJob}
              <Button variant="destructive" size="sm" onclick={() => onCancelJob?.(selectedJob.id)}>
                Cancel
              </Button>
            {/if}
          </div>

          {#if selectedJob.error}
            <div
              class="rounded-xl border border-destructive/40 bg-destructive/10 p-3 text-xs text-destructive"
            >
              <strong>Error:</strong>
              {selectedJob.error}
            </div>
          {/if}

          <div class="flex flex-col gap-2.5">
            <div
              class="flex items-center gap-1.5 text-xs font-semibold uppercase tracking-wider text-muted-foreground"
            >
              <Layers class="size-3.5" />
              <span>Stage Summary ({selectedJob.stages.length})</span>
            </div>
            <div class="rounded-xl border bg-card/60 p-3">
              {#if selectedJob.stages.length > 0}
                <JobPhases stages={selectedJob.stages} {nowSeconds} />
              {:else}
                <p class="py-2 text-xs text-muted-foreground">
                  No individual stage records available.
                </p>
              {/if}
            </div>
          </div>

          <div class="flex flex-col gap-2.5">
            <div class="flex items-center justify-between gap-2">
              <div
                class="flex items-center gap-1.5 text-xs font-semibold uppercase tracking-wider text-muted-foreground"
              >
                <Terminal class="size-3.5" />
                <span>Logs</span>
              </div>
              <Button
                variant="outline"
                size="sm"
                class="h-7 gap-1 px-2.5 text-xs"
                onclick={copyLogs}
                disabled={logLoading || !logText}
              >
                {#if copied}
                  <Check class="size-3" />
                  <span>Copied</span>
                {:else}
                  <Copy class="size-3" />
                  <span>Copy Logs</span>
                {/if}
              </Button>
            </div>
            <div
              class="max-h-[450px] min-h-[220px] overflow-auto rounded-xl border border-zinc-800 bg-zinc-950 p-4 font-mono text-xs leading-relaxed text-zinc-100"
            >
              {#if logLoading}
                <div class="flex items-center justify-center gap-2 py-6 text-zinc-400">
                  <Loader2 class="size-4 animate-spin" />
                  <span>Loading logs…</span>
                </div>
              {:else if logText}
                <pre class="whitespace-pre-wrap break-words">{logText}</pre>
              {:else}
                <p class="py-4 text-center text-zinc-500">No logs saved for this job.</p>
              {/if}
            </div>
          </div>
        </div>
      {:else if jobs.length === 0}
        <div
          class="flex flex-col items-center justify-center p-12 text-center text-muted-foreground"
        >
          <History class="mb-2 size-8 stroke-1 text-muted-foreground/50" />
          <p class="text-sm font-medium">No jobs in history yet</p>
          <p class="mt-1 text-xs text-muted-foreground">
            Rebuild resources or export hubs to see build logs and stage progress here.
          </p>
        </div>
      {:else}
        <div class="divide-y divide-border/60">
          {#each jobs as job (job.id)}
            <button
              type="button"
              class="group flex w-full items-center justify-between gap-4 p-4 text-left transition-colors hover:bg-muted/40"
              onclick={() => (selectedJobId = job.id)}
            >
              <div class="flex min-w-0 items-start gap-3">
                <span
                  class={`mt-0.5 flex size-7 shrink-0 items-center justify-center rounded-full ${statusColorClasses(job.status)}`}
                >
                  {#if job.status === 'running'}
                    <Loader2 class="size-3.5 animate-spin" />
                  {:else if job.status === 'done'}
                    <CheckCircle2 class="size-3.5" />
                  {:else if job.status === 'failed'}
                    <XCircle class="size-3.5" />
                  {:else}
                    <Clock class="size-3.5" />
                  {/if}
                </span>
                <div class="min-w-0">
                  <div class="flex items-center gap-2">
                    <span
                      class="truncate text-sm font-semibold transition-colors group-hover:text-primary"
                    >
                      {formatJobTitle(job)}
                    </span>
                    <span
                      class={`rounded-md border px-1.5 py-0.2 font-mono text-[10px] capitalize ${statusBadgeClasses(job.status)}`}
                    >
                      {job.status}
                    </span>
                  </div>
                  <div
                    class="mt-1 flex flex-wrap items-center gap-x-2 gap-y-0.5 text-xs text-muted-foreground"
                  >
                    <span>{formatTime(job.created_at)}</span>
                    <span>·</span>
                    <span class="font-mono">{formatDuration(job.elapsed_seconds)}</span>
                    {#if job.stages.length}
                      <span>·</span>
                      <span>{job.stages.length} stages</span>
                    {/if}
                  </div>
                  {#if job.error}
                    <p class="mt-1 truncate text-xs text-destructive">{job.error}</p>
                  {/if}
                </div>
              </div>
              <ChevronRight
                class="size-4 shrink-0 text-muted-foreground/50 transition-all group-hover:translate-x-0.5 group-hover:text-foreground"
              />
            </button>
          {/each}
        </div>
      {/if}
    </div>
  </SheetContent>
</Sheet>
