<script lang="ts">
  import { onMount } from 'svelte';
  import * as Tooltip from '$lib/components/ui/tooltip/index.js';
  import { createPoller } from '$lib/features/admin/polling';
  import type { components } from '$lib/api/types.generated';

  type Status = components['schemas']['ServingStatus'];

  /** How busy the service is, refreshed every 15 s while the page is visible. */
  let status = $state<Status | null>(null);
  let unreachable = $state(false);

  const poller = createPoller(
    async (signal) => {
      const response = await fetch('/app-api/status', { signal });
      if (!response.ok) throw new Error(`Status ${response.status}`);
      status = (await response.json()) as Status;
      unreachable = false;
    },
    () => 15_000,
    () => (unreachable = true),
  );

  onMount(() => {
    const sync = () => (document.visibilityState === 'visible' ? poller.start() : poller.stop());
    // A background tab reads the status once, so it never shows only "…"; it polls once seen.
    if (document.visibilityState === 'visible') poller.start();
    else void poller.refresh();
    document.addEventListener('visibilitychange', sync);
    return () => {
      document.removeEventListener('visibilitychange', sync);
      poller.stop();
    };
  });

  const level = $derived<Status['state'] | 'unreachable' | 'unknown'>(
    unreachable ? 'unreachable' : (status?.state ?? 'unknown'),
  );
  const label = $derived(
    {
      ok: 'OK',
      busy: 'Busy',
      queued: 'Queued',
      unreachable: 'Offline',
      unknown: '…',
    }[level],
  );
  const percent = (value: number | null | undefined) =>
    value == null ? 'unknown' : `${Math.round(value * 100)}%`;
</script>

<Tooltip.Provider delayDuration={150}>
  <Tooltip.Root>
    <Tooltip.Trigger
      class="inline-flex h-8 items-center gap-1.5 rounded-md px-2 text-xs text-muted-foreground hover:text-foreground focus-visible:ring-2 focus-visible:ring-ring focus-visible:outline-none"
      aria-label={`Service status: ${label}`}
    >
      <span class={`status-dot ${level}`} aria-hidden="true"></span>
      <span class="hidden sm:inline">{label}</span>
    </Tooltip.Trigger>
    <Tooltip.Content side="bottom" sideOffset={6} class="max-w-xs text-xs">
      {#if unreachable}
        The service did not answer; results may not load.
      {:else if status}
        <div class="space-y-0.5 tabular-nums">
          <p class="font-medium">
            {level === 'queued'
              ? 'Requests are waiting for a free query slot'
              : level === 'busy'
                ? 'The service is busy; responses may be slower'
                : 'The service is responsive'}
          </p>
          <p>
            Queries: {status.queries.running} of {status.queries.slots} slots in use, {status
              .queries.waiting} waiting
          </p>
          <p>
            Last minute: {status.requestsLastMinute} requests{status.medianResponseMs != null
              ? `, median ${Math.round(status.medianResponseMs)} ms`
              : ''}
          </p>
          <p>CPU load {percent(status.cpuLoad)} · memory {percent(status.memory)}</p>
        </div>
      {:else}
        Checking the service…
      {/if}
    </Tooltip.Content>
  </Tooltip.Root>
</Tooltip.Provider>

<style>
  .status-dot {
    width: 0.5rem;
    height: 0.5rem;
    border-radius: 9999px;
    background: var(--muted-foreground);
    opacity: 0.5;
  }
  .status-dot.ok {
    background: #22c55e;
    opacity: 1;
  }
  .status-dot.busy {
    background: #f59e0b;
    opacity: 1;
  }
  .status-dot.queued,
  .status-dot.unreachable {
    background: #ef4444;
    opacity: 1;
  }
</style>
