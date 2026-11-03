<script lang="ts">
  import { browser } from '$app/environment';
  import { goto } from '$app/navigation';
  import { page } from '$app/stores';
  import { onMount } from 'svelte';
  import { compareResourceVersions } from '$lib/utils/versions';
  import AdminCard from '$lib/components/admin/AdminCard.svelte';
  import AdminHistorySheet from '$lib/components/admin/AdminHistorySheet.svelte';
  import JobPhases from '$lib/components/admin/JobPhases.svelte';
  import LogDialog from '$lib/components/admin/LogDialog.svelte';
  import ResourceVersionsDialog from '$lib/components/admin/ResourceVersionsDialog.svelte';
  import ExploreBrowserShell from '$lib/components/explore/ExploreBrowserShell.svelte';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Input } from '$lib/components/ui/input/index.js';
  import { History } from '@lucide/svelte';
  import { clearAdminPassword, setAdminPassword } from '$lib/admin/auth';

  import type {
    AdminStatus,
    AvailableSource,
    Job,
    BuiltResource,
    ResolutionStats,
    FileInfo,
  } from '$lib/features/admin/types';
  import {
    fetchAdminStatus,
    fetchAdminSources,
    submitAdminJob,
    cancelAdminJob,
    fetchAdminLog,
    removeResourceVersion,
    AdminRequestError,
  } from '$lib/features/admin/operations';
  import { createPoller } from '$lib/features/admin/polling';

  let status = $state<AdminStatus | null>(null);
  let sources = $state<AvailableSource[]>([]);
  let error = $state('');
  let busy = $state(false);
  let maxRecords = $state('');
  let parallelJobs = $state('1');
  let resourceVersion = $state('');
  let inputRef = $state<HTMLInputElement | null>(null);
  let draftQuery = $state('');
  let nowSeconds = $state(Date.now() / 1000);
  let authChecked = $state(false);
  let locked = $state(false);
  let password = $state('');
  let authError = $state('');
  let logOpen = $state(false);
  let logTitle = $state('Build log');
  let logText = $state('');
  let logLoading = $state(false);
  let logKind = $state('');
  let logName = $state('');
  let logVersion = $state('');

  let activeJobIds = $state<string[]>([]);
  let dismissedJobId = $state<string | null>(null);
  let historyOpen = $state(false);
  let historySelectedJobId = $state<string | null>(null);

  let versionsOpen = $state(false);
  let selectedResourceName = $state('');

  const rawTab = $derived($page.url.searchParams.get('tab') || 'resolver');
  const tab = $derived(rawTab === 'resources' ? 'resources' : 'resolver');
  const query = $derived($page.url.searchParams.get('q') || '');
  const jobRunning = $derived(
    status?.job?.status === 'running' || status?.job?.status === 'queued',
  );
  const canBuild = $derived(!busy && !jobRunning);
  const allJobs = $derived<Job[]>(status?.jobs ?? (status?.job ? [status.job] : []));
  const showJobBanner = $derived.by(() => {
    const job = status?.job;
    if (!job || job.id === dismissedJobId) return false;
    return jobRunning || activeJobIds.includes(job.id);
  });

  const runningStage = $derived(
    status?.job?.stages.find((stage) => stage.status === 'running') ?? null,
  );

  type ResolverCard = {
    id: string;
    title: string;
    file: FileInfo;
  };

  type ResourceCardModel = {
    id: string;
    title: string;
    built: BuiltResource | null;
    versions: BuiltResource[];
    datasets: number;
  };

  const resolverCards = $derived.by((): ResolverCard[] =>
    (status?.resolver.hubs ?? []).map((file) => ({ id: file.name, title: file.name, file })),
  );

  function builtStamp(item: BuiltResource): number {
    return Math.max(item.entities.mtime || 0, item.relations.mtime || 0, item.payloads.mtime || 0);
  }

  function builtUsable(item: BuiltResource): boolean {
    return Boolean(item.entities.exists && item.relations.exists);
  }

  function resolutionSummary(stats: ResolutionStats | null | undefined): string {
    if (!stats) return '';
    const lookupEntities = Math.max(
      0,
      stats.lookup_entities ?? stats.input_entities - stats.not_applicable_entities,
    );
    const structure = stats.structure_entities
      ? ` · ${stats.structure_entities.toLocaleString()} by structure`
      : '';
    return `${stats.resolved_entities.toLocaleString()} resolved / ${lookupEntities.toLocaleString()} lookup${structure} · ${stats.unresolved_entities.toLocaleString()} unresolved`;
  }

  function resourceSubtitle(built: BuiltResource | null, versionCount: number): string {
    if (!built) return 'Not built';
    const size = fmtBytes(
      built.entities.size_bytes + built.relations.size_bytes + built.payloads.size_bytes,
    );
    const version = versionCount > 1 ? `${versionCount} versions` : built.version;
    const resolution = resolutionSummary(built.resolution_stats);
    return `${size} parquet · ${version}${resolution ? ` · ${resolution}` : ''}`;
  }

  function preferBuilt(next: BuiltResource, current: BuiltResource | undefined): boolean {
    if (!current) return true;
    const nextOk = builtUsable(next);
    const currentOk = builtUsable(current);
    if (nextOk !== currentOk) return nextOk;
    return (
      compareResourceVersions(
        next.version,
        current.version,
        builtStamp(next),
        builtStamp(current),
      ) > 0
    );
  }

  const resourceCards = $derived.by((): ResourceCardModel[] => {
    const builtByName = new Map<string, BuiltResource>();
    const allVersionsByName = new Map<string, BuiltResource[]>();
    for (const item of status?.resources.built ?? []) {
      const list = allVersionsByName.get(item.resource) ?? [];
      list.push(item);
      allVersionsByName.set(item.resource, list);
      const current = builtByName.get(item.resource);
      if (preferBuilt(item, current)) builtByName.set(item.resource, item);
    }
    const names = [
      ...sources.map((item) => item.source),
      ...[...allVersionsByName.keys()].filter(
        (name) => !sources.some((item) => item.source === name),
      ),
    ];
    return names.map((name) => ({
      id: name,
      title: name,
      built: builtByName.get(name) ?? null,
      versions: allVersionsByName.get(name) ?? [],
      datasets: sources.find((item) => item.source === name)?.datasets.length ?? 0,
    }));
  });

  const selectedResourceVersions = $derived.by((): BuiltResource[] => {
    if (!selectedResourceName) return [];
    const card = resourceCards.find((c) => c.id === selectedResourceName);
    return card?.versions ?? [];
  });

  const selectedResourceActiveVersion = $derived.by((): string | null => {
    if (!selectedResourceName) return null;
    const card = resourceCards.find((c) => c.id === selectedResourceName);
    return card?.built?.version ?? null;
  });

  const searchedResolver = $derived.by(() => {
    const needle = query.trim().toLowerCase();
    if (!needle) return resolverCards;
    return resolverCards.filter(
      (card) => card.title.toLowerCase().includes(needle) || card.id.toLowerCase().includes(needle),
    );
  });

  const searchedResources = $derived.by(() => {
    const needle = query.trim().toLowerCase();
    if (!needle) return resourceCards;
    return resourceCards.filter((card) => card.title.toLowerCase().includes(needle));
  });

  const maxHubRows = $derived(Math.max(1, ...resolverCards.map((card) => card.file.rows || 0)));
  const maxEntities = $derived(
    Math.max(1, ...resourceCards.map((card) => card.built?.entities.rows || 0)),
  );
  const maxRelations = $derived(
    Math.max(1, ...resourceCards.map((card) => card.built?.relations.rows || 0)),
  );

  $effect(() => {
    draftQuery = query;
  });

  function setTab(next: string) {
    const url = new URL($page.url);
    url.searchParams.set('tab', next);
    goto(url, { replaceState: true, keepFocus: true, noScroll: true });
  }

  function setQuery(next: string) {
    const url = new URL($page.url);
    if (next.trim()) url.searchParams.set('q', next.trim());
    else url.searchParams.delete('q');
    goto(url, { replaceState: true, keepFocus: true, noScroll: true });
  }

  function submitSearch() {
    setQuery(draftQuery);
  }

  function fmtBytes(value: number) {
    if (!value) return '0 B';
    if (value < 1024) return `${value} B`;
    if (value < 1024 ** 2) return `${(value / 1024).toFixed(1)} KB`;
    if (value < 1024 ** 3) return `${(value / 1024 ** 2).toFixed(1)} MB`;
    return `${(value / 1024 ** 3).toFixed(1)} GB`;
  }

  function parsedMaxRecords() {
    const trimmed = maxRecords.trim();
    if (!trimmed) return null;
    const value = Number(trimmed);
    return Number.isFinite(value) && value >= 0 ? value : null;
  }

  function parsedParallel() {
    const trimmed = parallelJobs.trim();
    if (!trimmed) return 1;
    const value = Number(trimmed);
    return Number.isFinite(value) && value >= 1 ? Math.min(32, Math.floor(value)) : 1;
  }

  let refreshController: AbortController | null = null;
  async function refresh(parentSignal?: AbortSignal) {
    refreshController?.abort();
    refreshController = new AbortController();
    const signal = parentSignal
      ? AbortSignal.any([parentSignal, refreshController.signal])
      : refreshController.signal;
    try {
      const next = await fetchAdminStatus(signal);
      if (signal?.aborted) return false;
      authChecked = true;
      status = next;
      if (status?.job && (status.job.status === 'running' || status.job.status === 'queued')) {
        if (!activeJobIds.includes(status.job.id)) {
          activeJobIds = [...activeJobIds, status.job.id];
        }
      }
      error = '';
      locked = false;
      authError = '';
      return true;
    } catch (err) {
      if (signal?.aborted) return false;
      if (err instanceof AdminRequestError && err.status === 401) {
        authChecked = true;
        locked = true;
        status = null;
        return false;
      }
      authChecked = true;
      error = err instanceof Error ? err.message : String(err);
      return false;
    }
  }

  async function loadSources() {
    try {
      sources = await fetchAdminSources();
    } catch (err) {
      if (err instanceof AdminRequestError && err.status === 401) locked = true;
      else error = err instanceof Error ? err.message : String(err);
    }
  }

  async function trigger(body: Record<string, unknown>) {
    busy = true;
    error = '';
    try {
      await submitAdminJob(body);
      await refresh();
      if (status?.job?.id && !activeJobIds.includes(status.job.id)) {
        activeJobIds = [...activeJobIds, status.job.id];
      }
    } catch (err) {
      if (err instanceof AdminRequestError && err.status === 401) locked = true;
      error = err instanceof Error ? err.message : String(err);
    } finally {
      busy = false;
    }
  }

  async function cancelJobById(id: string) {
    if (!id) return;
    await cancelAdminJob(id);
    await refresh();
  }

  async function cancelJob() {
    const id = status?.job?.id;
    if (!id) return;
    await cancelJobById(id);
  }

  let logController: AbortController | null = null;
  async function loadLogs(signal?: AbortSignal) {
    if (!logKind || !logName) return;
    logController?.abort();
    const controller = new AbortController();
    logController = controller;
    const requestSignal = signal ? AbortSignal.any([signal, controller.signal]) : controller.signal;
    logLoading = true;
    try {
      const text = await fetchAdminLog(logKind, logName, logVersion, requestSignal);
      if (!requestSignal.aborted) logText = text;
    } catch (err) {
      if (!requestSignal.aborted) logText = err instanceof Error ? err.message : String(err);
    } finally {
      if (!requestSignal.aborted) logLoading = false;
    }
  }

  function openLogs(kind: string, name: string, title: string, version: string = '') {
    logKind = kind;
    logName = name;
    logVersion = version;
    logTitle = title;
    logOpen = true;
    void loadLogs();
  }

  async function deleteResourceVersion(resource: string, version: string) {
    await removeResourceVersion(resource, version);
    await refresh();
  }

  async function cleanupOlderVersions(resource: string) {
    const card = resourceCards.find((c) => c.id === resource);
    if (!card) return;
    const activeVer = card.built?.version;
    const toDelete = card.versions.filter((v) => v.version !== activeVer);
    for (const item of toDelete) {
      await deleteResourceVersion(resource, item.version);
    }
  }

  function rebuildHub(name: string) {
    return trigger({
      action: 'export_hubs',
      hubs: [name],
      max_records: parsedMaxRecords(),
      parallel: parsedParallel(),
    });
  }

  function validateResourceVersion() {
    if (/^[0-9]+(?:\.[0-9]+)*$/.test(resourceVersion.trim())) return true;
    error = 'Enter a new resource version, for example 1.0.0';
    return false;
  }

  function rebuildResource(name: string) {
    if (!validateResourceVersion()) return;
    return trigger({
      action: 'build_resources',
      version: resourceVersion.trim(),
      sources: [name],
      max_records: parsedMaxRecords(),
      parallel: parsedParallel(),
    });
  }

  function buildAll() {
    if (tab === 'resolver') {
      return trigger({
        action: 'export_hubs',
        max_records: parsedMaxRecords(),
        parallel: parsedParallel(),
      });
    }
    if (!validateResourceVersion()) return;
    const names = resourceCards.map((card) => card.id);
    if (!names.length) {
      error = 'No resources to build';
      return;
    }
    return trigger({
      action: 'build_resources',
      version: resourceVersion.trim(),
      sources: names,
      max_records: parsedMaxRecords(),
      parallel: parsedParallel(),
    });
  }

  function jobTouchesHub(name: string) {
    const job = status?.job;
    if (!job || !jobRunning) return false;
    if (job.action !== 'export_hubs') return false;
    const hubs = job.params.hubs as string[] | undefined;
    if (!hubs?.length) return true;
    return hubs.includes(name);
  }

  function jobTouchesResource(name: string) {
    const job = status?.job;
    if (!job || !jobRunning) return false;
    if (job.action !== 'build_resource' && job.action !== 'build_resources') return false;
    const names =
      (job.params.sources as string[] | undefined) ||
      (job.params.source ? [String(job.params.source)] : []);
    if (!names.length) return true;
    return (
      names.includes(name) ||
      runningStage?.id?.includes(`:${name}`) ||
      runningStage?.id === `build:${name}`
    );
  }

  async function signIn(event: Event) {
    event.preventDefault();
    authError = '';
    setAdminPassword(password);
    const ok = await refresh();
    if (ok) {
      password = '';
      void loadSources();
      return;
    }
    if (locked) {
      clearAdminPassword();
      authError = 'Wrong password';
    }
  }

  const statusPolling = createPoller(
    async (signal) => {
      if (locked) return;
      const first = !authChecked;
      if ((await refresh(signal)) && first && !signal.aborted) await loadSources();
    },
    () => (jobRunning ? 500 : 4000),
  );
  const logPolling = createPoller(loadLogs, () => 1000);
  onMount(() => {
    statusPolling.start();
    return () => {
      statusPolling.stop();
      logPolling.stop();
      logController?.abort();
      refreshController?.abort();
    };
  });
  $effect(() => {
    if (!browser || !logOpen || !jobRunning) return;
    logPolling.start();
    return logPolling.stop;
  });

  $effect(() => {
    if (!browser || !jobRunning) return;
    nowSeconds = Date.now() / 1000;
    const timer = setInterval(() => {
      nowSeconds = Date.now() / 1000;
    }, 200);
    return () => clearInterval(timer);
  });

  $effect(() => {
    if (!browser) return;
    const onKeyDown = (event: KeyboardEvent) => {
      const target = event.target as HTMLElement | null;
      const tagName = target?.tagName;
      const isTypingTarget =
        tagName === 'INPUT' ||
        tagName === 'TEXTAREA' ||
        tagName === 'SELECT' ||
        target?.isContentEditable;
      if (event.key === '/' && !isTypingTarget) {
        event.preventDefault();
        inputRef?.focus();
        inputRef?.select();
      }
    };
    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  });
</script>

<svelte:head>
  <title>OmniPath Admin</title>
</svelte:head>

{#if !authChecked}
  <div class="py-16 text-center text-sm text-muted-foreground">Loading…</div>
{:else if locked}
  <div class="flex justify-center py-16">
    <form
      class="w-full max-w-sm space-y-4 rounded-[1.25rem] border border-border/50 bg-card/70 p-6"
      onsubmit={signIn}
    >
      <div>
        <h1 class="text-lg font-semibold tracking-tight">Admin</h1>
        <p class="mt-1 text-sm text-muted-foreground">Enter the admin password to continue.</p>
      </div>
      <Input
        type="password"
        bind:value={password}
        placeholder="Password"
        autocomplete="current-password"
        autofocus
      />
      {#if authError}
        <p class="text-sm text-destructive">{authError}</p>
      {/if}
      <Button type="submit" class="w-full">Sign in</Button>
    </form>
  </div>
{:else}
  <ExploreBrowserShell
    {query}
    {draftQuery}
    onDraftQueryChange={(value) => (draftQuery = value)}
    onSubmitSearch={submitSearch}
    {tab}
    onTabChange={setTab}
    tabs={[
      { value: 'resolver', label: 'resolver', badge: resolverCards.length },
      { value: 'resources', label: 'resources', badge: resourceCards.length },
    ]}
    searchPlaceholder={tab === 'resolver' ? 'Search hubs…' : 'Search resources…'}
    explanationTitle={tab === 'resolver' ? 'Build identifier lookups' : 'Build parquet resources'}
    explanationText={tab === 'resolver'
      ? 'Export identifier hubs, then combine them into unique lookups plus nested alias tables (target + matching ids). Rebuild one hub or run Build all.'
      : status?.pipeline_available === false
        ? 'Resource rebuilds need pypath in the API process. Hub rebuilds on the Resolver tab work without it.'
        : 'Each resource is three nested parquet files. Rebuild one source or run Build all. Empty row cap means uncapped.'}
    bind:searchInputRef={inputRef}
  >
    {#snippet summarySlot()}
      <div class="flex flex-col gap-3">
        <div class="flex flex-wrap items-center justify-between gap-3">
          <div class="flex flex-wrap items-center gap-3">
            <div class="flex items-center gap-1.5">
              {#if tab === 'resources'}
                <label class="flex items-center gap-2 text-xs text-muted-foreground"
                  >Resource version
                  <Input class="w-28" bind:value={resourceVersion} placeholder="1.0.0" />
                </label>
              {/if}
              <Input
                class="w-32"
                bind:value={maxRecords}
                placeholder="row cap"
                inputmode="numeric"
              />
              <span class="text-xs text-muted-foreground">cap</span>
            </div>
            <div class="flex items-center gap-1.5">
              <Input
                class="w-20"
                bind:value={parallelJobs}
                placeholder="1"
                inputmode="numeric"
                type="number"
                min="1"
                max="32"
              />
              <span class="text-xs text-muted-foreground">workers</span>
            </div>
          </div>
          <div class="flex flex-wrap items-center gap-2">
            <Button
              variant="outline"
              size="sm"
              onclick={() => {
                historySelectedJobId = null;
                historyOpen = true;
              }}
            >
              <History class="mr-1.5 size-3.5" />
              History
              {#if allJobs.length > 0}
                <span
                  class="ml-1.5 rounded-full bg-muted px-1.5 py-0.2 text-[10px] font-medium text-foreground"
                >
                  {allJobs.length}
                </span>
              {/if}
            </Button>
            <Button size="sm" onclick={buildAll} disabled={!canBuild}>Build all</Button>
          </div>
        </div>
        {#if error}
          <p
            class="rounded-md border border-destructive/40 bg-destructive/10 px-3 py-2 text-sm text-destructive"
          >
            {error}
          </p>
        {/if}
        {#if showJobBanner && status?.job}
          {@const job = status.job}
          {@const parallelCount =
            typeof job.params?.parallel === 'number' ? job.params.parallel : 1}
          <div
            class="flex flex-col gap-3 rounded-[1.25rem] border border-border/50 bg-card/70 px-4 py-3"
          >
            <div class="flex flex-wrap items-center justify-between gap-3">
              <div class="min-w-0">
                <div class="text-sm font-medium">
                  {job.action}
                  <span class="ml-2 text-xs font-normal capitalize text-muted-foreground">
                    {job.status} · {parallelCount > 1
                      ? `${parallelCount} workers · `
                      : ''}{job.elapsed_seconds}s
                  </span>
                </div>
                {#if job.error}
                  <p class="truncate text-xs text-destructive">{job.error}</p>
                {/if}
              </div>
              <div class="flex flex-wrap items-center gap-2">
                {#if jobRunning}
                  <Button variant="destructive" size="sm" onclick={cancelJob}>Cancel</Button>
                {/if}
                {#if job.id}
                  <Button
                    variant="outline"
                    size="sm"
                    onclick={() => {
                      historySelectedJobId = job.id;
                      historyOpen = true;
                    }}
                  >
                    Logs & Details
                  </Button>
                {/if}
                {#if !jobRunning}
                  <Button variant="ghost" size="sm" onclick={() => (dismissedJobId = job.id)}>
                    Dismiss
                  </Button>
                {/if}
              </div>
            </div>
            {#if job.stages.length}
              <JobPhases stages={job.stages} {nowSeconds} />
            {/if}
          </div>
        {/if}
      </div>
    {/snippet}

    {#snippet content()}
      <div class="pb-4">
        {#if tab === 'resolver'}
          {#if !status}
            <div
              class="rounded-[1.25rem] border border-border/60 bg-card px-6 py-14 text-center text-muted-foreground shadow-sm"
            >
              Loading hubs…
            </div>
          {:else if searchedResolver.length > 0}
            <div class="grid grid-cols-1 items-start gap-4 xl:grid-cols-2 2xl:grid-cols-3">
              {#each searchedResolver as card (card.id)}
                {@const active = jobTouchesHub(card.id)}
                <AdminCard
                  title={card.title}
                  subtitle={card.file.exists
                    ? `${fmtBytes(card.file.size_bytes)} parquet`
                    : 'Not built'}
                  actionLabel={card.file.exists ? 'Rebuild' : 'Build'}
                  onAction={() => rebuildHub(card.id)}
                  disabled={!canBuild}
                  busy={active}
                  statusMessage={active ? runningStage?.message || '' : ''}
                  primaryLabel="Rows"
                  primaryValue={card.file.rows}
                  primaryMax={maxHubRows}
                  files={card.file.exists
                    ? [
                        {
                          name: `${card.id}.parquet`,
                          size_bytes: card.file.size_bytes,
                          rows: card.file.rows,
                        },
                      ]
                    : []}
                  hasLogs={Boolean(card.file.log?.exists) || active}
                  onOpenLogs={() => openLogs('hub', card.id, `${card.title} log`)}
                />
              {/each}
            </div>
          {:else}
            <div
              class="rounded-[1.25rem] border border-border/60 bg-card px-6 py-14 text-center text-muted-foreground shadow-sm"
            >
              {#if query.trim()}
                No hubs matched the current search.
              {:else}
                No resolver hubs found.
              {/if}
            </div>
          {/if}
        {:else if !status}
          <div
            class="rounded-[1.25rem] border border-border/60 bg-card px-6 py-14 text-center text-muted-foreground shadow-sm"
          >
            Loading resources…
          </div>
        {:else if searchedResources.length > 0}
          <div class="grid grid-cols-1 items-start gap-4 xl:grid-cols-2 2xl:grid-cols-3">
            {#each searchedResources as card (card.id)}
              {@const built = card.built}
              {@const active = jobTouchesResource(card.id)}
              <AdminCard
                title={card.title}
                subtitle={resourceSubtitle(built, card.versions.length)}
                actionLabel={built ? 'Rebuild' : 'Build'}
                onAction={() => rebuildResource(card.id)}
                disabled={!canBuild}
                busy={active}
                statusMessage={active ? runningStage?.message || '' : ''}
                primaryLabel="Entities"
                primaryValue={built?.entities.rows || 0}
                primaryMax={maxEntities}
                secondaryLabel="Relations"
                secondaryValue={built?.relations.rows || 0}
                secondaryMax={maxRelations}
                versionsCount={card.versions.length}
                onManageVersions={() => {
                  selectedResourceName = card.id;
                  versionsOpen = true;
                }}
                files={built
                  ? [
                      {
                        name: 'entity.parquet',
                        size_bytes: built.entities.size_bytes,
                        rows: built.entities.rows,
                      },
                      {
                        name: 'relation.parquet',
                        size_bytes: built.relations.size_bytes,
                        rows: built.relations.rows,
                      },
                      {
                        name: 'evidence_payloads.parquet',
                        size_bytes: built.payloads.size_bytes,
                        rows: built.payloads.rows,
                      },
                    ]
                  : []}
                hasLogs={Boolean(status?.resources.logs?.[card.id]?.exists || built?.log?.exists) ||
                  active}
                onOpenLogs={() => openLogs('resource', card.id, `${card.title} log`)}
              />
            {/each}
          </div>
        {:else}
          <div
            class="rounded-[1.25rem] border border-border/60 bg-card px-6 py-14 text-center text-muted-foreground shadow-sm"
          >
            {#if !resourceCards.length}
              No parquet resources are available yet.
            {:else}
              No resources matched the current search.
            {/if}
          </div>
        {/if}
      </div>
    {/snippet}
  </ExploreBrowserShell>

  <LogDialog bind:open={logOpen} title={logTitle} text={logText} loading={logLoading} />
  <AdminHistorySheet
    bind:open={historyOpen}
    jobs={allJobs}
    bind:selectedJobId={historySelectedJobId}
    {nowSeconds}
    onRefresh={() => void refresh()}
    onCancelJob={(id) => void cancelJobById(id)}
  />
  <ResourceVersionsDialog
    bind:open={versionsOpen}
    resourceName={selectedResourceName}
    versions={selectedResourceVersions}
    activeVersion={selectedResourceActiveVersion}
    onDeleteVersion={deleteResourceVersion}
    onCleanupOlder={cleanupOlderVersions}
    onViewLogs={(res, ver) => openLogs('resource', res, `${res}/${ver} log`, ver)}
  />
{/if}
