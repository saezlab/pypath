<script lang="ts">
  import './layout.css';
  import { page } from '$app/state';
  import { workspaceWidth } from '$lib/stores/workspace-width.svelte';
  import { ModeWatcher } from 'mode-watcher';
  import { Toaster } from 'svelte-sonner';
  import AppHeader from '$lib/components/layout/AppHeader.svelte';

  let { children } = $props();
  const isWorkspace = $derived(
    ['/explore', '/selection', '/resources'].includes(page.url.pathname),
  );
</script>

<svelte:head>
  <link rel="icon" type="image/svg+xml" href="/omnipath-logo-gradient.svg" />
  <title>OmniPath Explorer</title>
  <meta
    name="description"
    content="Explore molecular interactions, pathways, and biological annotations"
  />
</svelte:head>

<ModeWatcher defaultMode="system" />
<div class="flex h-svh w-full flex-col overflow-hidden bg-background text-foreground">
  <div class="mx-auto w-full max-w-7xl shrink-0 px-4 pt-2 md:px-6">
    <AppHeader />
  </div>
  <main
    id="main-content"
    class={`min-h-0 w-full flex-1 ${isWorkspace ? 'overflow-hidden' : 'overflow-y-auto'}`}
    tabindex="-1"
  >
    <div
      class={`mx-auto w-full px-4 md:px-6 ${isWorkspace ? 'flex h-full min-h-0 flex-col pb-2' : 'pb-5'}`}
      style:max-width={`${isWorkspace ? workspaceWidth.value : 1280}px`}
    >
      {@render children()}
    </div>
  </main>
</div>
<Toaster richColors />
