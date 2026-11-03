<script lang="ts">
  import { onMount } from 'svelte';
  import { workspaceWidth } from '$lib/stores/workspace-width.svelte';
  import { page } from '$app/state';
  import { releaseUrl } from '$lib/api/release';
  import { mode, setMode } from 'mode-watcher';
  import { Check, Info, Menu, Moon, Settings2, Sun } from '@lucide/svelte';
  import {
    DropdownMenu,
    DropdownMenuCheckboxItem,
    DropdownMenuContent,
    DropdownMenuItem,
    DropdownMenuLabel,
    DropdownMenuSeparator,
    DropdownMenuSub,
    DropdownMenuSubTrigger,
    DropdownMenuSubContent,
    DropdownMenuRadioGroup,
    DropdownMenuRadioItem,
    DropdownMenuTrigger,
  } from '$lib/components/ui/dropdown-menu/index.js';
  import { getUiPreferences } from '$lib/stores/ui-preferences.svelte';
  import { getSelectionStore } from '$lib/stores/selection.svelte';
  import { cn } from '$lib/utils.js';
  import ServiceStatus from '$lib/components/layout/ServiceStatus.svelte';

  onMount(() => workspaceWidth.load());
  const selection = getSelectionStore();
  const navigationItems = $derived([
    { title: 'Explore', url: '/explore' },
    { title: 'Selection', url: '/selection', count: selection.totalSelectionCount },
    { title: 'Resources', url: '/resources' },
    { title: 'Skills', url: '/skills' },
    { title: 'API Docs', url: '/api/docs', external: true },
  ]);

  const uiPreferences = getUiPreferences();
  let darkMode = $state(false);

  $effect(() => {
    darkMode = mode.current === 'dark';
  });

  function isPathActive(url: string) {
    return page.url.pathname === url || page.url.pathname.startsWith(`${url}/`);
  }

  function selectRelease(version: string) {
    const url = new URL(window.location.href);
    url.searchParams.set('release', version);
    // A fresh page also discards query caches and details from the previous release.
    window.location.assign(url);
  }
</script>

<a
  href="#main-content"
  class="sr-only fixed left-3 top-3 z-50 rounded-md bg-background px-3 py-2 text-sm font-medium shadow-lg focus:not-sr-only"
>
  Skip to content
</a>

<header class="z-40 shrink-0 bg-background">
  <div class="flex h-12 items-center gap-6">
    <div class="flex min-w-0 items-center gap-2 justify-self-start">
      <a
        href={releaseUrl('/explore')}
        class="flex shrink-0 items-center gap-2 justify-self-start rounded-md focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
        aria-label="OmniPath home"
      >
        <img src="/omnipath-logo-gradient.svg" alt="" width="30" height="30" />
        <span class="text-base font-semibold tracking-tight text-foreground">OmniPath</span>
      </a>
    </div>

    <nav aria-label="Primary" class="hidden h-full min-w-0 items-center gap-1 lg:flex">
      {#each navigationItems.filter((item) => !item.external) as item (item.title)}
        <a
          href={item.external ? item.url : releaseUrl(item.url)}
          target={item.external ? '_blank' : undefined}
          rel={item.external ? 'noopener noreferrer' : undefined}
          aria-current={isPathActive(item.url) ? 'page' : undefined}
          class={cn(
            'flex h-full items-center px-2.5 text-[13px] transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-ring',
            isPathActive(item.url)
              ? 'font-medium text-foreground'
              : 'text-muted-foreground hover:text-foreground',
          )}
        >
          <span>{item.title}</span>
          {#if item.count}<span
              class="ml-1.5 inline-grid h-[18px] min-w-[18px] place-items-center rounded-full bg-primary px-1.5 text-[11px] font-semibold text-primary-foreground tabular-nums"
              aria-label={`${item.count} selected`}>{item.count}</span
            >{/if}
        </a>
      {/each}
    </nav>

    <div class="ml-auto flex items-center gap-2">
      <ServiceStatus />
      <a
        href="/api/docs"
        target="_blank"
        rel="noopener noreferrer"
        class="hidden h-8 items-center gap-1 rounded-md px-2 text-xs text-muted-foreground hover:text-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring lg:inline-flex"
        >Documentation ↗</a
      >
      <label for="release-select" class="hidden text-xs text-muted-foreground/70 sm:inline"
        >Release</label
      >
      <select
        id="release-select"
        aria-label="OmniPath release"
        class="max-w-32 cursor-pointer rounded-md border border-transparent bg-transparent px-1.5 py-1 text-xs text-muted-foreground hover:text-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
        value={page.data.selectedRelease || 'latest'}
        disabled={page.data.releasesUnavailable}
        onchange={(event) => selectRelease(event.currentTarget.value)}
      >
        {#if page.data.releasesUnavailable}
          <option value={page.data.selectedRelease || 'latest'}>Versions unavailable</option>
        {:else}
          <option value="latest">Latest</option>
          {#each page.data.releases || [] as release}
            <option value={release.version}>{release.version}</option>
          {/each}
        {/if}
      </select>
      <DropdownMenu>
        <DropdownMenuTrigger
          class="inline-flex size-9 items-center justify-center rounded-md text-muted-foreground transition-colors hover:bg-muted hover:text-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring lg:hidden"
          aria-label="Open navigation"
        >
          <Menu class="h-5 w-5" />
        </DropdownMenuTrigger>
        <DropdownMenuContent align="end" class="w-56">
          <DropdownMenuLabel>Navigate</DropdownMenuLabel>
          <DropdownMenuSeparator />
          {#each navigationItems as item (item.title)}
            <DropdownMenuItem>
              {#snippet child({ props })}
                <a
                  href={item.external ? item.url : releaseUrl(item.url)}
                  target={item.external ? '_blank' : undefined}
                  rel={item.external ? 'noopener noreferrer' : undefined}
                  aria-current={isPathActive(item.url) ? 'page' : undefined}
                  class="flex w-full items-center gap-2"
                  {...props}
                >
                  <span class={isPathActive(item.url) ? 'text-primary' : undefined}
                    >{item.title}</span
                  >
                  {#if item.count}<span class="ml-auto text-xs text-muted-foreground tabular-nums"
                      >{item.count}</span
                    >{/if}
                </a>
              {/snippet}
            </DropdownMenuItem>
          {/each}
        </DropdownMenuContent>
      </DropdownMenu>

      <DropdownMenu>
        <DropdownMenuTrigger
          class="inline-flex size-9 items-center justify-center rounded-md text-muted-foreground transition-colors hover:bg-muted hover:text-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
          aria-label="Open display preferences"
        >
          <Settings2 class="h-4 w-4" />
        </DropdownMenuTrigger>
        <DropdownMenuContent align="end" class="w-64">
          <DropdownMenuLabel>Display preferences</DropdownMenuLabel>
          <DropdownMenuSub>
            <DropdownMenuSubTrigger>Content width</DropdownMenuSubTrigger>
            <DropdownMenuSubContent>
              <DropdownMenuRadioGroup
                value={String(workspaceWidth.value)}
                onValueChange={(value) => workspaceWidth.set(Number(value))}
              >
                <DropdownMenuRadioItem value="1280">Default</DropdownMenuRadioItem>
                <DropdownMenuRadioItem value="1600">Wide</DropdownMenuRadioItem>
                <DropdownMenuRadioItem value="1920">Extra wide</DropdownMenuRadioItem>
                <DropdownMenuRadioItem value="3840">Full width</DropdownMenuRadioItem>
              </DropdownMenuRadioGroup>
            </DropdownMenuSubContent>
          </DropdownMenuSub>
          <DropdownMenuSeparator />
          <DropdownMenuCheckboxItem
            checked={uiPreferences.showExplanations}
            onCheckedChange={(checked) => uiPreferences.setShowExplanations(checked)}
          >
            <Info class="h-4 w-4" />
            Show explanations
          </DropdownMenuCheckboxItem>
          <DropdownMenuSeparator />
          <DropdownMenuItem onclick={() => setMode('light')}>
            <Sun class="h-4 w-4" />
            Light
            {#if !darkMode}<Check class="ml-auto h-4 w-4" />{/if}
          </DropdownMenuItem>
          <DropdownMenuItem onclick={() => setMode('dark')}>
            <Moon class="h-4 w-4" />
            Dark
            {#if darkMode}<Check class="ml-auto h-4 w-4" />{/if}
          </DropdownMenuItem>
        </DropdownMenuContent>
      </DropdownMenu>
    </div>
  </div>
</header>
