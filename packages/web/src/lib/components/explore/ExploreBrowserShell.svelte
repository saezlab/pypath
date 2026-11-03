<script lang="ts">
  import { Info, Search } from '@lucide/svelte';
  import { getUiPreferences } from '$lib/stores/ui-preferences.svelte';

  interface Props {
    fitViewport?: boolean;
    query: string;
    draftQuery: string;
    onDraftQueryChange: (value: string) => void;
    onSubmitSearch: () => void;
    tab?: string;
    onTabChange?: (tab: string) => void;
    tabs?: Array<{ value: string; label: string; badge?: string | number | null }>;
    filterRail?: boolean;
    content: import('svelte').Snippet;
    searchActions?: import('svelte').Snippet;
    searchPlaceholder: string;
    searchInputRef?: HTMLInputElement | null;
    footerCta?: import('svelte').Snippet;
    summarySlot?: import('svelte').Snippet;
    explanationTitle?: string;
    explanationText?: string;
  }

  let {
    fitViewport = false,
    draftQuery,
    onDraftQueryChange,
    onSubmitSearch,
    tab = '',
    onTabChange = () => {},
    tabs = [],
    content,
    searchActions,
    searchPlaceholder,
    searchInputRef = $bindable(null),
    footerCta,
    summarySlot,
    explanationTitle = '',
    explanationText = '',
    filterRail = false,
  }: Props = $props();

  const uiPreferences = getUiPreferences();

  function handleKeyDown(event: KeyboardEvent) {
    if (event.key === 'Enter') {
      event.preventDefault();
      onSubmitSearch();
    }
  }
</script>

<div
  class={`relative flex w-full flex-col ${fitViewport ? 'h-full min-h-0 gap-2 overflow-hidden' : 'gap-4'}`}
>
  <div class="sticky top-0 z-20 -mx-4 shrink-0 bg-background px-4 pt-2 pb-0 md:-mx-6 md:px-6">
    <div class="flex w-full flex-wrap items-center gap-2 sm:flex-nowrap">
      {#if tabs.length}
        <div
          class="flex w-full shrink-0 items-center gap-0.5 overflow-x-auto rounded-md bg-muted/40 p-0.5 sm:w-auto [scrollbar-width:none] [&::-webkit-scrollbar]:hidden"
        >
          {#each tabs as item}
            <button
              type="button"
              onclick={() => onTabChange(item.value)}
              aria-pressed={tab === item.value}
              class={`inline-flex h-8 shrink-0 items-center gap-1.5 rounded px-3 text-[13px] font-medium capitalize whitespace-nowrap transition-colors ${
                tab === item.value
                  ? 'bg-background text-foreground shadow-xs'
                  : 'text-muted-foreground hover:bg-muted/40 hover:text-foreground'
              }`}
            >
              <span>{item.label}</span>
              {#if item.badge}
                <span class="text-xs opacity-70">({item.badge})</span>
              {/if}
            </button>
          {/each}
        </div>
      {/if}
      <div
        class="flex min-w-0 flex-1 items-center gap-2 rounded-md border border-border/60 bg-background px-3 focus-within:border-ring"
      >
        <Search class="h-4 w-4 shrink-0 text-muted-foreground" />
        <input
          type="text"
          bind:this={searchInputRef}
          value={draftQuery}
          oninput={(e) => onDraftQueryChange(e.currentTarget.value)}
          onkeydown={handleKeyDown}
          placeholder={searchPlaceholder}
          aria-label={searchPlaceholder}
          class="h-9 min-w-0 flex-1 bg-transparent text-sm text-foreground placeholder:text-muted-foreground focus:outline-none"
        />
        <button
          type="button"
          onclick={onSubmitSearch}
          aria-label="Search"
          title="Search (Enter)"
          class="inline-flex h-6 shrink-0 items-center rounded border border-border/70 px-1.5 font-mono text-[11px] text-muted-foreground transition-colors hover:border-ring hover:text-foreground"
        >
          ↵
        </button>
      </div>
      {#if searchActions}{@render searchActions()}{/if}
    </div>
  </div>

  {#if summarySlot}
    <div class={filterRail ? 'md:pr-72' : ''}>
      {@render summarySlot()}
    </div>
  {/if}

  {#if uiPreferences.showExplanations && explanationText}
    <div class={`px-1 ${filterRail ? 'md:pr-72' : ''}`}>
      <div class="flex gap-3">
        <Info class="mt-0.5 size-4 shrink-0 text-primary" />
        <div class="min-w-0 space-y-1.5">
          {#if explanationTitle}
            <h2 class="text-sm font-semibold leading-5 text-foreground">{explanationTitle}</h2>
          {/if}
          <p class="text-sm leading-6 text-muted-foreground">{explanationText}</p>
        </div>
      </div>
    </div>
  {/if}

  <div class={fitViewport ? 'flex min-h-0 flex-1 flex-col overflow-hidden' : ''}>
    {@render content()}
  </div>

  {#if footerCta}
    {@render footerCta()}
  {/if}
</div>
