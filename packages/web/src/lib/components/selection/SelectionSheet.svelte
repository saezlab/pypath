<script lang="ts">
  import { ArrowRight, Database, Settings2, Tag, Trash2, X } from '@lucide/svelte';
  import { Badge } from '$lib/components/ui/badge/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Switch } from '$lib/components/ui/switch/index.js';
  import {
    Sheet,
    SheetContent,
    SheetHeader,
    SheetTrigger,
  } from '$lib/components/ui/sheet/index.js';
  import { buildSelectionUrl } from '$lib/navigation/url-codecs';
  import { getSelectionStore } from '$lib/stores/selection.svelte';
  import { getSelectionScopeSettings } from '$lib/stores/selection-scope.svelte';

  interface Props {
    open?: boolean;
    triggerClass?: string;
    /** Without a trigger the sheet is opened through `open` only. */
    showTrigger?: boolean;
  }

  let { open = $bindable(false), triggerClass = '', showTrigger = true }: Props = $props();

  const selection = getSelectionStore();
  const scopeSettings = getSelectionScopeSettings();
  const totalSelectionCount = $derived(
    selection.selectedEntities.length + selection.selectedAnnotations.length,
  );
  const selectionHref = $derived(
    buildSelectionUrl({
      entityIds: selection.entityIds,
      annotationIds: selection.annotationIds,
    }),
  );
</script>

{#if totalSelectionCount > 0}
  <Sheet bind:open>
    {#if showTrigger}
      <SheetTrigger>
        {#snippet child({ props })}
          <Button {...props} variant="outline" class={triggerClass}>
            <Settings2 class="size-4" />
            Manage selection
          </Button>
        {/snippet}
      </SheetTrigger>
    {/if}
    <SheetContent side="right" class="w-[92vw] gap-0 overflow-hidden p-0 sm:max-w-xl">
      <SheetHeader class="border-b px-6 py-4 pr-12">
        <Button href={selectionHref} variant="outline" size="sm" class="w-fit gap-1.5">
          Selection page
          <ArrowRight class="size-3.5" />
        </Button>
      </SheetHeader>
      <div class="flex min-h-0 flex-1 flex-col">
        <div class="border-b bg-muted/20 px-6 py-4">
          <div class="grid grid-cols-2 gap-3">
            <div class="rounded-xl border bg-background/70 p-3">
              <div class="flex items-center gap-2 text-sm font-medium">
                <Database class="size-4 text-primary" />
                Entities
              </div>
              <div class="mt-1 text-2xl font-semibold tabular-nums">
                {selection.selectedEntities.length}
              </div>
            </div>
            <div class="rounded-xl border bg-background/70 p-3">
              <div class="flex items-center gap-2 text-sm font-medium">
                <Tag class="size-4 text-primary" />
                CV terms
              </div>
              <div class="mt-1 text-2xl font-semibold tabular-nums">
                {selection.selectedAnnotations.length}
              </div>
            </div>
          </div>

          <div class="mt-4 rounded-xl border bg-background/70 p-3">
            <div class="mb-3 flex items-center gap-2 text-sm font-medium">
              <Settings2 class="size-4 text-primary" />
              Selection settings
            </div>

            <div class="space-y-3">
              <div class="space-y-1.5">
                <div class="text-xs font-medium text-muted-foreground">Scope operator</div>
                <div class="inline-flex w-full rounded-md border bg-muted/30 p-0.5">
                  <Button
                    type="button"
                    variant={scopeSettings.mode === 'union' ? 'default' : 'ghost'}
                    size="sm"
                    class="h-8 flex-1 rounded-sm px-3"
                    onclick={() => scopeSettings.setMode('union')}
                  >
                    Union
                  </Button>
                  <Button
                    type="button"
                    variant={scopeSettings.mode === 'intersection' ? 'default' : 'ghost'}
                    size="sm"
                    class="h-8 flex-1 rounded-sm px-3"
                    onclick={() => scopeSettings.setMode('intersection')}
                  >
                    Intersection
                  </Button>
                </div>
              </div>

              <label
                class="flex items-center justify-between gap-3 rounded-md border bg-muted/20 px-3 py-2"
              >
                <span class="min-w-0 text-sm">
                  <span class="block font-medium"
                    >Include associated entities, members, and participants</span
                  >
                </span>
                <Switch
                  size="sm"
                  checked={scopeSettings.expandSelection}
                  onCheckedChange={(checked) => scopeSettings.setExpandSelection(checked)}
                />
              </label>
            </div>
          </div>
        </div>

        <div class="min-h-0 flex-1 overflow-y-auto px-6 py-5">
          <div class="mb-5 flex items-center justify-between gap-3">
            <div>
              <h3 class="text-sm font-semibold">Selected items</h3>
              <p class="text-xs text-muted-foreground">
                Remove individual entries or clear the whole selection.
              </p>
            </div>
            <Button
              variant="outline"
              size="sm"
              onclick={selection.clearSelection}
              class="h-8 shrink-0 gap-1.5"
            >
              <Trash2 class="size-3.5" />
              Clear all
            </Button>
          </div>

          <div class="space-y-6">
            {#if selection.selectedEntities.length > 0}
              <section class="space-y-2.5">
                <div
                  class="flex items-center gap-2 text-xs font-semibold uppercase tracking-wide text-muted-foreground"
                >
                  <Database class="size-3.5" />
                  Entities
                </div>
                <div class="space-y-2">
                  {#each selection.selectedEntities as entity}
                    <div
                      class="group flex items-center gap-3 rounded-xl border bg-card p-3 shadow-sm transition-colors hover:bg-muted/30"
                    >
                      <div
                        class="flex size-9 shrink-0 items-center justify-center rounded-full bg-primary/10 text-sm font-semibold text-primary"
                      >
                        {(entity.name || entity.id).slice(0, 2).toUpperCase()}
                      </div>
                      <div class="min-w-0 flex-1">
                        <div class="truncate text-sm font-medium">
                          {entity.name || entity.id}
                        </div>
                        <div class="mt-1 flex items-center gap-2">
                          {#if entity.type}
                            <Badge
                              variant="secondary"
                              class="h-5 rounded-md px-1.5 text-[10px] uppercase"
                            >
                              {entity.type}
                            </Badge>
                          {/if}
                          <span class="truncate font-mono text-xs text-muted-foreground">
                            {entity.id}
                          </span>
                        </div>
                      </div>
                      <Button
                        type="button"
                        variant="ghost"
                        size="icon"
                        aria-label={`Remove ${entity.name || entity.id}`}
                        onclick={() => selection.removeEntity(entity.id)}
                        class="size-8 shrink-0 text-muted-foreground hover:text-foreground"
                      >
                        <X class="size-4" />
                      </Button>
                    </div>
                  {/each}
                </div>
              </section>
            {/if}

            {#if selection.selectedAnnotations.length > 0}
              <section class="space-y-2.5">
                <div
                  class="flex items-center gap-2 text-xs font-semibold uppercase tracking-wide text-muted-foreground"
                >
                  <Tag class="size-3.5" />
                  CV terms
                </div>
                <div class="space-y-2">
                  {#each selection.selectedAnnotations as annotation}
                    <div
                      class="group flex items-center gap-3 rounded-xl border bg-card p-3 shadow-sm transition-colors hover:bg-muted/30"
                    >
                      <div
                        class="flex size-9 shrink-0 items-center justify-center rounded-full bg-muted text-sm font-semibold text-muted-foreground"
                      >
                        {annotation.namespace?.slice(0, 2).toUpperCase() || 'CV'}
                      </div>
                      <div class="min-w-0 flex-1">
                        <div class="truncate text-sm font-medium">
                          {annotation.label || annotation.id}
                        </div>
                        <div class="mt-1 flex items-center gap-2">
                          {#if annotation.namespace}
                            <Badge
                              variant="outline"
                              class="h-5 rounded-md px-1.5 font-mono text-[10px] uppercase"
                            >
                              {annotation.namespace}
                            </Badge>
                          {/if}
                          <span class="truncate font-mono text-xs text-muted-foreground">
                            {annotation.id}
                          </span>
                        </div>
                      </div>
                      <Button
                        type="button"
                        variant="ghost"
                        size="icon"
                        aria-label={`Remove ${annotation.label || annotation.id}`}
                        onclick={() => selection.removeAnnotation(annotation.id)}
                        class="size-8 shrink-0 text-muted-foreground hover:text-foreground"
                      >
                        <X class="size-4" />
                      </Button>
                    </div>
                  {/each}
                </div>
              </section>
            {/if}
          </div>
        </div>
      </div>
    </SheetContent>
  </Sheet>
{/if}
