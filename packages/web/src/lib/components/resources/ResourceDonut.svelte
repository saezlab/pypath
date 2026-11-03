<script lang="ts">
  import type { Snippet } from 'svelte';
  import type { ResourceRecord } from '$lib/resources/types';

  let {
    resource,
    size = 26,
    ring = 6,
    children,
  }: { resource: ResourceRecord; size?: number; ring?: number; children?: Snippet } = $props();

  const total = $derived(resource.entity_count + resource.interaction_count);
  const entityPercent = $derived(total > 0 ? (resource.entity_count / total) * 100 : 0);
</script>

<!-- The share of entity (cyan) and relation (green) records. -->
<span
  class="donut"
  style:width={`${size}px`}
  style:height={`${size}px`}
  role="img"
  aria-label={`${entityPercent.toFixed(0)}% entity records, ${(total ? 100 - entityPercent : 0).toFixed(0)}% relation records`}
>
  <span
    class="ring"
    style:--ring={`${ring}px`}
    style:background={total
      ? `conic-gradient(var(--resource-entities) 0 ${entityPercent}%, var(--resource-relations) ${entityPercent}% 100%)`
      : 'var(--border)'}
  ></span>
  {#if children}<span class="center">{@render children()}</span>{/if}
</span>

<style>
  .donut {
    position: relative;
    flex: none;
    display: inline-grid;
    place-items: center;
  }
  .ring {
    position: absolute;
    inset: 0;
    border-radius: 50%;
    /* A mask cuts the hole, so the ring sits on any background. */
    mask: radial-gradient(
      farthest-side,
      transparent calc(100% - var(--ring)),
      #000 calc(100% - var(--ring) + 0.5px)
    );
  }
  .center {
    position: relative;
  }
</style>
