<script lang="ts">
  import { ArrowDown, ArrowUp } from '@lucide/svelte';
  import { getRelationPredicateLabel } from '$lib/domain/display';
  import { relationPhrase, type PhraseRelation } from '$lib/domain/relation-phrase';

  // `responsive`: inside a size container (`@container`), the short phrase in a narrow one and
  // the full phrase from 32rem (the full phrase is also the tooltip).
  let { relation, responsive = false }: { relation: PhraseRelation; responsive?: boolean } =
    $props();

  const phrase = $derived(
    relationPhrase(
      relation,
      getRelationPredicateLabel(relation.predicate, relation.sign ?? 0) || '—',
    ),
  );
</script>

<span
  class="inline-flex max-w-full flex-col items-center text-center leading-tight"
  title={responsive ? [phrase.text, phrase.detail].filter(Boolean).join(' ') : undefined}
>
  <!-- The arrow is inline so it stays with the first word when the phrase wraps. A word too
       long for the column breaks rather than running into the next one. -->
  <span class="max-w-full text-sm text-balance text-foreground/80 [overflow-wrap:anywhere]">
    {#if phrase.direction === 'up'}<ArrowUp
        class="mr-0.5 inline-block size-3.5 align-[-2px] text-green-600 dark:text-green-400"
        aria-label="increased"
      />{:else if phrase.direction === 'down'}<ArrowDown
        class="mr-0.5 inline-block size-3.5 align-[-2px] text-red-600 dark:text-red-400"
        aria-label="decreased"
      />{/if}{#if responsive}<span class="@lg:hidden">{phrase.short}</span><span
        class="hidden @lg:inline">{phrase.text}</span
      >{:else}{phrase.text}{/if}
  </span>
  {#if phrase.detail}<span
      class="max-w-full text-xs text-muted-foreground [overflow-wrap:anywhere] {responsive
        ? 'hidden @lg:inline'
        : ''}">{phrase.detail}</span
    >{/if}
</span>
