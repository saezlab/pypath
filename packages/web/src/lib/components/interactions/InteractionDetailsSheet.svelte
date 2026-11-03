<script lang="ts">
  import {
    Dialog,
    DialogContent,
    DialogHeader,
    DialogTitle,
  } from '$lib/components/ui/dialog/index.js';
  import InteractionDetails from './InteractionDetails.svelte';
  import { fetchRelationEvidence } from '$lib/api/client';
  import type { InteractionDetailsData, InteractionListRow } from '$lib/types/interactions';
  import type { SearchFilters } from '$lib/types/search';
  import type { ProductSummary } from '$lib/utils/molecular-presentation';

  interface Props {
    open: boolean;
    onOpenChange: (open: boolean) => void;
    interaction: InteractionListRow | null;
    /** Narrows the evidence to records matching these filters, e.g. one product. */
    filters?: SearchFilters;
    /** Known products, so molecular forms show their names instead of a bare reference. */
    products?: readonly ProductSummary[];
  }

  let { open, onOpenChange, interaction, filters, products = [] }: Props = $props();

  let dialogOpen = $state(false);
  let evidence = $state<InteractionDetailsData['evidence']>([]);
  let evidenceLoading = $state(false);
  let evidenceError = $state<string | null>(null);

  const relationPk = $derived(interaction?.relation.relationPk);

  $effect(() => {
    dialogOpen = open;
  });

  $effect(() => {
    if (dialogOpen !== open) onOpenChange(dialogOpen);
  });

  $effect(() => {
    evidence = [];
    if (!open || !relationPk) return;
    let cancelled = false;
    evidenceLoading = true;
    evidenceError = null;
    fetchRelationEvidence(relationPk, filters)
      .then((data) => {
        if (!cancelled) evidence = data.evidence;
      })
      .catch((err) => {
        if (!cancelled)
          evidenceError = err instanceof Error ? err.message : 'Failed to load evidence';
      })
      .finally(() => {
        if (!cancelled) evidenceLoading = false;
      });
    return () => {
      cancelled = true;
    };
  });

  const details = $derived<InteractionDetailsData | InteractionListRow | null>(
    interaction ? { ...interaction, evidence } : null,
  );
</script>

<Dialog bind:open={dialogOpen}>
  <DialogContent
    class="grid max-h-[88vh] grid-rows-[auto_minmax(0,1fr)] gap-0 overflow-hidden p-0 sm:max-w-4xl"
  >
    <DialogHeader class="border-b px-6 py-4 pr-14">
      <DialogTitle class="text-lg">Relation details</DialogTitle>
    </DialogHeader>

    <div class="min-h-0 overflow-y-auto overscroll-contain">
      {#if evidenceError}
        <div class="p-6 text-sm text-destructive">{evidenceError}</div>
      {:else}
        <InteractionDetails selectedInteraction={details} {evidenceLoading} {products} />
      {/if}
    </div>
  </DialogContent>
</Dialog>
