<script lang="ts">
  import { Check, Link, LoaderCircle, Network, Plus } from '@lucide/svelte';
  import { truncateTitle } from '$lib/actions/truncate-title';
  import EntityDetailsDialog from './EntityDetailsDialog.svelte';
  import { getSelectionStore } from '$lib/stores/selection.svelte';
  import {
    getEntityDisplayName,
    getEntityPublicId,
    getEntitySecondaryName,
    getEntityTypeLabel,
    getIdentifierTypeLabel,
    getIdentifierDisplayTypeForValue,
    isChemicalEntity,
    isCvTermEntity,
    isUnresolvedEntity,
  } from '$lib/domain/display';
  import { getEntityTypeEmoji, getEntityTypeStyle } from '$lib/utils/entity-types';
  import { releaseFetch } from '$lib/api/release';
  import { formatNumber } from '$lib/utils/format';
  import { formatTaxonomy } from '$lib/utils/taxonomy';
  import type { EntityWithIdentifiers as EntityResult } from '$lib/types/entities';
  let { result, onopen }: { result: EntityResult; onopen?: (entity: EntityResult) => void } =
    $props();
  const selection = getSelectionStore();
  let detailsOpen = $state(false);
  let selectingGroup = $state(false);
  let selectionError = $state('');
  async function toggleGroup() {
    if (selectingGroup) return;
    selectingGroup = true;
    selectionError = '';
    try {
      const keys = [...(result.groupMemberKeys || [])];
      let cursor = result.groupMemberCursor;
      while (cursor) {
        const response = await releaseFetch('/app-api/entities/groups', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            strategy: result.groupStrategy || 'chemical_connectivity',
            group_key: result.entityPk,
            query: result.groupQuery,
            filters: result.groupFilters,
            resources: result.groupResources,
            member_cursor: cursor,
            member_limit: 100,
          }),
        });
        if (!response.ok) throw new Error('Unable to select the group. Please retry.');
        const group = (await response.json()).groups[0];
        if (!group) throw new Error('This group is no longer available.');
        keys.push(...group.entity.groupMemberKeys);
        cursor = group.nextMemberCursor;
      }
      selection.toggleEntityGroup(keys, result.entityPk);
    } catch (error) {
      selectionError = (error as Error).message;
    } finally {
      selectingGroup = false;
    }
  }
  function openDetails(entity: EntityResult) {
    if (onopen) onopen(entity);
    else detailsOpen = true;
  }
  function handleDetailsKeydown(event: KeyboardEvent, entity: EntityResult) {
    if (event.key !== 'Enter' && event.key !== ' ') return;
    event.preventDefault();
    openDetails(entity);
  }
  function getEntityOntologyTerm(entity: EntityResult) {
    if (entity.ontologyHierarchy) return entity.ontologyHierarchy;
    if (!isCvTermEntity(entity)) return null;
    const id = entity.canonicalIdentifier || '';
    if (!id) return null;
    const prefix = id.includes(':') ? id.split(':')[0] : id.includes('-') ? id.split('-')[0] : null;
    return {
      termId: id,
      ontologyPrefix: prefix,
      label: getEntityDisplayName(entity),
      definition: null,
      ontologyId: null,
      childCount: 0,
      parentCount: 0,
    };
  }

  function isPlainNameIdentifierType(identifierType: string | undefined): boolean {
    const label = getIdentifierTypeLabel(identifierType).toLowerCase();
    const normalized = [identifierType || '', label].join(' ').toLowerCase().replace(/[_-]/g, ' ');
    return (
      label === 'name' ||
      normalized.includes('gene name') ||
      normalized.includes('gene symbol') ||
      normalized.includes('recommended name') ||
      normalized.includes('entry name')
    );
  }

  function firstHint(values: string[] | undefined): string | undefined {
    return values?.map((value) => value.trim()).find(Boolean);
  }

  function chemicalHint(entity: EntityResult): string | undefined {
    if (!isChemicalEntity(entity)) return undefined;
    const hints = entity.entityFacetHints;
    return (
      firstHint(hints?.structuralSpecificities) ||
      firstHint(hints?.chemicalClasses) ||
      firstHint(hints?.metabolicDomains)
    );
  }

  function isHashLike(value: string): boolean {
    return /^[a-f0-9-]{24,}$/i.test(value) && /[a-f]/i.test(value);
  }

  /** The identifier shown after the name: the secondary identifier, else the canonical one. */
  function rowIdentifier(
    entity: EntityResult,
    displayName: string,
  ): { label: string; value: string } | null {
    if (entity.groupStrategy === 'chemical_connectivity') return null;
    for (const raw of [getEntitySecondaryName(entity), entity.canonicalIdentifier]) {
      const value = raw?.trim();
      if (!value || value === displayName || isHashLike(value)) continue;
      const type = getIdentifierDisplayTypeForValue(entity, value);
      if (type && isPlainNameIdentifierType(type)) continue;
      const label = type ? getIdentifierTypeLabel(type) : '';
      // "CHEBI:15954" already names its namespace.
      const prefixed = label && value.toLowerCase().startsWith(`${label.toLowerCase()}:`);
      return { label: prefixed ? '' : label, value };
    }
    return null;
  }

  function groupSummary(entity: EntityResult): string | undefined {
    const count = entity.groupMemberCount ?? 0;
    if (count < 2) return undefined;
    return entity.groupStrategy === 'gene_reference' ? `${count} records` : `${count} forms`;
  }
</script>

{#snippet card()}
  {@const publicId = getEntityPublicId(result)}
  {@const displayName = getEntityDisplayName(result)}
  {@const identifier = rowIdentifier(result, displayName)}
  {@const identifierText = identifier
    ? [identifier.label, identifier.value].filter(Boolean).join(' ')
    : ''}
  {@const entityTypeLabel = getEntityTypeLabel(result)}
  {@const selected = result.groupMemberKeys
    ? selection.selectedEntities.filter((entity) => entity.groupKey === publicId).length ===
        result.groupMemberCount ||
      (!result.groupMemberCursor &&
        result.groupMemberKeys.every((key) => selection.isSelected(key)))
    : selection.isSelected(publicId)}
  {@const typeKey = result.entityType || entityTypeLabel}
  {@const entityTypeStyle = getEntityTypeStyle(typeKey)}
  {@const unresolved = isUnresolvedEntity(result)}
  {@const showOntologyHint = !!getEntityOntologyTerm(result)}
  {@const meta = [
    entityTypeLabel,
    formatTaxonomy(result.taxonomyId, result.taxonomyName),
    groupSummary(result),
    chemicalHint(result),
    unresolved ? 'Unresolved' : null,
  ].filter(Boolean) as string[]}
  <div
    class={`group/row grid w-full min-w-0 cursor-pointer grid-cols-[2rem_minmax(0,1fr)_auto] items-center gap-x-3 px-3 py-2 transition-colors outline-none focus-visible:ring-2 focus-visible:ring-inset focus-visible:ring-ring ${
      unresolved
        ? 'bg-amber-500/[0.05] outline-1 -outline-offset-1 outline-dashed outline-amber-500/60 hover:bg-amber-500/[0.10]'
        : entityTypeStyle.rowClass
    } ${selected ? 'ring-1 ring-inset ring-green-600 dark:ring-green-500' : ''}`}
    role="button"
    tabindex="0"
    aria-label={`${displayName}, ${meta.join(', ')}`}
    onclick={() => openDetails(result)}
    onkeydown={(event) => handleDetailsKeydown(event, result)}
  >
    <!-- Same emoji as the entity type filter. -->
    <span
      class={`grid size-8 place-items-center rounded-lg text-base leading-none ${
        unresolved ? 'bg-amber-500/15' : entityTypeStyle.iconClass
      }`}
      aria-hidden="true">{getEntityTypeEmoji(typeKey)}</span
    >
    <div class="min-w-0">
      <div class="flex min-w-0 items-baseline gap-2">
        <span
          class="min-w-0 shrink truncate text-sm font-semibold text-foreground"
          use:truncateTitle={displayName}>{displayName}</span
        >
        {#if identifier}<span
            class="min-w-[3ch] flex-1 truncate font-mono text-[11px] text-muted-foreground/80 max-sm:hidden"
            use:truncateTitle={identifierText}>{identifierText}</span
          >{/if}
      </div>
      <p class="truncate text-xs text-muted-foreground" use:truncateTitle={meta.join(' · ')}>
        {meta.join(' · ')}
      </p>
    </div>
    <div class="flex shrink-0 items-center gap-2">
      {#if showOntologyHint}
        <span title="Has an ontology hierarchy" class="text-muted-foreground">
          <Network class="size-3.5" />
        </span>
      {/if}
      {#if (result.relationCount || 0) > 0}
        <span
          class="inline-flex items-center gap-1 text-xs text-muted-foreground tabular-nums"
          title="Relations"
        >
          <Link class="size-3.5" />
          <span class="font-semibold text-foreground"
            >{formatNumber(result.relationCount || 0)}</span
          >
        </span>
      {/if}
      <button
        type="button"
        disabled={selectingGroup}
        aria-label={selected ? 'Remove from selection' : 'Add to selection'}
        title={selected ? 'Remove from selection' : 'Add to selection'}
        onkeydown={(event) => event.stopPropagation()}
        onclick={(event) => {
          event.stopPropagation();
          if (result.groupMemberKeys) {
            void toggleGroup();
          } else if (selected) {
            selection.removeEntity(publicId);
          } else {
            selection.addEntity({
              id: publicId,
              entityId: publicId,
              entityPk: result.entityPk,
              name: displayName,
              type: entityTypeLabel,
              fullResult: result,
            });
          }
        }}
        class={`grid size-7 place-items-center rounded-md transition focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring ${
          selected
            ? 'bg-green-600 text-white hover:bg-green-700'
            : selectingGroup
              ? 'text-muted-foreground'
              : 'text-muted-foreground opacity-0 hover:bg-foreground/10 hover:text-foreground focus-visible:opacity-100 group-hover/row:opacity-100 group-focus-within/row:opacity-100 [@media(hover:none)]:opacity-100'
        }`}
      >
        {#if selectingGroup}
          <LoaderCircle class="size-4 animate-spin" />
        {:else if selected}
          <Check class="size-4" />
        {:else}
          <Plus class="size-4" />
        {/if}
      </button>
    </div>
  </div>
{/snippet}
{@render card()}
{#if selectionError}<p role="alert" class="text-sm text-destructive">{selectionError}</p>{/if}
{#if !onopen}<EntityDetailsDialog bind:open={detailsOpen} entity={result} />{/if}
