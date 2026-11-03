<script lang="ts">
  import { AlertTriangle } from '@lucide/svelte';
  import { getEntityTypeStyle } from '$lib/utils/entity-types';
  import { getEntityTypeIcon } from '$lib/utils/entity-type-icons';

  function isResolverBackedEntityType(entityType: string | undefined) {
    if (!entityType) return false;
    const typeName = entityType.includes(':') ? entityType.split(':')[0] : entityType;
    const normalized = typeName.toLowerCase().replace(/[\s_]/g, '');
    return (
      normalized === 'protein' ||
      normalized === 'chemical' ||
      normalized === 'smallmolecule' ||
      normalized === 'compound' ||
      normalized === 'metabolite' ||
      normalized === 'drug' ||
      normalized === 'lipid'
    );
  }

  interface Props {
    displayName: string;
    canonicalIdentifier: string;
    entityType?: string;
    resolutionStatus?: string | null;
  }

  let { displayName, canonicalIdentifier, entityType, resolutionStatus }: Props = $props();

  const typeConfig = $derived(getEntityTypeStyle(entityType));
  const TypeIcon = $derived(getEntityTypeIcon(entityType));
  const isDuplicate = $derived(
    displayName === canonicalIdentifier || /^\d+$/.test(canonicalIdentifier),
  );
  const isUnresolved = $derived(
    isResolverBackedEntityType(entityType) &&
      resolutionStatus?.trim().toLowerCase() === 'unresolved',
  );
  const rootClass = $derived(
    isUnresolved
      ? 'border-dashed border-amber-300/90 bg-amber-50/45 dark:border-amber-700/70 dark:bg-amber-950/20'
      : `bg-gradient-to-br ${typeConfig.bgColor} ${typeConfig.borderColor}`,
  );
</script>

<div class="relative">
  <div
    class="relative {rootClass} backdrop-blur-sm border rounded-md px-2 py-1 shadow-sm min-w-[80px] w-full"
  >
    <div class="flex items-center gap-1.5 min-h-[32px]">
      <TypeIcon class="h-4 w-4 {typeConfig.color} shrink-0" />
      <div class="flex flex-col justify-center flex-1 min-w-0">
        <span
          class="truncate text-xs font-medium leading-tight {isUnresolved
            ? 'text-slate-700 dark:text-slate-200'
            : 'text-slate-900 dark:text-slate-100'}"
          title={`${typeConfig.label}: ${displayName || canonicalIdentifier}`}
        >
          {displayName || canonicalIdentifier}
        </span>
        {#if displayName && !isDuplicate}
          <span
            class="text-[10px] font-mono text-slate-500 dark:text-slate-400 truncate leading-none"
            title={canonicalIdentifier}
          >
            {canonicalIdentifier}
          </span>
        {/if}
      </div>
      {#if isUnresolved}
        <span
          class="inline-flex shrink-0 items-center gap-1 rounded-md border border-amber-300/80 bg-amber-100/70 px-1.5 py-0.5 text-[10px] font-medium leading-none text-amber-800 dark:border-amber-700/70 dark:bg-amber-950/60 dark:text-amber-200"
          title="Unresolved entity"
        >
          <AlertTriangle class="size-3" />
          Unresolved
        </span>
      {/if}
    </div>
  </div>
</div>
