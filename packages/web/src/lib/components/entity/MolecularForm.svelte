<script lang="ts">
  import type { MolecularFormRecord } from '$lib/types/molecular';
  import {
    featureRange,
    productLabel,
    type ProductSummary,
  } from '$lib/utils/molecular-presentation';
  let {
    form,
    products = [],
    onSelectProduct,
  }: {
    form?: MolecularFormRecord | null;
    products?: readonly ProductSummary[];
    onSelectProduct?: (key: string) => void;
  } = $props();

  function identifierLabel(value: unknown): string | undefined {
    if (!value || typeof value !== 'object') return undefined;
    const identifier = value as Record<string, unknown>;
    if (typeof identifier.ns !== 'string' || typeof identifier.id !== 'string') return undefined;
    return `${identifier.ns}:${identifier.id}`;
  }

  /** A modification's name: its source description when the term is an accession. */
  function modificationLabel(item: Record<string, unknown>): string {
    const term = typeof item.term === 'string' ? item.term : '';
    const description = typeof item.description === 'string' ? item.description : '';
    if (/^MOD:\d+$/.test(term) && description) return `${description} · ${term}`;
    return term || description || 'Unspecified modification';
  }

  function coordinateLabel(value: unknown): string {
    if (!value || typeof value !== 'object') return 'Coordinates unspecified';
    const reference = value as Record<string, unknown>;
    const system = reference.coordinate_system;
    const labels = [
      typeof system === 'string' && system !== 'unknown'
        ? `${system} coordinates`
        : 'Coordinate system unspecified',
    ];
    if (reference.position_base === 0 || reference.position_base === 1) {
      labels.push(`${reference.position_base}-based positions`);
    }
    labels.push(identifierLabel(reference.identifier) || 'Sequence unspecified');
    return labels.join(' · ');
  }
</script>

{#if form}
  <dl class="space-y-2 text-sm">
    {#if form.protein_entity_key}
      {@const key = form.protein_entity_key}
      {@const label = productLabel(key, products)}
      <div>
        <dt class="text-muted-foreground">Protein reference</dt>
        <dd>
          {#if onSelectProduct}<button class="underline" onclick={() => onSelectProduct?.(key)}
              >{label || 'View protein'}</button
            >{:else}{label || 'Protein reference reported'}{/if}
        </dd>
      </div>{/if}
    {#if form.transcript_entity_key}
      {@const key = form.transcript_entity_key}
      {@const label = productLabel(key, products)}
      <div>
        <dt class="text-muted-foreground">Transcript reference</dt>
        <dd>
          {#if onSelectProduct}<button class="underline" onclick={() => onSelectProduct?.(key)}
              >{label || 'View transcript'}</button
            >{:else}{label || 'Transcript reference reported'}{/if}
        </dd>
      </div>{/if}
    {#if form.isoform_identifier}<div>
        <dt class="text-muted-foreground">Reported isoform</dt>
        <dd>{form.isoform_identifier.ns}:{form.isoform_identifier.id}</dd>
      </div>{/if}
    {#if form.sequence_identifiers?.length}<div>
        <dt class="text-muted-foreground">Sequence references</dt>
        <dd>{form.sequence_identifiers.map((i) => `${i.ns}:${i.id}`).join(', ')}</dd>
      </div>{/if}
    {#if form.modifications?.length}<div>
        <dt class="text-muted-foreground">Reported modifications</dt>
        {#each form.modifications as item}<dd class="break-words">
            {modificationLabel(item)}
            {String(item.residue || '')}{featureRange(
              item.position,
              item.end_position,
            )}{#if item.coordinate_reference}<span class="block text-xs text-muted-foreground"
                >{coordinateLabel(item.coordinate_reference)}</span
              >{/if}
          </dd>{/each}
      </div>{/if}
    {#if form.variants?.length}<div>
        <dt class="text-muted-foreground">Reported variants</dt>
        {#each form.variants as item}<dd class="break-words">
            {String(item.description || '')}
            {String(item.reference || '')}{featureRange(
              item.position,
              item.end_position,
            )}{item.alternate === ''
              ? ' deleted'
              : String(item.alternate || '')}{#if item.identifier}<span class="block text-xs"
                >{identifierLabel(item.identifier) || 'Variant identifier reported'}</span
              >{/if}{#if item.coordinate_reference}<span class="block text-xs text-muted-foreground"
                >{coordinateLabel(item.coordinate_reference)}</span
              >{/if}
          </dd>{/each}
      </div>{/if}
    {#if form.regions?.length}<div>
        <dt class="text-muted-foreground">Reported regions</dt>
        {#each form.regions as item}<dd class="break-words">
            {String(item.type || 'Unspecified region')}
            {featureRange(item.position, item.end_position)}{#if item.identifier}<span
                class="block text-xs">{identifierLabel(item.identifier)}</span
              >{/if}{#if item.description}<span class="block text-xs"
                >{String(item.description)}</span
              >{/if}{#if item.coordinate_reference}<span class="block text-xs text-muted-foreground"
                >{coordinateLabel(item.coordinate_reference)}</span
              >{/if}
          </dd>{/each}
      </div>{/if}
  </dl>
{:else}<p class="text-sm text-muted-foreground">Molecular form unspecified.</p>{/if}
