<script lang="ts">
  import MolecularForm from '$lib/components/entity/MolecularForm.svelte';
  import EntityBadge from '$lib/components/entity/EntityBadge.svelte';
  import RelationLabel from './RelationLabel.svelte';
  import EntityDetailsDialog from '$lib/components/entity/EntityDetailsDialog.svelte';
  import type { EntityLike } from '$lib/domain/display';
  import { Search, ExternalLink } from '@lucide/svelte';
  import {
    getEntityBadgeIdentifier,
    getEntityDisplayName,
    getEntityTypeLabel,
    getIdentifierTypeLabel,
  } from '$lib/domain/display';
  import type {
    InteractionDetailsData,
    InteractionListRow,
    ParsedAnnotation,
  } from '$lib/types/interactions';
  import { measurementPresentation } from '$lib/utils/measurements';
  import { groupPublications, isPublicationTerm } from '$lib/utils/publications';
  import {
    annotationDefaultUnit,
    annotationLabelFor,
    annotationValueText,
    collapseListedAnnotations,
    safeSourceUrl,
  } from '$lib/utils/annotation-presentation';
  import { observationSummary, type ProductSummary } from '$lib/utils/molecular-presentation';

  interface Props {
    selectedInteraction: InteractionDetailsData | InteractionListRow | null;
    evidenceLoading?: boolean;
    products?: readonly ProductSummary[];
    onSelectProduct?: (key: string) => void;
  }

  let {
    selectedInteraction,
    evidenceLoading = false,
    products = [],
    onSelectProduct,
  }: Props = $props();

  let entityDetailsOpen = $state(false);
  let detailsEntity = $state<EntityLike | null>(null);

  function openEntityDetails(entity: EntityLike) {
    detailsEntity = entity;
    entityDetailsOpen = true;
  }

  function parseAnnotationTerm(value: string): { term: string; termId?: string } {
    const text = value.trim();
    if (!text) return { term: '' };

    // 1. Direct or CURIE match in CV_TERM_LABELS dictionary (e.g. MI:0364 -> Inferred by Curator)
    const label = getIdentifierTypeLabel(text);
    if (
      label &&
      label.toLowerCase() !== text.toLowerCase() &&
      /^[a-z0-9_-]+:[a-z0-9_-]+$/i.test(text)
    ) {
      return {
        term: label,
        termId: text.toUpperCase(),
      };
    }

    const parts = text.split(':');

    // 2. New format: Label:PREFIX:NUMBER (e.g. Ncbi Tax Id:OM:0205)
    if (parts.length >= 3 && /^[A-Z][A-Z0-9_-]*$/.test(parts[parts.length - 2])) {
      return {
        term: parts.slice(0, -2).join(':').trim(),
        termId: `${parts[parts.length - 2]}:${parts[parts.length - 1]}`,
      };
    }

    // 3. Legacy format: PREFIX:NUMBER:Label (e.g., OM:1228:Source, MI:0840:positive)
    if (parts.length >= 3) {
      return {
        termId: `${parts[0]}:${parts[1]}`,
        term: parts.slice(2).join(':').trim(),
      };
    }

    // 4. 2-part CURIE: PREFIX:NUMBER (e.g., MI:0364, OM:0007)
    if (parts.length === 2 && /^[a-z0-9_-]+$/i.test(parts[0]) && /^[a-z0-9_-]+$/i.test(parts[1])) {
      const curieLabel = getIdentifierTypeLabel(text);
      if (curieLabel && curieLabel.toLowerCase() !== text.toLowerCase()) {
        return {
          term: curieLabel,
          termId: text.toUpperCase(),
        };
      }
      return {
        term: text,
      };
    }

    return { term: text };
  }

  function normalizeAnnotationArray(value: unknown): ParsedAnnotation[] {
    if (!Array.isArray(value)) return [];

    return value.flatMap((item) => {
      if (!item || typeof item !== 'object') return [];
      const record = item as Record<string, unknown>;
      if (typeof record.term !== 'string' || !record.term.trim()) return [];

      const term = parseAnnotationTerm(record.term);
      const unit = typeof record.unit === 'string' ? parseAnnotationTerm(record.unit) : undefined;

      const measurement = measurementPresentation(record);
      return [
        {
          term: term.term,
          termId: term.termId,
          value:
            measurement?.value ?? (typeof record.value === 'string' ? record.value : undefined),
          sourceField: measurement?.sourceField,
          unit: measurement?.unit ?? unit?.term ?? annotationDefaultUnit(record.term),
          unitId: unit?.termId,
        },
      ];
    });
  }

  function dedupeAnnotations(values: ParsedAnnotation[]): ParsedAnnotation[] {
    const seen = new Set<string>();
    const deduped: ParsedAnnotation[] = [];

    for (const value of values) {
      const key = [
        value.term,
        value.termId,
        value.sourceField,
        value.value,
        value.unit,
        value.unitId,
      ].join('|');
      if (seen.has(key)) continue;
      seen.add(key);
      deduped.push(value);
    }

    return collapseListedAnnotations(deduped).sort((a, b) => {
      const byTerm = a.term.localeCompare(b.term);
      if (byTerm !== 0) return byTerm;
      return (a.value || '').localeCompare(b.value || '');
    });
  }

  function splitPubmedAnnotations(annotations: ParsedAnnotation[]): {
    pubmedIds: string[];
    annotations: ParsedAnnotation[];
  } {
    const publications = groupPublications(
      annotations.map((a) => ({ term: a.term, value: a.value || '' })),
    );
    return {
      // Historical field name retained inside this component; values are CURIEs.
      pubmedIds: publications.map((p) => p.id),
      annotations: annotations.filter((a) => !isPublicationTerm(a.term)),
    };
  }

  function extractComments(annotations: ParsedAnnotation[]): {
    comments: string[];
    nonComments: ParsedAnnotation[];
  } {
    const comments: string[] = [];
    const nonComments: ParsedAnnotation[] = [];
    for (const ann of annotations) {
      if (ann.term.toLowerCase() === 'comment' && ann.value && ann.value !== 'None') {
        comments.push(ann.value);
      } else {
        nonComments.push(ann);
      }
    }
    return { comments, nonComments };
  }

  function buildEvidenceRows(evidence: InteractionDetailsData['evidence']) {
    return evidence.map((row, index) => {
      const subjectSplit = splitPubmedAnnotations(normalizeAnnotationArray(row.subjectAttributes));
      const relationSplit = splitPubmedAnnotations([
        ...normalizeAnnotationArray(row.recordAttributes),
        ...normalizeAnnotationArray(row.evidence),
      ]);
      const objectSplit = splitPubmedAnnotations(normalizeAnnotationArray(row.objectAttributes));

      const { comments: subjectComments, nonComments: subjectAnns } = extractComments(
        subjectSplit.annotations,
      );
      const { comments: relationComments, nonComments: relationAnns } = extractComments(
        relationSplit.annotations,
      );
      const { comments: objectComments, nonComments: objectAnns } = extractComments(
        objectSplit.annotations,
      );
      const allComments = Array.from(
        new Set([...subjectComments, ...relationComments, ...objectComments]),
      );

      return {
        key: `${row.source || 'unknown'}-${index}`,
        source: row.source,
        subjectMolecularForm: row.subjectMolecularForm,
        objectMolecularForm: row.objectMolecularForm,
        pubmedIds: Array.from(
          new Set([
            ...subjectSplit.pubmedIds,
            ...relationSplit.pubmedIds,
            ...objectSplit.pubmedIds,
          ]),
        ).sort((a, b) => a.localeCompare(b, undefined, { numeric: true })),
        subjectAnnotations: dedupeAnnotations(subjectAnns),
        relationAnnotations: dedupeAnnotations(relationAnns),
        objectAnnotations: dedupeAnnotations(objectAnns),
        comments: allComments,
      };
    });
  }

  const detailedInteraction = $derived(
    selectedInteraction && 'evidence' in selectedInteraction ? selectedInteraction : null,
  );
  const evidence = $derived(detailedInteraction?.evidence ?? []);
  const evidenceRows = $derived(buildEvidenceRows(evidence));
  const subjectEntity = $derived(selectedInteraction?.subjectEntity);
  const objectEntity = $derived(selectedInteraction?.objectEntity);
  const loadedSummary = $derived(
    selectedInteraction
      ? observationSummary(
          selectedInteraction.relation,
          !evidenceLoading && 'evidence' in selectedInteraction
            ? selectedInteraction.evidence
            : undefined,
        )
      : undefined,
  );
  const knownProducts = $derived([
    ...products,
    ...(subjectEntity ? [subjectEntity] : []),
    ...(objectEntity ? [objectEntity] : []),
  ]);
</script>

{#snippet annotationRows(title: string, annotations: ParsedAnnotation[])}
  {#if annotations.length}
    <section class="min-w-0 space-y-2">
      <h4 class="text-sm font-medium">{title}</h4>
      <dl class="divide-y border-y">
        {#each annotations as annotation}
          <div class="grid gap-1 py-2.5 sm:grid-cols-[minmax(0,1fr)_minmax(0,2fr)] sm:gap-6">
            <dt class="break-words text-sm text-muted-foreground">
              {annotationLabelFor(annotation.sourceField || annotation.term, annotation.value)}
              {#if annotation.termId && annotation.termId.toLowerCase() !== annotation.term.toLowerCase()}<span
                  class="mt-0.5 block font-mono text-[11px]">{annotation.termId}</span
                >{/if}
            </dt>
            <dd class="min-w-0 break-words text-sm leading-6 [overflow-wrap:anywhere]">
              {#if safeSourceUrl(annotation.value)}
                <a
                  href={safeSourceUrl(annotation.value)}
                  target="_blank"
                  rel="noreferrer"
                  class="text-foreground underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                  >{annotation.value} <ExternalLink class="inline size-3" /></a
                >
              {:else if annotation.value && annotation.value !== 'None' && annotation.value !== 'null'}
                {annotationValueText(annotation.term, annotation.value)}{annotation.unit
                  ? ` ${annotation.unit}`
                  : ''}
              {:else}<span class="text-muted-foreground">Recorded without a value</span>{/if}
            </dd>
          </div>
        {/each}
      </dl>
    </section>
  {/if}
{/snippet}

{#if !selectedInteraction || !subjectEntity || !objectEntity}
  <div class="p-8 text-center">
    <Search class="mx-auto mb-4 size-8 text-muted-foreground" />
    <p class="text-sm text-muted-foreground">Select a relation to view details.</p>
  </div>
{:else}
  <div class="space-y-7 p-6">
    <section aria-label="Interaction summary" class="space-y-5">
      <div class="grid gap-4 sm:grid-cols-[minmax(0,1fr)_auto_minmax(0,1fr)] sm:items-center">
        <button
          type="button"
          class="min-w-0 rounded-md text-left transition-opacity hover:opacity-80 focus-visible:outline-2 focus-visible:outline-offset-4 focus-visible:outline-ring"
          aria-label={`View subject: ${getEntityDisplayName(subjectEntity)}`}
          onclick={() => openEntityDetails(subjectEntity)}
        >
          <EntityBadge
            displayName={getEntityDisplayName(subjectEntity)}
            canonicalIdentifier={getEntityBadgeIdentifier(subjectEntity)}
            entityType={getEntityTypeLabel(subjectEntity)}
            resolutionStatus={subjectEntity.resolutionStatus}
          />
        </button>
        <div class="min-w-0 space-y-1 px-3 py-2 text-center sm:max-w-52">
          <p class="text-sm font-medium leading-snug">
            <RelationLabel relation={selectedInteraction.relation} />
          </p>
          {#if selectedInteraction.relation.isDirected}<p class="text-xs text-muted-foreground">
              Directed
            </p>{/if}
        </div>
        <button
          type="button"
          class="min-w-0 rounded-md text-left transition-opacity hover:opacity-80 focus-visible:outline-2 focus-visible:outline-offset-4 focus-visible:outline-ring"
          aria-label={`View object: ${getEntityDisplayName(objectEntity)}`}
          onclick={() => openEntityDetails(objectEntity)}
        >
          <EntityBadge
            displayName={getEntityDisplayName(objectEntity)}
            canonicalIdentifier={getEntityBadgeIdentifier(objectEntity)}
            entityType={getEntityTypeLabel(objectEntity)}
            resolutionStatus={objectEntity.resolutionStatus}
          />
        </button>
      </div>
      <dl class="grid gap-4 border-y py-4 text-sm sm:grid-cols-3">
        <div>
          <dt class="text-xs text-muted-foreground">Source observations</dt>
          <dd class="mt-1 tabular-nums">
            {loadedSummary?.evidenceCount.toLocaleString()}
          </dd>
        </div>
        <div>
          <dt class="text-xs text-muted-foreground">Sources</dt>
          <dd class="mt-1 break-words">
            {loadedSummary?.sources.join(', ') || 'Unknown'}
          </dd>
        </div>
        <div>
          <dt class="text-xs text-muted-foreground">Relation category</dt>
          <dd class="mt-1 break-words">{selectedInteraction.relation.relationCategory}</dd>
        </div>
      </dl>
    </section>
    <section class="space-y-5" aria-label="Evidence">
      <div>
        <h3 class="font-semibold">Evidence</h3>
      </div>
      {#if evidenceLoading}
        <p role="status" class="py-6 text-sm text-muted-foreground">Loading evidence…</p>
      {:else if !evidenceRows.length}
        <p class="py-6 text-sm text-muted-foreground">No evidence is available.</p>
      {:else}
        {#each evidenceRows as row, index (row.key)}
          <article class="space-y-5 rounded-lg border p-4 sm:p-5">
            <header class="flex flex-wrap items-baseline justify-between gap-2">
              <h4 class="font-semibold">{row.source || 'Unknown source'}</h4>
              <span class="text-xs tabular-nums text-muted-foreground"
                >Observation {index + 1} of {evidenceRows.length}</span
              >
            </header>
            <div class="grid gap-4 sm:grid-cols-2">
              <section>
                <h4 class="mb-2 text-sm font-medium">Subject molecular form</h4>
                <MolecularForm
                  form={row.subjectMolecularForm}
                  products={knownProducts}
                  {onSelectProduct}
                />
              </section>
              <section>
                <h4 class="mb-2 text-sm font-medium">Object molecular form</h4>
                <MolecularForm
                  form={row.objectMolecularForm}
                  products={knownProducts}
                  {onSelectProduct}
                />
              </section>
            </div>
            {@render annotationRows('Relation measurements and context', row.relationAnnotations)}
            {#if row.subjectAnnotations.length || row.objectAnnotations.length}
              <div class="space-y-5">
                {@render annotationRows('Subject context', row.subjectAnnotations)}
                {@render annotationRows('Object context', row.objectAnnotations)}
              </div>
            {/if}
            {#if row.comments.length}
              <section class="space-y-2">
                <h4 class="text-sm font-medium">Notes</h4>
                {#each row.comments as comment}<p
                    class="break-words border-l-2 pl-3 text-sm leading-6"
                  >
                    {comment}
                  </p>{/each}
              </section>
            {/if}
            {#if row.pubmedIds.length}
              <section class="space-y-2">
                <h4 class="text-sm font-medium">Publications</h4>
                <ul class="flex flex-wrap gap-x-5 gap-y-2 text-sm">
                  {#each row.pubmedIds as pubmedId}
                    {@const publicationUrl = groupPublications([
                      { term: 'publications', value: pubmedId },
                    ])[0]?.url}
                    <li class="min-w-0 break-all">
                      {#if publicationUrl}<a
                          href={publicationUrl}
                          target="_blank"
                          rel="noreferrer"
                          class="text-foreground underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                          >{pubmedId} <ExternalLink class="inline size-3" /></a
                        >{:else}{pubmedId}{/if}
                    </li>
                  {/each}
                </ul>
              </section>
            {/if}
          </article>
        {/each}
      {/if}
    </section>
  </div>
{/if}

<EntityDetailsDialog bind:open={entityDetailsOpen} entity={detailsEntity} />
