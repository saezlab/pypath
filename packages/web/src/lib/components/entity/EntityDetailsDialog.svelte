<script lang="ts">
  import EntityRelationsTab from './EntityRelationsTab.svelte';
  import IdentifierBadge from './IdentifierBadge.svelte';
  import { entityFromWire } from '$lib/api/adapters';
  import { fetchEntitiesByPublicIds, fetchRelationsSearch } from '$lib/api/client';
  import { untrack } from 'svelte';
  import { page } from '$app/state';
  import { formatTaxonomy } from '$lib/utils/taxonomy';
  import { measurementPresentation } from '$lib/utils/measurements';
  import {
    annotationDefaultUnit,
    annotationLabelFor,
    annotationValueText,
    plainText,
    isNarrativeAnnotation,
    ontologyTermHref,
  } from '$lib/utils/annotation-presentation';
  import { releaseFetch } from '$lib/api/release';
  import { groupPublications, isPublicationTerm } from '$lib/utils/publications';
  import { getEntityTypeEmoji, getEntityTypeStyle } from '$lib/utils/entity-types';
  import { formatNumber } from '$lib/utils/format';
  import {
    ALL_SCOPE,
    identifiersOrgHref,
    keyIdentifiers,
    relationScopes,
    type MolecularContextLike,
    type RelationScope,
    type ScopeEntity,
  } from '$lib/utils/entity-overview';

  import { ArrowLeft, ExternalLink } from '@lucide/svelte';
  import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogHeader,
    DialogTitle,
  } from '$lib/components/ui/dialog/index.js';
  import { Button } from '$lib/components/ui/button/index.js';
  import { Input } from '$lib/components/ui/input/index.js';
  import * as Table from '$lib/components/ui/table/index.js';
  import * as Tabs from '$lib/components/ui/tabs/index.js';
  import MoleculeStructure from '$lib/components/entity/MoleculeStructure.svelte';
  import OntologyHierarchyBrowser from '$lib/components/explore/OntologyHierarchyBrowser.svelte';
  import {
    getEntityDisplayName,
    getEntityIdentifiers,
    getEntityPrimaryIdentifierBadge,
    getEntityPublicId,
    getEntitySecondaryName,
    getEntitySmiles,
    getEntityTypeLabel,
    getEntityTypeValue,
    getIdentifierTypeLabel,
    isChemicalEntity,
    isCvTermEntity,
    type EntityLike,
  } from '$lib/domain/display';

  interface Props {
    open: boolean;
    entity: EntityLike | null;
  }

  type AnnotationRow = {
    identity: string;
    term: string;
    label: string;
    value: string;
    unit: string;
    source: string;
  };

  type EntityIdentifierLike = {
    key: string;
    value: string;
  };

  type IdentifierRow = {
    section: 'Name' | 'Synonym' | 'Identifier';
    type: string;
    value: string;
    href: string | null;
  };

  type IdentifierGroup = {
    type: string;
    rows: IdentifierRow[];
  };

  type AnnotationGroup = {
    label: string;
    rows: AnnotationRow[];
  };

  // Collections a tab pages on its own (the API's detail_field).
  type DetailField = 'identifiers' | 'entityAttributes';

  type EntityOntologyHierarchyLike = {
    termId: string;
    ontologyPrefix: string | null;
    label: string | null;
    definition: string | null;
    ontologyId?: string | null;
    childCount?: number;
    parentCount?: number;
  };

  let { open = $bindable(false), entity }: Props = $props();
  let hydratedEntity = $state<EntityLike | null>(null);
  let loadingDetails = $state(false);
  let detailsError = $state<string | null>(null);

  let loadingField = $state<DetailField | null>(null);
  let identifierFilter = $state('');
  let annotationFilter = $state('');
  // Groups whose values are all shown, by tab and group label.
  let expandedGroups = $state<Record<string, boolean>>({});
  let descriptionExpanded = $state(false);
  let relationScopeList = $state<RelationScope[]>([]);
  let relationScopeId = $state(ALL_SCOPE);
  let relationTotals = $state<Record<string, number>>({});
  // Entities visited from the relations table, so the header can offer a way back.
  let history = $state<Array<{ entity: EntityLike; tab: string }>>([]);
  let restoreTab: string | null = null;
  let detailAttempt = $state(0);
  const detailKey = $derived(
    JSON.stringify([
      open,
      entity ? getEntityPublicId(entity) : null,
      page.data.selectedRelease,
      entity?.groupQuery,
      entity?.groupFilters,
      detailAttempt,
    ]),
  );
  let activeTab = $state('overview');

  $effect(() => {
    if (!open || !entity || loadingDetails) return;
    const current = hydratedEntity ?? entity;
    const rows = normalizeAnnotationRows(
      Array.isArray(current.entityAttributes) ? current.entityAttributes : [],
    );
    const available = [
      'overview',
      'relations',
      ...(getEntityIdentifierTotal(current) || normalizeIdentifierEntries(current).length
        ? ['identifiers']
        : []),
      ...(rows.some((row) => !isPublicationTerm(row.term) && !isNarrative(row.term)) ||
      nextCursor(current, 'entityAttributes')
        ? ['annotations']
        : []),
      ...(groupPublications(rows).length ? ['publications'] : []),
      ...(getOntologyHierarchy(current) ? ['ontology'] : []),
    ];
    if (available.length && !available.includes(activeTab)) activeTab = available[0];
  });

  // Descriptions lead the attribute pages; an annotations tab opened before any
  // annotation arrived keeps loading pages until one does.
  $effect(() => {
    if (activeTab !== 'annotations' || !hydratedEntity || loadingField || detailsError) return;
    const rows = normalizeAnnotationRows(
      Array.isArray(hydratedEntity.entityAttributes) ? hydratedEntity.entityAttributes : [],
    );
    if (rows.some((row) => !isPublicationTerm(row.term) && !isNarrative(row.term))) return;
    if (nextCursor(hydratedEntity, 'entityAttributes')) void loadMore('entityAttributes');
  });

  const FIRST_PAGE_SIZE = 50;
  const NEXT_PAGE_SIZE = 200;
  // Values shown per identifier type or annotation before "+N more".
  const GROUP_PREVIEW = 8;
  const DESCRIPTION_PREVIEW = 2;

  function normalizeIdentifierEntries(
    entityLike: EntityLike | null | undefined,
  ): EntityIdentifierLike[] {
    if (!entityLike) return [];
    const identifiers = new Map<string, EntityIdentifierLike>();
    for (const identifier of getEntityIdentifiers(entityLike)) {
      const key = identifier.key?.trim();
      const value = identifier.value?.trim();
      if (!key || !value) continue;
      identifiers.set(`${key.toLowerCase()}\u0000${value.toLowerCase()}`, { key, value });
    }
    return Array.from(identifiers.values());
  }

  // The first detail page, with identifiers and sources the search result already had.
  function mergeEntityDetails(fallbackEntity: EntityLike, nextEntity: EntityLike): EntityLike {
    const merged = { ...fallbackEntity, ...nextEntity } as EntityLike & {
      identifiers?: EntityIdentifierLike[];
      sources?: string[];
    };
    const identifiers = new Map<string, EntityIdentifierLike>();
    for (const entityLike of [fallbackEntity, nextEntity]) {
      for (const identifier of normalizeIdentifierEntries(entityLike)) {
        identifiers.set(
          `${identifier.key.toLowerCase()}\u0000${identifier.value.toLowerCase()}`,
          identifier,
        );
      }
    }
    if (identifiers.size > 0) merged.identifiers = Array.from(identifiers.values());
    if (!Array.isArray(nextEntity.entityAttributes))
      merged.entityAttributes = fallbackEntity.entityAttributes;
    const sources = [fallbackEntity, nextEntity].flatMap((entityLike) =>
      Array.isArray(entityLike.sources) ? entityLike.sources.filter(Boolean) : [],
    );
    if (sources.length > 0) merged.sources = Array.from(new Set(sources));
    return merged;
  }

  async function fetchDetailPage(
    fallback: EntityLike,
    signal?: AbortSignal,
    next?: { field: DetailField; offset: number },
  ) {
    const publicId = getEntityPublicId(fallback);
    const detail = {
      detail_limit: next ? NEXT_PAGE_SIZE : FIRST_PAGE_SIZE,
      detail_offset: next?.offset ?? 0,
      ...(next ? { detail_field: next.field } : {}),
    };
    const response = fallback.groupMemberKeys
      ? await releaseFetch('/app-api/entities/groups', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          signal,
          body: JSON.stringify({
            strategy: fallback.groupStrategy || 'chemical_connectivity',
            group_key: publicId,
            query: fallback.groupQuery,
            filters: fallback.groupFilters,
            resources: fallback.groupResources,
            include_details: true,
            ...detail,
          }),
        })
      : await releaseFetch(
          `/app-api/entities/${encodeURIComponent(publicId)}?${new URLSearchParams({
            includeRelationships: 'false',
            ...Object.fromEntries(Object.entries(detail).map(([k, v]) => [k, String(v)])),
          })}`,
          { signal },
        );
    if (!response.ok) throw new Error('Entity details could not be loaded.');
    const body = await response.json();
    const page = fallback.groupMemberKeys ? body.groups?.[0]?.entity : body.entity;
    if (page?.entityPk !== publicId)
      throw new Error('The response did not match the requested entity.');
    return page as Record<string, unknown>;
  }

  function nextCursor(entityLike: EntityLike | null, field: DetailField): string | null {
    const cursor = (entityLike as Record<string, unknown> | null)?.[`${field}NextCursor`];
    return typeof cursor === 'string' ? cursor : null;
  }

  function fieldTotal(entityLike: EntityLike, field: DetailField): number {
    const total = Number((entityLike as Record<string, unknown>)[`${field}Total`]);
    const rows = (entityLike as Record<string, unknown>)[field];
    return Number.isFinite(total) ? total : Array.isArray(rows) ? rows.length : 0;
  }

  async function loadMore(field: DetailField) {
    const cursor = nextCursor(hydratedEntity, field);
    if (!entity || loadingField || !cursor) return;
    const key = detailKey;
    loadingField = field;
    detailsError = null;
    try {
      const page = await fetchDetailPage(entity, undefined, { field, offset: Number(cursor) });
      if (key !== detailKey || !hydratedEntity) return;
      const current = hydratedEntity as Record<string, unknown>;
      const rows = page[field];
      hydratedEntity = {
        ...hydratedEntity,
        [field]: [
          ...(Array.isArray(current[field]) ? (current[field] as unknown[]) : []),
          ...(Array.isArray(rows) ? rows : []),
        ],
        [`${field}NextCursor`]: page[`${field}NextCursor`] ?? null,
      } as EntityLike;
    } catch (error) {
      if (key === detailKey) detailsError = (error as Error).message;
    } finally {
      if (key === detailKey) loadingField = null;
    }
  }

  function toggleGroup(id: string) {
    expandedGroups = { ...expandedGroups, [id]: !expandedGroups[id] };
  }

  function matches(filter: string, ...values: string[]) {
    const needle = filter.trim().toLowerCase();
    return !needle || values.some((value) => value.toLowerCase().includes(needle));
  }

  $effect(() => {
    const [isOpen, publicId] = JSON.parse(detailKey);
    const fallback = untrack(() => $state.snapshot(entity));
    hydratedEntity = null;
    detailsError = null;
    descriptionExpanded = false;
    activeTab = restoreTab ?? 'overview';
    restoreTab = null;
    loadingField = null;
    identifierFilter = '';
    annotationFilter = '';
    expandedGroups = {};
    loadingDetails = false;
    if (!isOpen || !fallback || !publicId) return;
    const controller = new AbortController();
    loadingDetails = true;
    untrack(() => fetchDetailPage(fallback, controller.signal))
      .then((next) => {
        if (!controller.signal.aborted)
          hydratedEntity = mergeEntityDetails(fallback, entityFromWire(next as never));
      })
      .catch((error) => {
        if (!controller.signal.aborted) detailsError = error.message;
      })
      .finally(() => {
        if (!controller.signal.aborted) loadingDetails = false;
      });
    return () => controller.abort();
  });

  function isGeneLike(entityLike: EntityLike): boolean {
    const type = (getEntityTypeValue(entityLike) || '').toLowerCase();
    return (
      entityLike.groupStrategy === 'gene_reference' ||
      Boolean(entityLike.geneReferenceKeys?.length) ||
      type === 'gene' ||
      type === 'protein'
    );
  }

  async function loadRelationScopes(current: EntityLike, signal: AbortSignal) {
    const scopeEntity: ScopeEntity = {
      entityPk: getEntityPublicId(current),
      label: getEntityDisplayName(current),
    };
    const memberKeys = current.groupMemberKeys ?? [];
    let context: MolecularContextLike | null = null;
    if (isGeneLike(current)) {
      const response = await releaseFetch(
        `/app-api/entities/${encodeURIComponent(scopeEntity.entityPk)}/molecular-context?view=reference&limit=1`,
        { signal },
      );
      if (response.ok) context = await response.json();
    }
    const members =
      !context?.referenceEntityKey && memberKeys.length > 1
        ? (await fetchEntitiesByPublicIds(memberKeys.slice(0, 24))).entities.map((member) => ({
            entityPk: member.entityPk,
            // The record label tells members apart ("Beta-D-Glucose") where display names agree.
            label: member.label || getEntityDisplayName(member),
            canonicalIdentifier: member.canonicalIdentifier,
          }))
        : [];
    return relationScopes(scopeEntity, { context, memberKeys, members });
  }

  // Structure groups only list their members once details load; gene groups scope by reference.
  const scopeKey = $derived(
    JSON.stringify([
      detailKey,
      entity && !isGeneLike(entity) ? ((hydratedEntity ?? entity).groupMemberKeys ?? null) : null,
    ]),
  );

  $effect(() => {
    const [key, loadedMemberKeys] = JSON.parse(scopeKey);
    const [isOpen, publicId] = JSON.parse(key);
    const snapshot = untrack(() => $state.snapshot(entity)) as EntityLike | null;
    const current =
      snapshot && loadedMemberKeys ? { ...snapshot, groupMemberKeys: loadedMemberKeys } : snapshot;
    relationScopeList = [];
    relationTotals = {};
    relationScopeId = ALL_SCOPE;
    if (!isOpen || !current || !publicId) return;
    if (current.groupMemberCount && !current.groupMemberKeys?.length) return;
    const controller = new AbortController();
    const countScopes = (scopes: RelationScope[]) => {
      for (const scope of scopes.slice(0, 25)) {
        fetchRelationsSearch({ filters: scope.filters, limit: 1, offset: 0 }, controller.signal)
          .then((data) => {
            if (!controller.signal.aborted)
              relationTotals = { ...relationTotals, [scope.id]: Number(data.total ?? 0) };
          })
          .catch(() => {});
      }
    };
    loadRelationScopes(current, controller.signal)
      .catch(() =>
        relationScopes(
          { entityPk: publicId, label: getEntityDisplayName(current) },
          { memberKeys: current.groupMemberKeys ?? [] },
        ),
      )
      .then(({ scopes, initial }) => {
        if (controller.signal.aborted) return;
        relationScopeList = scopes;
        relationScopeId = initial;
        countScopes(scopes);
      });
    return () => controller.abort();
  });

  $effect(() => {
    if (!open) history = [];
  });

  function showEntity(next: EntityLike) {
    if (entity) history = [...history, { entity, tab: activeTab }];
    entity = next;
  }

  function goBack() {
    const previous = history.at(-1);
    if (!previous) return;
    history = history.slice(0, -1);
    restoreTab = previous.tab;
    entity = previous.entity;
  }

  function showRelations(scopeId = ALL_SCOPE) {
    relationScopeId = scopeId;
    activeTab = 'relations';
  }

  function groupLabel(entityLike: EntityLike): string | null {
    const count = entityLike.groupMemberCount;
    if (!count) return null;
    const kind = entityLike.groupStrategy === 'gene_reference' ? 'Gene group' : 'Structure group';
    const types = entityLike.memberEntityTypes?.length
      ? ` (${entityLike.memberEntityTypes.join(', ')})`
      : '';
    return `${kind} · ${count} ${count === 1 ? 'entity' : 'entities'}${types}`;
  }

  function isNarrative(term: string) {
    return isNarrativeAnnotation(term);
  }

  function getDescriptionSections(entity: EntityLike) {
    const sections = new Map<
      string,
      { label: string; items: { text: string; source: string }[] }
    >();
    for (const row of normalizeAnnotationRows(
      Array.isArray(entity.entityAttributes) ? entity.entityAttributes : [],
    )) {
      if (!isNarrative(row.term) || !row.value) continue;
      const section = sections.get(row.label) ?? { label: row.label, items: [] };
      section.items.push({ text: plainText(row.value), source: row.source });
      sections.set(row.label, section);
    }
    return Array.from(sections.values());
  }

  function normalizeAnnotationRows(attributes: unknown[]): AnnotationRow[] {
    const rows = attributes.flatMap((attribute) => {
      if (typeof attribute !== 'object' || attribute === null || !('term' in attribute)) return [];
      const row = attribute as Record<string, unknown>;
      const measurement = measurementPresentation(row);
      const term = typeof row.term === 'string' ? row.term.trim() : '';
      const value = measurement?.value ?? (typeof row.value === 'string' ? row.value.trim() : '');
      const unit =
        measurement?.unit ??
        (typeof row.unit === 'string' && row.unit.trim()
          ? row.unit.trim()
          : (annotationDefaultUnit(term) ?? ''));
      const source = typeof row.source === 'string' ? row.source.trim() : '';
      const label = annotationLabelFor(measurement?.sourceField || term, value);
      if (!term) return [];

      return [
        {
          identity: JSON.stringify([term, row.value, row.quantity, unit]),
          term,
          label,
          value,
          unit,
          source,
        },
      ];
    });

    // The same statement from several sources is one row listing its sources.
    const uniqueRows = new Map<string, AnnotationRow>();
    for (const row of rows) {
      const existing = uniqueRows.get(row.identity);
      if (!existing) uniqueRows.set(row.identity, row);
      else if (row.source && !existing.source.split(', ').includes(row.source))
        existing.source = existing.source ? `${existing.source}, ${row.source}` : row.source;
    }
    return Array.from(uniqueRows.values());
  }

  function annotationSortRank(labelValue: string): number {
    const label = labelValue.toLowerCase();
    if (label === 'function') return 0;
    if (label === 'subcellular location') return 1;
    if (label === 'disease involvement' || label === 'disease') return 2;
    return 10;
  }

  function sortAnnotations(rows: AnnotationRow[]): AnnotationRow[] {
    return [...rows].sort(
      (a, b) =>
        annotationSortRank(a.label) - annotationSortRank(b.label) ||
        a.label.localeCompare(b.label) ||
        a.value.localeCompare(b.value),
    );
  }

  function formatAnnotationValue(row: AnnotationRow, entityLike: EntityLike): string {
    const taxon = /^ncbitaxon:(\d+)$/i.exec(row.value)?.[1];
    if (taxon && row.term.replace(/^biolink:/, '') === 'in_taxon') {
      const own = String(entityLike.taxonomyId ?? '') === taxon;
      return formatTaxonomy(taxon, own ? entityLike.taxonomyName : null) ?? row.value;
    }
    const rawValue = row.value
      ? `${annotationValueText(row.term, row.value)}${row.unit ? ` ${row.unit}` : ''}`
      : '';
    if (!rawValue) return 'No value';

    const withoutPrefix = rawValue.replace(
      new RegExp(`^${row.label.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}:\\s*`, 'i'),
      '',
    );
    const cleaned = plainText(withoutPrefix).replace(/\s+/g, ' ').trim();
    return cleaned || rawValue;
  }

  function groupAnnotations(rows: AnnotationRow[]): AnnotationGroup[] {
    const groups = new Map<string, AnnotationGroup>();
    for (const row of rows) {
      const group = groups.get(row.label) ?? { label: row.label, rows: [] };
      group.rows.push(row);
      groups.set(row.label, group);
    }
    return Array.from(groups.values());
  }

  function getEntityIdentifierTotal(entity: EntityLike): number {
    const total = (entity as { identifiersTotal?: unknown }).identifiersTotal;
    const numericTotal = Number(total);
    return Number.isFinite(numericTotal) ? numericTotal : normalizeIdentifierEntries(entity).length;
  }

  function identifierTypeText(identifier: EntityIdentifierLike): string {
    return `${identifier.key} ${getIdentifierTypeLabel(identifier.key)}`
      .toLowerCase()
      .replace(/[_-]/g, ' ');
  }

  function identifierNameSection(identifier: EntityIdentifierLike): 'Name' | 'Synonym' | null {
    const text = identifierTypeText(identifier);
    if (text.includes('synonym')) return 'Synonym';
    if (text.includes('gene name primary')) return 'Name';
    if (text.includes('recommended name')) return 'Name';
    if (text.includes('entry name')) return 'Name';
    if (getIdentifierTypeLabel(identifier.key).toLowerCase() === 'name') return 'Name';
    return null;
  }

  function getIdentifierRows(identifiers: EntityIdentifierLike[]): IdentifierRow[] {
    const nameSynonymRows = new Map<string, IdentifierRow>();
    const identifierRows: IdentifierRow[] = [];

    for (const identifier of identifiers) {
      const section = identifierNameSection(identifier);
      const type = getIdentifierTypeLabel(identifier.key);
      if (section) {
        const key = `${section}\u0000${type.toLowerCase()}\u0000${identifier.value.toLowerCase()}`;
        nameSynonymRows.set(key, { section, type, value: identifier.value, href: null });
      } else {
        identifierRows.push({
          section: 'Identifier',
          type,
          value: identifier.value,
          href: identifiersOrgHref(identifierTypeText(identifier), identifier.value),
        });
      }
    }

    return [
      ...Array.from(nameSynonymRows.values()).sort(
        (a, b) =>
          a.section.localeCompare(b.section) ||
          a.type.localeCompare(b.type) ||
          a.value.localeCompare(b.value),
      ),
      ...identifierRows.sort(
        (a, b) => a.type.localeCompare(b.type) || a.value.localeCompare(b.value),
      ),
    ];
  }

  function identifierSectionRank(section: IdentifierRow['section']): number {
    if (section === 'Name') return 0;
    if (section === 'Synonym') return 1;
    return 2;
  }

  function getIdentifierGroups(rows: IdentifierRow[]): IdentifierGroup[] {
    const groups = new Map<string, IdentifierGroup>();
    for (const row of rows) {
      const group = groups.get(row.type) ?? { type: row.type, rows: [] };
      group.rows.push(row);
      groups.set(row.type, group);
    }

    return Array.from(groups.values())
      .map((group) => ({
        ...group,
        rows: group.rows.sort(
          (a, b) =>
            identifierSectionRank(a.section) - identifierSectionRank(b.section) ||
            a.value.localeCompare(b.value),
        ),
      }))
      .sort(
        (a, b) =>
          Math.min(...a.rows.map((row) => identifierSectionRank(row.section))) -
            Math.min(...b.rows.map((row) => identifierSectionRank(row.section))) ||
          a.type.localeCompare(b.type),
      );
  }

  function getOntologyHierarchy(entityLike: EntityLike): EntityOntologyHierarchyLike | null {
    const hierarchy = (entityLike as { ontologyHierarchy?: unknown }).ontologyHierarchy;
    if (hierarchy && typeof hierarchy === 'object') {
      const row = hierarchy as Partial<EntityOntologyHierarchyLike>;
      if (row.termId) {
        return {
          termId: row.termId,
          ontologyPrefix: row.ontologyPrefix ?? null,
          label: row.label ?? null,
          definition: row.definition ?? null,
          ontologyId: row.ontologyId ?? null,
          childCount: Number(row.childCount || 0),
          parentCount: Number(row.parentCount || 0),
        };
      }
    }
    if (!isCvTermEntity(entityLike)) return null;
    const id = entityLike.canonicalIdentifier || '';
    if (!id) return null;
    const prefix = id.includes(':') ? id.split(':')[0] : id.includes('-') ? id.split('-')[0] : null;
    return {
      termId: id,
      ontologyPrefix: prefix,
      label: getEntityDisplayName(entityLike),
      definition: null,
      ontologyId: null,
      childCount: 0,
      parentCount: 0,
    };
  }
</script>

{#snippet moreAttributes(detailEntity: EntityLike)}
  {#if nextCursor(detailEntity, 'entityAttributes')}<div
      class="mt-4 flex flex-wrap items-center gap-3"
    >
      <Button
        variant="outline"
        size="sm"
        disabled={loadingField !== null}
        onclick={() => loadMore('entityAttributes')}
        >{loadingField === 'entityAttributes' ? 'Loading…' : 'Load more'}</Button
      >
      <span class="text-xs text-muted-foreground tabular-nums"
        >{formatNumber(
          (Array.isArray(detailEntity.entityAttributes) ? detailEntity.entityAttributes : [])
            .length,
          { compact: false },
        )} of {formatNumber(fieldTotal(detailEntity, 'entityAttributes'), { compact: false })} source
        annotations loaded</span
      >
    </div>{/if}
{/snippet}

<Dialog bind:open>
  <DialogContent
    class="flex h-[min(720px,90dvh)] max-h-[90dvh] flex-col gap-0 overflow-hidden p-0 sm:max-w-5xl"
  >
    {#if entity}
      {@const detailEntity = hydratedEntity ?? entity}
      {@const displayName = getEntityDisplayName(detailEntity)}
      {@const detailSections = getDescriptionSections(detailEntity)}
      {@const detailIdentifiers = normalizeIdentifierEntries(detailEntity)}
      {@const detailIdentifierRows = getIdentifierRows(detailIdentifiers)}
      {@const detailIdentifierGroups = getIdentifierGroups(detailIdentifierRows)}
      {@const detailIdentifierTotal = getEntityIdentifierTotal(detailEntity)}
      {@const detailSmiles = getEntitySmiles(detailEntity)}
      {@const headerIdentifiers = keyIdentifiers(
        getEntityPrimaryIdentifierBadge(detailEntity),
        detailIdentifiers,
      )}
      {@const secondaryName = getEntitySecondaryName(detailEntity)}
      {@const subtitle =
        secondaryName &&
        secondaryName !== detailEntity.canonicalIdentifier &&
        !detailIdentifiers.some(
          (identifier) =>
            identifier.value === secondaryName && identifierNameSection(identifier) === null,
        )
          ? secondaryName
          : null}
      {@const typeLabel = getEntityTypeLabel(detailEntity)}
      {@const typeStyle = getEntityTypeStyle(getEntityTypeValue(detailEntity))}
      {@const detailGroupLabel = groupLabel(detailEntity)}
      {@const detailTaxonomy = formatTaxonomy(detailEntity.taxonomyId, detailEntity.taxonomyName)}
      {@const ontologyHierarchy = getOntologyHierarchy(detailEntity)}
      {@const showChemicalStructure =
        !detailEntity.groupMemberKeys && isChemicalEntity(detailEntity) && detailSmiles}
      {@const detailAttributes = Array.isArray(detailEntity.entityAttributes)
        ? detailEntity.entityAttributes
        : []}
      {@const detailAnnotationRows = normalizeAnnotationRows(detailAttributes)}
      {@const detailPublications = groupPublications(detailAnnotationRows)}
      {@const detailNonPubmedAnnotations = sortAnnotations(
        detailAnnotationRows.filter(
          (row) => !isPublicationTerm(row.term) && !isNarrative(row.term),
        ),
      )}
      {@const detailAnnotationGroups = groupAnnotations(detailNonPubmedAnnotations)}
      {@const moreAttributesPending = nextCursor(detailEntity, 'entityAttributes') !== null}
      {@const relationTotal = relationTotals[ALL_SCOPE]}
      {@const partScopes = relationScopeList.filter((scope) => scope.kind !== 'all')}
      {@const longDescription = detailSections.some((section) =>
        section.items.some((item) => item.text.length > 420),
      )}

      <DialogHeader class="shrink-0 gap-3 px-6 pt-5 pb-4 pr-14 text-left">
        {#if history.length}
          <button
            type="button"
            class="inline-flex w-fit items-center gap-1 rounded-sm text-xs text-muted-foreground hover:text-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
            onclick={goBack}
          >
            <ArrowLeft class="size-3.5" />
            Back to {getEntityDisplayName(history[history.length - 1].entity)}
          </button>
        {/if}
        <div class="space-y-1">
          <DialogTitle class="break-words text-2xl font-semibold leading-tight tracking-tight"
            >{displayName}</DialogTitle
          >
          {#if subtitle}<DialogDescription class="break-words">{subtitle}</DialogDescription>{/if}
        </div>
        <div class="flex flex-wrap items-center gap-1.5">
          <span
            class={`inline-flex items-center gap-1 rounded-md border px-2 py-0.5 text-xs font-medium ${typeStyle.chipClass}`}
            ><span aria-hidden="true"
              >{getEntityTypeEmoji(getEntityTypeValue(detailEntity) || '')}</span
            >{typeLabel}</span
          >
          {#if detailGroupLabel}<span
              class="rounded-md border px-2 py-0.5 text-xs font-medium text-muted-foreground"
              >{detailGroupLabel}</span
            >{/if}
          {#if detailTaxonomy}<span
              class="inline-flex items-center gap-1 rounded-md border px-2 py-0.5 text-xs font-medium text-muted-foreground"
              ><span aria-hidden="true">🌿</span>{detailTaxonomy}</span
            >{/if}
        </div>
        {#if headerIdentifiers.length}
          <div class="flex flex-wrap gap-1.5">
            {#each headerIdentifiers as identifier (`${identifier.key}:${identifier.value}`)}
              {@const href = identifiersOrgHref(
                `${identifier.key} ${getIdentifierTypeLabel(identifier.key)}`,
                identifier.value,
              )}
              {#if href}
                <a
                  {href}
                  target="_blank"
                  rel="noreferrer"
                  class="inline-flex max-w-full items-center gap-1 rounded-md focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring [&:hover>span]:border-foreground/40"
                  title="Open on identifiers.org"
                >
                  <IdentifierBadge
                    identifierType={identifier.key}
                    value={identifier.value}
                    variant="subtle"
                  />
                </a>
              {:else}
                <IdentifierBadge
                  identifierType={identifier.key}
                  value={identifier.value}
                  variant="subtle"
                />
              {/if}
            {/each}
          </div>
        {/if}
      </DialogHeader>
      {#if detailsError}
        <div
          role="alert"
          class="mx-6 mb-3 flex flex-wrap items-center gap-3 rounded-md border border-destructive/40 px-3 py-2 text-sm text-destructive"
        >
          <span>{detailsError}</span>
          <Button variant="outline" size="sm" onclick={() => detailAttempt++}>Try again</Button>
        </div>
      {/if}
      <Tabs.Root bind:value={activeTab} class="min-h-0 flex-1 gap-0">
        <div class="shrink-0 overflow-x-auto border-b px-6 pb-3">
          <Tabs.List class="justify-start" aria-label="Entity details">
            <Tabs.Trigger value="overview" class="flex-none px-3">Overview</Tabs.Trigger>
            <Tabs.Trigger value="relations" class="flex-none px-3"
              >Relations {#if relationTotal !== undefined}<span
                  class="text-xs text-muted-foreground tabular-nums"
                  >{formatNumber(relationTotal)}</span
                >{/if}</Tabs.Trigger
            >
            {#if detailIdentifierRows.length > 0 || detailIdentifierTotal > 0}<Tabs.Trigger
                value="identifiers"
                class="flex-none px-3"
                >Identifiers <span class="text-xs text-muted-foreground tabular-nums"
                  >{formatNumber(detailIdentifierTotal)}</span
                ></Tabs.Trigger
              >{/if}
            {#if detailNonPubmedAnnotations.length > 0 || moreAttributesPending}<Tabs.Trigger
                value="annotations"
                class="flex-none px-3"
                >Annotations <span class="text-xs text-muted-foreground tabular-nums"
                  >{detailNonPubmedAnnotations.length}{moreAttributesPending ? '+' : ''}</span
                ></Tabs.Trigger
              >{/if}
            {#if detailPublications.length > 0}<Tabs.Trigger
                value="publications"
                class="flex-none px-3"
                >Publications <span class="text-xs text-muted-foreground tabular-nums"
                  >{detailPublications.length}{moreAttributesPending ? '+' : ''}</span
                ></Tabs.Trigger
              >{/if}
            {#if ontologyHierarchy}<Tabs.Trigger value="ontology" class="flex-none px-3"
                >Ontology</Tabs.Trigger
              >{/if}
          </Tabs.List>
        </div>
        <Tabs.Content value="overview" class="min-h-0 overflow-y-auto overscroll-contain p-6">
          <div
            class={`grid gap-8 ${showChemicalStructure ? 'md:grid-cols-[minmax(0,3fr)_minmax(0,2fr)]' : 'max-w-3xl'}`}
          >
            <div class="min-w-0 space-y-6">
              {#if detailSections.length}
                {#each detailSections as section (section.label)}
                  {@const expanded = expandedGroups[`description:${section.label}`]}
                  <section class="space-y-1.5">
                    <h3 class="text-xs font-semibold uppercase tracking-wide text-muted-foreground">
                      {section.label}
                      {#if section.items.length > 1}<span
                          class="ml-1 font-normal normal-case tracking-normal tabular-nums"
                          >{section.items.length}</span
                        >{/if}
                    </h3>
                    <div class="space-y-3">
                      {#each expanded ? section.items : section.items.slice(0, DESCRIPTION_PREVIEW) as item}
                        <div>
                          <p
                            class={`break-words text-sm leading-6 ${descriptionExpanded ? '' : 'line-clamp-5'}`}
                          >
                            {item.text}
                          </p>
                          {#if item.source}<p class="text-xs text-muted-foreground">
                              Source: {item.source}
                            </p>{/if}
                        </div>
                      {/each}
                    </div>
                    {#if section.items.length > DESCRIPTION_PREVIEW}<Button
                        variant="link"
                        size="sm"
                        class="h-auto px-0 text-xs"
                        onclick={() => toggleGroup(`description:${section.label}`)}
                        >{expanded
                          ? 'Show fewer'
                          : `Show ${section.items.length - DESCRIPTION_PREVIEW} more`}</Button
                      >{/if}
                  </section>
                {/each}
                {#if longDescription}<Button
                    variant="link"
                    size="sm"
                    class="h-auto px-0"
                    onclick={() => (descriptionExpanded = !descriptionExpanded)}
                    >{descriptionExpanded ? 'Show less' : 'Show more'}</Button
                  >{/if}
              {:else if loadingDetails && !hydratedEntity}
                <p role="status" class="text-sm text-muted-foreground">Loading details…</p>
              {:else}
                <p class="text-sm text-muted-foreground">
                  No description is available for this entity.
                </p>
              {/if}

              {#if partScopes.length}
                <section class="space-y-2">
                  <h3 class="text-xs font-semibold uppercase tracking-wide text-muted-foreground">
                    {partScopes[0].kind === 'product' ? 'Products' : 'Grouped entities'}
                    <span class="font-normal normal-case tracking-normal tabular-nums"
                      >{partScopes.length}</span
                    >
                  </h3>
                  <div class="divide-y overflow-hidden rounded-lg border">
                    {#each partScopes as scope (scope.id)}
                      <button
                        type="button"
                        class="flex w-full items-center justify-between gap-3 px-3 py-2 text-left text-sm transition-colors hover:bg-muted/50 focus-visible:bg-muted/50 focus-visible:outline-none"
                        onclick={() => showRelations(scope.id)}
                      >
                        <span class="min-w-0">
                          <span class="font-medium">{scope.title ?? scope.label}</span>
                          {#if scope.identifier && scope.identifier !== scope.title}<span
                              class="ml-2 break-all font-mono text-xs text-muted-foreground"
                              >{scope.identifier}</span
                            >{/if}
                        </span>
                        <span class="shrink-0 text-xs text-muted-foreground tabular-nums">
                          {relationTotals[scope.id] !== undefined
                            ? `${formatNumber(relationTotals[scope.id])} relations`
                            : 'Relations'} →
                        </span>
                      </button>
                    {/each}
                  </div>
                </section>
              {/if}
            </div>

            {#if showChemicalStructure}
              <section
                class="flex min-w-0 justify-center self-start overflow-x-auto rounded-lg border p-4"
              >
                <MoleculeStructure
                  smiles={detailSmiles}
                  width={300}
                  height={220}
                  renderOnClick={false}
                />
              </section>
            {/if}
          </div>
        </Tabs.Content>
        <Tabs.Content value="relations" class="min-h-0 overflow-y-auto overscroll-contain p-6">
          {#if activeTab === 'relations'}
            {#if relationScopeList.length}
              <EntityRelationsTab
                scopes={relationScopeList}
                bind:scopeId={relationScopeId}
                totals={relationTotals}
                onEntitySelect={showEntity}
              />
            {:else}
              <p role="status" class="text-sm text-muted-foreground">Loading relations…</p>
            {/if}
          {/if}
        </Tabs.Content>
        <Tabs.Content value="identifiers" class="min-h-0 overflow-y-auto overscroll-contain p-6">
          {#if detailIdentifierRows.length}
            {@const groups = detailIdentifierGroups
              .map((group) => ({
                ...group,
                rows: group.rows.filter((row) => matches(identifierFilter, group.type, row.value)),
              }))
              .filter((group) => group.rows.length)}
            {#if detailIdentifierRows.length > GROUP_PREVIEW}<Input
                bind:value={identifierFilter}
                placeholder="Filter identifiers"
                aria-label="Filter identifiers"
                class="mb-4 h-8 max-w-xs"
              />{/if}
            {#if groups.length}
              <dl class="divide-y rounded-lg border text-sm">
                {#each groups as group (group.type)}
                  {@const expanded = expandedGroups[`id:${group.type}`] || identifierFilter.trim()}
                  {@const shown = expanded ? group.rows : group.rows.slice(0, GROUP_PREVIEW)}
                  <div
                    class="grid gap-x-4 gap-y-1.5 px-3 py-2.5 sm:grid-cols-[10rem_minmax(0,1fr)]"
                  >
                    <dt class="text-xs font-medium text-muted-foreground sm:pt-0.5">
                      {group.type}
                      <span class="ml-1 font-normal tabular-nums">{group.rows.length}</span>
                    </dt>
                    <dd class="flex min-w-0 flex-wrap items-center gap-1.5">
                      {#each shown as row}
                        {#if row.href}
                          <a
                            href={row.href}
                            target="_blank"
                            rel="noreferrer"
                            class={`inline-flex max-w-full items-center gap-1 rounded-md border px-1.5 py-0.5 text-xs text-foreground hover:bg-muted ${row.section === 'Identifier' ? 'break-all font-mono' : 'break-words'}`}
                          >
                            {row.value}
                            <ExternalLink class="size-3 shrink-0 text-muted-foreground" />
                          </a>
                        {:else}
                          <span
                            class={`max-w-full rounded-md bg-muted px-1.5 py-0.5 text-xs text-foreground ${row.section === 'Identifier' ? 'break-all font-mono' : 'break-words'}`}
                            >{row.value}</span
                          >
                        {/if}
                      {/each}
                      {#if group.rows.length > GROUP_PREVIEW && !identifierFilter.trim()}
                        <Button
                          variant="link"
                          size="sm"
                          class="h-auto px-1 py-0 text-xs"
                          onclick={() => toggleGroup(`id:${group.type}`)}
                          >{expanded
                            ? 'Show less'
                            : `+${group.rows.length - GROUP_PREVIEW} more`}</Button
                        >
                      {/if}
                    </dd>
                  </div>
                {/each}
              </dl>
            {:else}<p class="text-sm text-muted-foreground">
                No identifiers match the filter.
              </p>{/if}
          {:else}<p class="text-muted-foreground">No identifiers are available.</p>{/if}
          {#if nextCursor(detailEntity, 'identifiers')}<div
              class="mt-4 flex flex-wrap items-center gap-3"
            >
              <Button
                variant="outline"
                size="sm"
                disabled={loadingField !== null}
                onclick={() => loadMore('identifiers')}
                >{loadingField === 'identifiers' ? 'Loading…' : 'Load more identifiers'}</Button
              >
              <span class="text-xs text-muted-foreground tabular-nums"
                >{formatNumber(detailIdentifiers.length, { compact: false })} of {formatNumber(
                  detailIdentifierTotal,
                  { compact: false },
                )} loaded</span
              >
            </div>{/if}
        </Tabs.Content>
        <Tabs.Content value="annotations" class="min-h-0 overflow-y-auto overscroll-contain p-6">
          {#if detailAnnotationGroups.length}
            {@const groups = detailAnnotationGroups
              .map((group) => ({
                ...group,
                rows: group.rows.filter((row) =>
                  matches(
                    annotationFilter,
                    group.label,
                    formatAnnotationValue(row, detailEntity),
                    row.source,
                  ),
                ),
              }))
              .filter((group) => group.rows.length)}
            {#if detailNonPubmedAnnotations.length > GROUP_PREVIEW}<Input
                bind:value={annotationFilter}
                placeholder="Filter annotations"
                aria-label="Filter annotations"
                class="mb-4 h-8 max-w-xs"
              />{/if}
            {#if groups.length}
              <Table.Root>
                <Table.Header class="sticky top-0 z-10 bg-muted/80 backdrop-blur">
                  <Table.Row>
                    <Table.Head class="w-48">Annotation</Table.Head>
                    <Table.Head>Value</Table.Head>
                    <Table.Head class="w-36">Source</Table.Head>
                  </Table.Row>
                </Table.Header>
                <Table.Body>
                  {#each groups as group (group.label)}
                    {@const expanded =
                      expandedGroups[`annotation:${group.label}`] || annotationFilter.trim()}
                    {@const shown = expanded ? group.rows : group.rows.slice(0, GROUP_PREVIEW)}
                    {@const more = group.rows.length - shown.length}
                    {@const toggle = group.rows.length > GROUP_PREVIEW && !annotationFilter.trim()}
                    {#each shown as annotation, index}
                      <Table.Row class={index < shown.length - 1 || toggle ? 'border-b-0' : ''}>
                        <Table.Cell class="whitespace-normal align-top font-medium text-foreground">
                          {#if index === 0}{group.label}{#if group.rows.length > 1}<span
                                class="ml-1.5 text-xs font-normal text-muted-foreground tabular-nums"
                                >{group.rows.length}</span
                              >{/if}{/if}
                        </Table.Cell>
                        <Table.Cell
                          class="max-w-lg whitespace-normal break-words align-top text-foreground"
                        >
                          {#if annotation.term.replace(/^biolink:/, '') === 'has_biological_sequence'}
                            <details>
                              <summary>Sequence ({annotation.value.length} residues)</summary>
                              <pre
                                class="whitespace-pre-wrap break-all text-xs">{annotation.value}</pre>
                            </details>
                          {:else if ontologyTermHref(annotation.value)}<a
                              href={ontologyTermHref(annotation.value)}
                              target="_blank"
                              rel="noreferrer"
                              class="inline-flex items-center gap-1 font-mono text-xs underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                              >{annotation.value}<ExternalLink class="size-3 shrink-0" /></a
                            >
                          {:else}{formatAnnotationValue(annotation, detailEntity)}{/if}
                        </Table.Cell>
                        <Table.Cell class="whitespace-normal align-top text-muted-foreground">
                          {annotation.source || 'Unknown'}
                        </Table.Cell>
                      </Table.Row>
                    {/each}
                    {#if toggle}
                      <Table.Row class="hover:bg-transparent">
                        <Table.Cell></Table.Cell>
                        <Table.Cell colspan={2} class="pt-0">
                          <Button
                            variant="link"
                            size="sm"
                            class="h-auto px-0 py-0 text-xs"
                            onclick={() => toggleGroup(`annotation:${group.label}`)}
                            >{expanded ? 'Show less' : `Show ${more} more`}</Button
                          >
                        </Table.Cell>
                      </Table.Row>
                    {/if}
                  {/each}
                </Table.Body>
              </Table.Root>
            {:else}<p class="text-sm text-muted-foreground">
                No annotations match the filter.
              </p>{/if}
          {:else if moreAttributesPending}<p role="status" class="text-sm text-muted-foreground">
              Loading annotations…
            </p>
          {:else}<p class="text-muted-foreground">
              No annotations are available for this entity.
            </p>{/if}
          {@render moreAttributes(detailEntity)}
        </Tabs.Content>
        <Tabs.Content value="publications" class="min-h-0 overflow-y-auto overscroll-contain p-6">
          {#if detailPublications.length}
            <Table.Root>
              <Table.Header class="sticky top-0 z-10 bg-muted/80 backdrop-blur">
                <Table.Row>
                  <Table.Head class="w-40">Publication</Table.Head>
                  <Table.Head>Sources</Table.Head>
                </Table.Row>
              </Table.Header>
              <Table.Body>
                {#each detailPublications as publication}
                  {@const publicationSources = publication.sources}
                  <Table.Row>
                    <Table.Cell class="align-top">
                      {#if publication.url}
                        <a
                          href={publication.url}
                          target="_blank"
                          rel="noreferrer"
                          class="inline-flex items-center gap-1 font-mono text-foreground underline decoration-muted-foreground/50 underline-offset-4 hover:decoration-foreground"
                        >
                          {publication.id}
                          <ExternalLink class="size-3" />
                        </a>
                      {:else}{publication.id}{/if}
                    </Table.Cell>
                    <Table.Cell
                      class="whitespace-normal break-words align-top text-muted-foreground"
                    >
                      {publicationSources.length ? publicationSources.join(', ') : 'Unknown'}
                    </Table.Cell>
                  </Table.Row>
                {/each}
              </Table.Body>
            </Table.Root>
          {:else}<p class="text-muted-foreground">No publications are available.</p>{/if}
          {@render moreAttributes(detailEntity)}
        </Tabs.Content>
        {#if ontologyHierarchy}
          <Tabs.Content value="ontology" class="min-h-0 overflow-y-auto overscroll-contain p-6"
            ><OntologyHierarchyBrowser term={ontologyHierarchy} /></Tabs.Content
          >
        {/if}
      </Tabs.Root>
    {/if}
  </DialogContent>
</Dialog>
