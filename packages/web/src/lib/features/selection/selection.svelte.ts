import { page } from '$app/stores';
import { goto } from '$app/navigation';
import { browser } from '$app/environment';
import type { SearchFilters } from '$lib/types/search';
import {
  parseEntityIdsParam,
  parseFiltersParam,
  parseSelectionTab,
  parseSearchMode,
  parseSearchType,
  parseEntityWorkflow,
  serializeEntityIdsParam,
  serializeFiltersParam,
  type SearchMode,
  type SearchType,
  type SelectionTab,
  type EntityWorkflow,
} from '$lib/navigation/url-codecs';

import {
  isEntityCache,
  isAnnotationCache,
  entityPrimaryKeyFromSelection,
  selectionFromUrl,
  type SelectedEntity,
  type SelectedAnnotation,
} from './model';
import { createUrlUpdater } from './navigation';
import { readStored, writeStored, isStringList } from './persistence';
export type { SelectedEntity, SelectedAnnotation } from './model';

/* ── reactive current URL from page store ── */
let currentUrl = $state<URL | null>(null);
if (browser) {
  page.subscribe((p) => {
    currentUrl = p.url;
  });
}

const updateSearchParams = createUrlUpdater(
  () => currentUrl || new URL('http://localhost'),
  (url) => goto(url, { replaceState: true, keepFocus: true, noScroll: true }),
  (url) => {
    currentUrl = url;
  },
);
function updateSearchParam(key: string, value: string | null) {
  updateSearchParams({ [key]: value });
}

export function getSearchUrlState(url?: URL) {
  const u = url || currentUrl || new URL('http://localhost');
  return {
    query: u.searchParams.get('q') || '',
    mode: parseSearchMode(u.searchParams.get('mode')),
    type: parseSearchType(u.searchParams.get('type')),
    species: u.searchParams.get('species') || '9606',
    entityWorkflow: parseEntityWorkflow(u.searchParams.get('entity_workflow')),
    filters: parseFiltersParam(u.searchParams.get('filters')),
  };
}

export function setSearchQuery(value: string) {
  updateSearchParam('q', value.trim() || null);
}

export function setSearchMode(value: SearchMode) {
  updateSearchParam('mode', value);
}

export function setSearchType(value: SearchType) {
  updateSearchParam('type', value);
}

export function setSpecies(value: string | null) {
  updateSearchParam('species', value);
}

export function setEntityWorkflow(value: EntityWorkflow) {
  updateSearchParam('entity_workflow', value === 'direct_lookup' ? null : value);
}

export function setSearchFilters(value: SearchFilters) {
  updateSearchParam('filters', serializeFiltersParam(value));
}

export function getInteractionsUrlState(url?: URL) {
  const u = url || currentUrl || new URL('http://localhost');
  const singleEntity = u.searchParams.get('entity');
  const multiEntities = u.searchParams.get('entities');
  const fromMany = parseEntityIdsParam(multiEntities);
  return {
    entityIds: u.searchParams.has('entities') ? fromMany : parseEntityIdsParam(singleEntity),
    filters: parseFiltersParam(u.searchParams.get('filters')),
  };
}

export function setInteractionsEntityIds(value: Array<string | number>) {
  updateSearchParams({
    entity: value.length === 1 ? String(value[0]) : null,
    entities: value.length === 1 ? null : value.map(String).join(','),
  });
}

export function setInteractionsFilters(value: SearchFilters) {
  updateSearchParam('filters', serializeFiltersParam(value));
}

export function getSelectionUrlState(url?: URL) {
  const u = url || currentUrl || new URL('http://localhost');
  return {
    tab: parseSelectionTab(u.searchParams.get('tab')),
    query: u.searchParams.get('q') || '',
    entityIds: parseEntityIdsParam(u.searchParams.get('entities')),
    filters: parseFiltersParam(u.searchParams.get('filters')),
  };
}

export function setSelectionTab(value: SelectionTab) {
  updateSearchParam('tab', value);
}

export function setSelectionQuery(value: string) {
  updateSearchParam('q', value.trim() || null);
}

export function setSelectionEntityIds(value: Array<string | number>) {
  updateSearchParam('entities', serializeEntityIdsParam(value) || '');
}

export function setSelectionFilters(value: SearchFilters) {
  updateSearchParam('filters', serializeFiltersParam(value));
}

/* ── localStorage-backed selection cache ── */

const SELECTION_STORAGE_KEY = 'omnipath-selection-entities';
const SELECTION_IDS_STORAGE_KEY = 'omnipath-selection-ids';
const SELECTION_ANNOTATIONS_STORAGE_KEY = 'omnipath-selection-annotations';
const SELECTION_ANNOTATION_IDS_STORAGE_KEY = 'omnipath-selection-annotation-ids';

type SelectionEntityCache = Record<string, SelectedEntity>;
type SelectionAnnotationCache = Record<string, SelectedAnnotation>;

function storage(): Storage | null {
  try {
    return browser ? localStorage : null;
  } catch {
    return null;
  }
}
function readSelectionCache() {
  return readStored(storage(), SELECTION_STORAGE_KEY, {}, isEntityCache) as SelectionEntityCache;
}
function writeSelectionCache(value: SelectionEntityCache) {
  writeStored(storage(), SELECTION_STORAGE_KEY, value);
}
function readSelectionIds() {
  return readStored(storage(), SELECTION_IDS_STORAGE_KEY, [], isStringList);
}
function writeSelectionIds(value: string[]) {
  writeStored(storage(), SELECTION_IDS_STORAGE_KEY, value);
}
function readSelectionAnnotationCache() {
  return readStored(
    storage(),
    SELECTION_ANNOTATIONS_STORAGE_KEY,
    {},
    isAnnotationCache,
  ) as SelectionAnnotationCache;
}
function writeSelectionAnnotationCache(value: SelectionAnnotationCache) {
  writeStored(storage(), SELECTION_ANNOTATIONS_STORAGE_KEY, value);
}
function readSelectionAnnotationIds() {
  return readStored(storage(), SELECTION_ANNOTATION_IDS_STORAGE_KEY, [], isStringList);
}
function writeSelectionAnnotationIds(value: string[]) {
  writeStored(storage(), SELECTION_ANNOTATION_IDS_STORAGE_KEY, value);
}

/* ── Global reactive selection state ── */

let entityCache = $state<SelectionEntityCache>(readSelectionCache());
let annotationCache = $state<SelectionAnnotationCache>(readSelectionAnnotationCache());
let fallbackEntityIds = $state<string[]>(readSelectionIds());
let fallbackAnnotationIds = $state<string[]>(readSelectionAnnotationIds());

export function getSelectionStore() {
  const ownedSelection = $derived(
    selectionFromUrl(currentUrl, fallbackEntityIds, fallbackAnnotationIds),
  );
  const entityIds = $derived(ownedSelection.entityIds);
  const annotationIds = $derived(ownedSelection.annotationIds);

  function syncSelectionUrl(entities: string[], annotations: string[]) {
    updateSearchParams({
      selection: '1',
      entities: entities.join(','),
      annotations: annotations.join(','),
    });
  }

  const selectedEntities = $derived<SelectedEntity[]>(
    entityIds.map((id) => entityCache[id] || { id, entityId: id, name: id }),
  );

  const selectedAnnotations = $derived<SelectedAnnotation[]>(
    annotationIds.map((id) => annotationCache[id] || { id, label: id }),
  );

  function setEntityIds(next: Array<string | number>) {
    const normalized = next.map(String);
    fallbackEntityIds = normalized;
    writeSelectionIds(normalized);
    syncSelectionUrl(normalized, annotationIds);
  }

  function setAnnotationIds(next: Array<string | number>) {
    const normalized = next.map(String);
    fallbackAnnotationIds = normalized;
    writeSelectionAnnotationIds(normalized);
    syncSelectionUrl(entityIds, normalized);
  }

  function addEntity(entity: SelectedEntity) {
    const id = String(entity.entityId ?? entity.id).trim();
    if (!id) return;

    const nextCache = {
      ...entityCache,
      [id]: {
        ...entityCache[id],
        ...entity,
        id,
        entityId: entity.entityId ?? id,
      },
    };
    entityCache = nextCache;
    writeSelectionCache(nextCache);

    if (entityIds.includes(id)) return;
    const next = [...entityIds, id];
    fallbackEntityIds = next;
    writeSelectionIds(next);
    syncSelectionUrl(next, annotationIds);
  }

  function toggleEntityGroup(keys: string[], groupKey?: string) {
    const remove = keys.every((key) => entityIds.includes(key));
    if (!remove) {
      const nextCache = { ...entityCache };
      for (const key of keys) {
        nextCache[key] = {
          ...nextCache[key],
          id: key,
          entityId: key,
          entityPk: key,
          groupKey,
          name: nextCache[key]?.name || key,
        };
      }
      entityCache = nextCache;
      writeSelectionCache(nextCache);
    }
    setEntityIds(
      remove
        ? entityIds.filter((key) => !keys.includes(key))
        : [...new Set([...entityIds, ...keys])],
    );
  }

  function addAnnotation(annotation: SelectedAnnotation) {
    const id = String(annotation.id).trim();
    if (!id) return;

    const nextCache = {
      ...annotationCache,
      [id]: {
        ...annotationCache[id],
        ...annotation,
        id,
        label: annotation.label || id,
      },
    };
    annotationCache = nextCache;
    writeSelectionAnnotationCache(nextCache);

    if (annotationIds.includes(id)) return;
    const next = [...annotationIds, id];
    fallbackAnnotationIds = next;
    writeSelectionAnnotationIds(next);
    syncSelectionUrl(entityIds, next);
  }

  function removeEntity(id: string) {
    const normalized = String(id).trim();
    if (!normalized) return;

    const { [normalized]: _, ...rest } = entityCache;
    entityCache = rest;
    writeSelectionCache(rest);

    const next = entityIds.filter((entry) => entry !== normalized);
    fallbackEntityIds = next;
    writeSelectionIds(next);
    syncSelectionUrl(next, annotationIds);
  }

  function removeAnnotation(id: string) {
    const normalized = String(id).trim();
    if (!normalized) return;

    const { [normalized]: _, ...rest } = annotationCache;
    annotationCache = rest;
    writeSelectionAnnotationCache(rest);

    const next = annotationIds.filter((entry) => entry !== normalized);
    fallbackAnnotationIds = next;
    writeSelectionAnnotationIds(next);
    syncSelectionUrl(entityIds, next);
  }

  function clearSelection() {
    entityCache = {};
    annotationCache = {};
    writeSelectionCache({});
    writeSelectionAnnotationCache({});
    fallbackEntityIds = [];
    fallbackAnnotationIds = [];
    writeSelectionIds([]);
    writeSelectionAnnotationIds([]);
    syncSelectionUrl([], []);
  }

  function isSelected(id: string) {
    return entityIds.includes(String(id).trim());
  }

  function isAnnotationSelected(id: string) {
    return annotationIds.includes(String(id).trim());
  }

  const selectedEntityPks = $derived(
    selectedEntities
      .map(entityPrimaryKeyFromSelection)
      .map((pk) => (pk == null ? '' : String(pk).trim()))
      .filter(Boolean),
  );

  return {
    get selectedEntities() {
      return selectedEntities;
    },
    get selectedAnnotations() {
      return selectedAnnotations;
    },
    get entityIds() {
      return entityIds;
    },
    get annotationIds() {
      return annotationIds;
    },
    get selectedEntityPks() {
      return selectedEntityPks;
    },
    setEntityIds,
    setAnnotationIds,
    addEntity,
    toggleEntityGroup,
    addAnnotation,
    removeEntity,
    removeAnnotation,
    clearSelection,
    isSelected,
    isAnnotationSelected,
    get selectionCount() {
      return entityIds.length;
    },
    get annotationCount() {
      return annotationIds.length;
    },
    get totalSelectionCount() {
      return entityIds.length + annotationIds.length;
    },
  };
}
