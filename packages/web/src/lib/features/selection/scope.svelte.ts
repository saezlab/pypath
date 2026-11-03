import { browser } from '$app/environment';
import { fetchSelectionScope, type SelectionScopeMode } from '$lib/api/client';

const SCOPE_SETTINGS_STORAGE_KEY = 'omnipath-selection-scope-settings';

export interface DerivedSelectionScope {
  selectedEntityIds: string[];
  selectedEntityPks: string[];
  selectedAnnotationIds: string[];
  annotationMatchedEntityIds: string[];
  scopedEntityIds: string[];
  scopedEntityPks: string[];
  termEntityPks: string[];
  ontologyTermIds: string[];
  criteriaCount: number;
  expandedEntityCount: number;
}

type SelectionScopeSettingsValue = {
  mode: SelectionScopeMode;
  expandSelection: boolean;
};

function readSelectionScopeSettings(): SelectionScopeSettingsValue {
  if (!browser) return { mode: 'union', expandSelection: true };
  try {
    const parsed = JSON.parse(
      localStorage.getItem(SCOPE_SETTINGS_STORAGE_KEY) || '{}',
    ) as Partial<SelectionScopeSettingsValue>;
    return {
      mode: parsed.mode === 'intersection' ? 'intersection' : 'union',
      expandSelection: parsed.expandSelection !== false,
    };
  } catch {
    return { mode: 'union', expandSelection: true };
  }
}

let scopeSettings = $state<SelectionScopeSettingsValue>(readSelectionScopeSettings());

function writeSelectionScopeSettings(next: SelectionScopeSettingsValue) {
  if (!browser) return;
  localStorage.setItem(SCOPE_SETTINGS_STORAGE_KEY, JSON.stringify(next));
}

export function getSelectionScopeSettings() {
  return {
    get mode() {
      return scopeSettings.mode;
    },
    get expandSelection() {
      return scopeSettings.expandSelection;
    },
    setMode(mode: SelectionScopeMode) {
      scopeSettings = { ...scopeSettings, mode };
      writeSelectionScopeSettings(scopeSettings);
    },
    setExpandSelection(expandSelection: boolean) {
      scopeSettings = { ...scopeSettings, expandSelection };
      writeSelectionScopeSettings(scopeSettings);
    },
  };
}

function normalizeIds(ids?: Array<string | number>) {
  return Array.from(new Set((ids || []).map((id) => String(id).trim()).filter(Boolean)));
}

export function deriveSelectionScope(input: {
  selectedEntityIds?: Array<string | number>;
  selectedEntityPks?: Array<string | number>;
  selectedAnnotationIds?: Array<string | number>;
  annotationMatchedEntityIds?: Array<string | number>;
  scopedEntityPks?: Array<string | number>;
  termEntityPks?: Array<string | number>;
  ontologyTermIds?: Array<string | number>;
  criteriaCount?: number;
  expandedEntityCount?: number;
}): DerivedSelectionScope {
  const selectedEntityIds = normalizeIds(input.selectedEntityIds);
  const selectedEntityPks = normalizeIds(input.selectedEntityPks);
  const selectedAnnotationIds = normalizeIds(input.selectedAnnotationIds);
  const annotationMatchedEntityIds = normalizeIds(input.annotationMatchedEntityIds);
  const scopedEntityPks = normalizeIds(input.scopedEntityPks);
  const termEntityPks = normalizeIds(input.termEntityPks);
  const ontologyTermIds = normalizeIds(input.ontologyTermIds);

  return {
    selectedEntityIds,
    selectedEntityPks,
    selectedAnnotationIds,
    annotationMatchedEntityIds,
    scopedEntityIds: selectedEntityIds,
    scopedEntityPks,
    termEntityPks,
    ontologyTermIds,
    criteriaCount: input.criteriaCount || 0,
    expandedEntityCount: input.expandedEntityCount || 0,
  };
}

export function getSelectionScope(
  selectedEntityIds: string[],
  selectedEntityPks: string[],
  selectedAnnotationIds: string[],
  options?: {
    resolveAnnotationEntities?: boolean;
    includeAssociatedEntities?: boolean;
    includeMembersParticipants?: boolean;
    mode?: SelectionScopeMode;
  },
) {
  const resolveAnnotationEntities = options?.resolveAnnotationEntities ?? true;
  const hasSelection = selectedEntityPks.length > 0 || selectedAnnotationIds.length > 0;
  const shouldResolve = browser && hasSelection;

  let annotationMatchedEntityIds = $state<string[]>([]);
  let scopedEntityPks = $state<string[]>(shouldResolve ? [] : selectedEntityPks);
  let termEntityPks = $state<string[]>([]);
  let ontologyTermIds = $state<string[]>([]);
  let criteriaCount = $state(0);
  let expandedEntityCount = $state(0);
  let isLoading = $state(shouldResolve);

  async function resolve() {
    if (!browser) {
      scopedEntityPks = selectedEntityPks;
      termEntityPks = [];
      ontologyTermIds = [];
      criteriaCount = 0;
      expandedEntityCount = 0;
      return;
    }
    if (selectedEntityPks.length === 0 && selectedAnnotationIds.length === 0) {
      annotationMatchedEntityIds = [];
      scopedEntityPks = [];
      termEntityPks = [];
      ontologyTermIds = [];
      criteriaCount = 0;
      expandedEntityCount = 0;
      return;
    }

    isLoading = true;
    try {
      const scope = await fetchSelectionScope({
        entityPks: selectedEntityPks,
        annotationTermIds: resolveAnnotationEntities ? selectedAnnotationIds : [],
        includeAssociatedEntities: options?.includeAssociatedEntities ?? true,
        includeMembersParticipants: options?.includeMembersParticipants ?? true,
        mode: options?.mode ?? 'union',
      });
      scopedEntityPks = scope.entityPks ?? [];
      termEntityPks = scope.termEntityPks ?? [];
      ontologyTermIds = scope.ontologyTermIds ?? [];
      criteriaCount = scope.criteriaCount ?? 0;
      expandedEntityCount = scope.expandedEntityCount ?? 0;
      annotationMatchedEntityIds = [];
    } catch (e) {
      console.error(e);
      annotationMatchedEntityIds = [];
      scopedEntityPks = selectedEntityPks;
      termEntityPks = [];
      ontologyTermIds = [];
      criteriaCount = 0;
      expandedEntityCount = 0;
    } finally {
      isLoading = false;
    }
  }

  $effect(() => {
    void resolve();
  });

  const scope = $derived(
    deriveSelectionScope({
      selectedEntityIds,
      selectedEntityPks,
      selectedAnnotationIds,
      annotationMatchedEntityIds: resolveAnnotationEntities ? annotationMatchedEntityIds : [],
      scopedEntityPks,
      termEntityPks,
      ontologyTermIds,
      criteriaCount,
      expandedEntityCount,
    }),
  );

  return {
    get selectedEntityIds() {
      return scope.selectedEntityIds;
    },
    get selectedEntityPks() {
      return scope.selectedEntityPks;
    },
    get selectedAnnotationIds() {
      return scope.selectedAnnotationIds;
    },
    get annotationMatchedEntityIds() {
      return scope.annotationMatchedEntityIds;
    },
    get scopedEntityIds() {
      return scope.scopedEntityIds;
    },
    get scopedEntityPks() {
      return scope.scopedEntityPks;
    },
    get termEntityPks() {
      return scope.termEntityPks;
    },
    get ontologyTermIds() {
      return scope.ontologyTermIds;
    },
    get criteriaCount() {
      return scope.criteriaCount;
    },
    get expandedEntityCount() {
      return scope.expandedEntityCount;
    },
    get isLoading() {
      return isLoading;
    },
  };
}
