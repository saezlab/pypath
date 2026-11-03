import { browser } from '$app/environment';

const UI_PREFERENCES_STORAGE_KEY = 'omnipath-ui-preferences';

type UiPreferencesValue = {
  showExplanations: boolean;
  /** Group entity results (genes with their products, chemicals by structure). */
  groupResults: boolean;
};

const DEFAULTS: UiPreferencesValue = { showExplanations: false, groupResults: true };

function readUiPreferences(): UiPreferencesValue {
  if (!browser) return DEFAULTS;
  try {
    const parsed = JSON.parse(
      localStorage.getItem(UI_PREFERENCES_STORAGE_KEY) || '{}',
    ) as Partial<UiPreferencesValue>;
    return {
      showExplanations: parsed.showExplanations === true,
      groupResults: parsed.groupResults !== false,
    };
  } catch {
    return DEFAULTS;
  }
}

let uiPreferences = $state<UiPreferencesValue>(readUiPreferences());

function writeUiPreferences(next: UiPreferencesValue) {
  if (!browser) return;
  try {
    localStorage.setItem(UI_PREFERENCES_STORAGE_KEY, JSON.stringify(next));
  } catch {
    // Private windows and blocked storage keep the setting for this visit only.
  }
}

export function getUiPreferences() {
  return {
    get showExplanations() {
      return uiPreferences.showExplanations;
    },
    setShowExplanations(showExplanations: boolean) {
      uiPreferences = { ...uiPreferences, showExplanations };
      writeUiPreferences(uiPreferences);
    },
    get groupResults() {
      return uiPreferences.groupResults;
    },
    setGroupResults(groupResults: boolean) {
      uiPreferences = { ...uiPreferences, groupResults };
      writeUiPreferences(uiPreferences);
    },
  };
}
