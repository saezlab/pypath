/** Persistence is optional: restricted storage must not disable selection. */
export function readStored<T>(
  storage: Pick<Storage, 'getItem'> | null,
  key: string,
  fallback: T,
  valid: (value: unknown) => value is T,
): T {
  try {
    const value: unknown = JSON.parse(storage?.getItem(key) || 'null');
    return valid(value) ? value : fallback;
  } catch {
    return fallback;
  }
}
export function writeStored(storage: Pick<Storage, 'setItem'> | null, key: string, value: unknown) {
  try {
    storage?.setItem(key, JSON.stringify(value));
  } catch (error) {
    console.warn('Selection could not be persisted', error);
  }
}
export const isStringList = (value: unknown): value is string[] =>
  Array.isArray(value) && value.every((item) => typeof item === 'string');
export const isRecord = (value: unknown): value is Record<string, unknown> =>
  !!value && typeof value === 'object' && !Array.isArray(value);
