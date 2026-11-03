const STORAGE_KEY = 'omnipath-admin-password';

// The admin secret lives only for the browser session; never persist it.
function forgetLegacyCopy(): void {
  try {
    if (typeof localStorage !== 'undefined') localStorage.removeItem(STORAGE_KEY);
  } catch {
    // Storage can be unavailable (private mode, blocked site data).
  }
}

export function getAdminPassword(): string {
  forgetLegacyCopy();
  try {
    return typeof sessionStorage === 'undefined' ? '' : sessionStorage.getItem(STORAGE_KEY) || '';
  } catch {
    return '';
  }
}

export function setAdminPassword(value: string): void {
  forgetLegacyCopy();
  try {
    if (typeof sessionStorage !== 'undefined') sessionStorage.setItem(STORAGE_KEY, value);
  } catch {
    // Without storage the secret must be re-entered; never fall back to localStorage.
  }
}

export function clearAdminPassword(): void {
  forgetLegacyCopy();
  try {
    if (typeof sessionStorage !== 'undefined') sessionStorage.removeItem(STORAGE_KEY);
  } catch {
    // Nothing stored.
  }
}

export function adminHeaders(init?: HeadersInit): Headers {
  const headers = new Headers(init);
  const password = getAdminPassword();
  if (password) headers.set('x-admin-secret', password);
  return headers;
}

export function adminFetch(input: RequestInfo | URL, init?: RequestInit): Promise<Response> {
  return fetch(input, { ...init, headers: adminHeaders(init?.headers) });
}
