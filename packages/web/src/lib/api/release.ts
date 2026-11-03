import { page } from '$app/state';
import { sharedFetch } from './shared-fetch';

/** Include the page's release on every query, independent of other browser tabs. */
export function releaseFetch(input: string | URL, init?: RequestInit): Promise<Response> {
  const headers = new Headers(init?.headers);
  headers.set('X-OmniPath-Release', page.data.selectedRelease || 'latest');
  const path = String(input).split('?')[0];
  const shared =
    !init?.method ||
    init.method === 'GET' ||
    [
      '/entities/search',
      '/entities/scoped-facets',
      '/relations/search',
      '/relations/scoped-facets',
      '/selection/scope',
    ].some((suffix) => path.endsWith(suffix));
  return (shared ? sharedFetch : fetch)(input, { ...init, headers });
}

export function releaseUrl(path: string): string {
  const separator = path.includes('?') ? '&' : '?';
  return `${path}${separator}release=${encodeURIComponent(page.data.selectedRelease || 'latest')}`;
}
