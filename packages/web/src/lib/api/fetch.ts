import { sharedFetch } from './shared-fetch';

/** Fetch from the app's API proxy; reads and the listed searches share identical requests. */
export function apiFetch(input: string | URL, init?: RequestInit): Promise<Response> {
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
  return (shared ? sharedFetch : fetch)(input, init);
}
