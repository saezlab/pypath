import { error } from '@sveltejs/kit';
import { API_SERVICE_URL } from '$lib/server/api-config';
import type { LayoutServerLoad } from './$types';

export const load: LayoutServerLoad = async ({ url, cookies }) => {
  const saved = url.searchParams.get('release') || cookies.get('omnipath_release');
  const requested = saved === 'working' ? 'latest' : saved;
  let catalog: {
    default: string;
    releases: Array<{ version: string; resources: Record<string, string> }>;
  };
  try {
    const base = API_SERVICE_URL.endsWith('/') ? API_SERVICE_URL : `${API_SERVICE_URL}/`;
    const response = await fetch(new URL('releases', base));
    if (!response.ok) throw new Error('Release catalog unavailable');
    catalog = await response.json();
  } catch {
    return { releases: [], selectedRelease: requested || 'latest', releasesUnavailable: true };
  }
  const selectedRelease = requested || catalog.default;
  if (
    selectedRelease !== 'latest' &&
    !catalog.releases.some((release) => release.version === selectedRelease)
  ) {
    error(404, `Unknown OmniPath release: ${selectedRelease}`);
  }
  cookies.set('omnipath_release', selectedRelease, { path: '/', httpOnly: false, sameSite: 'lax' });
  return { releases: catalog.releases, selectedRelease, releasesUnavailable: false };
};
