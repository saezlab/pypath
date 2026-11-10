import type { LayoutServerLoad } from './$types';

// One release is served (the newest); a choice remembered from when the explorer offered
// several would pin a release the API no longer serves.
export const load: LayoutServerLoad = async ({ cookies }) => {
  if (cookies.get('omnipath_release')) cookies.delete('omnipath_release', { path: '/' });
  return {};
};
