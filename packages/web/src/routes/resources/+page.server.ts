import { listResources, summarizeResources } from '$lib/server/resource';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async () => {
  // The list needs neither file records nor resolution statistics (they were half the
  // page): the details dialog loads a resource's files when it opens.
  const resources = (await listResources()).map(({ resolution_stats: _stats, ...resource }) => ({
    ...resource,
    files: [],
  }));
  return {
    resources,
    summary: summarizeResources(resources),
    resourcesUnavailable: resources.length === 0,
  };
};
