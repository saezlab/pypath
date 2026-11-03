import { listResources, summarizeResources } from '$lib/server/resource';
import type { PageServerLoad } from './$types';

export const load: PageServerLoad = async ({ parent }) => {
  const { selectedRelease } = await parent();
  const resources = await listResources(selectedRelease);
  return {
    resources,
    summary: summarizeResources(resources),
    resourcesUnavailable: resources.length === 0,
  };
};
