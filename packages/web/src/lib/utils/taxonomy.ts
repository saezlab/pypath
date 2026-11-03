/** Labels come from the release's taxonomy reference; IDs remain filter values. */
export function formatTaxonomy(id: string | null | undefined, name?: string | null): string | null {
  if (!id || id === '0') return null;
  return name ? `${name.charAt(0).toUpperCase() + name.slice(1)} (${id})` : `NCBI taxon ${id}`;
}
