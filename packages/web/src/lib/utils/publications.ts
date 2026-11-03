export type Publication = { id: string; url?: string; sources: string[] };

export function isPublicationTerm(term: string): boolean {
  return ['publications', 'pubmed', 'pmid', 'doi'].includes(
    term.toLowerCase().replace(/^biolink:/, ''),
  );
}

export function groupPublications(
  rows: { term: string; value: string; source?: string }[],
): Publication[] {
  const publications = new Map<string, Publication>();
  for (const row of rows) {
    if (!isPublicationTerm(row.term)) continue;
    const term = row.term.toLowerCase().replace(/^biolink:/, '');
    const values = ['pubmed', 'pmid'].includes(term)
      ? row.value.split(/[;,|\s]+/).filter(Boolean)
      : [row.value.trim()];
    for (const value of values) {
      if (!value) continue;
      const pmid =
        /^(?:PMID:|https?:\/\/pubmed\.ncbi\.nlm\.nih\.gov\/)(\d+)\/?$/i.exec(value)?.[1] ||
        (['pubmed', 'pmid'].includes(term) && /^\d+$/.test(value) ? value : undefined);
      const doi =
        /^(?:doi:|https?:\/\/(?:dx\.)?doi\.org\/)(10\.\d{4,9}\/.+)$/i.exec(value)?.[1] ||
        (term === 'doi' && /^10\.\d{4,9}\/.+/.test(value) ? value : undefined);
      const id = pmid ? `PMID:${pmid}` : doi ? `DOI:${doi}` : value;
      const key = id.toLowerCase();
      const entry = publications.get(key) || {
        id,
        url: pmid
          ? `https://pubmed.ncbi.nlm.nih.gov/${pmid}/`
          : doi
            ? `https://doi.org/${encodeURI(doi).replace(/#/g, '%23').replace(/\?/g, '%3F')}`
            : undefined,
        sources: [],
      };
      if (row.source && !entry.sources.includes(row.source)) entry.sources.push(row.source);
      publications.set(key, entry);
    }
  }
  return [...publications.values()]
    .map((p) => ({ ...p, sources: p.sources.sort() }))
    .sort((a, b) => a.id.localeCompare(b.id, undefined, { numeric: true }));
}
