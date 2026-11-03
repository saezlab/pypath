# analysis/

Reusable, read-only SQL evaluations of a built OmniPath DB. Run against any built
instance; they disable parallel gather so they work on a web PG too.

- **`gene_protein_resolution_report.sql`** — gene/protein entity resolution:
  unique genes; protein mentions grounded in a real gene vs virtual "unknown"
  gene vs unresolved; **why** the unresolved fail (by the id they carry × taxon);
  all **by resource**; and each resource's primary source id type + translation
  path. A written-up snapshot with interpretation lives in the spec cycle at
  `specs/003-curation-and-coverage/resolution-state-report.md` — re-run the script
  each phase to track movement (esp. the UniProt-coverage tail).

```bash
PGPASSWORD=… psql -h localhost -p <build-or-web-PG-port> -U omnipath -d omnipath \
  -f analysis/gene_protein_resolution_report.sql
```
