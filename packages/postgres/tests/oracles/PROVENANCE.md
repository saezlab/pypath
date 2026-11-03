# PostgreSQL parity oracles

All tests read tracked sources through the SHA256-checking loader in this directory.
The main 9f9bb709c764 source remains untouched under `legacy/postgres/omnipath_build`.
`cosmos_project_edges_951d500.sql` is the exact deployed source from commit
951d500b2535547076f292db949f436872b9fda2, original path
`packages/omnipath_postgres/src/omnipath_postgres/main_compat/cosmos/sql/project_edges.sql`.
SHA256: `987c6870938d52e1dd5ac6af1837f30ca273e3a1b0fd4f7d6032bf2d8b734071`.
It was copied once from that commit, preserving bytes; tests need no Git history.
These files must not be regenerated from candidate SQL.

`source_sha256.json` records all tracked Python/SQL/YAML source files in the frozen
main and CV-term oracles. Every imported reference source is checked before execution.
