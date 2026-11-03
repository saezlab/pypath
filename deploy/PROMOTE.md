# Build → web promotion (volume hand-off)

Build PGs and web (present) PGs are the **same image** (`omnipath-build-postgres:18`,
with `pg_roaringbitmap` + `rdkit`). They differ only in runtime settings (compose
`command:` flags — memory/WAL/parallelism), port, network, and the **data volume**
mounted. Postgres settings are startup flags, not baked into the data dir — so a
build is promoted by handing its **data volume** to a web stack, which mounts it
with web-tuned settings. No `pg_dump`/`pg_restore`.

## Standard flow
1. **Build into a dedicated result volume.** Run the build against a Postgres whose
   data volume is dedicated to that build (e.g. a dated volume), not a slot you
   want to keep. (The per-instance build slots `omnipath-build[-devN]_postgres18_data`
   are fine if you intend to hand that slot off.)
2. **Promote (hand-off).** Stop the build PG (one Postgres per data dir), point the
   target web stack's volume at it, restart:
   ```
   deploy/promote-build.sh <source-data-volume> <target-instance>
   ```
   Re-points `~/instances/<target>/.env:POSTGRES_DATA_VOLUME_NAME` and
   `systemctl --user restart omnipath-present@<target>`. Instant.
3. **Keep the source build slot?** Use `--copy` (raw volume copy to a dated
   `promoted-*` volume, source preserved + resumed). Minutes, not the ~hour a
   `pg_dump`/`pg_restore` of a multi-GB build takes.

## Notes
- Versions must match (build + web both Postgres 18).
- The old web volume is detached (not deleted) — `docker volume rm` it when done.
- systemd boot: build PGs `omnipath-build@<i>`; web stacks `omnipath-present@<i>`
  (enable + `loginctl enable-linger`). All PGs also carry `restart: unless-stopped`.
- For a one-off copy of an existing live build slot (can't stop it), `pg_dump -Fc`
  + `pg_restore -j` still works (see git history) — just slower.
