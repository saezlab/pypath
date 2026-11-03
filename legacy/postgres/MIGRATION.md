# PostgreSQL migration source

This directory retains the pre-migration PostgreSQL pipeline, its tests and
configuration. It is outside the active workspace so its `omnipath_build`
package cannot shadow the resolved-Parquet resource pipeline.

Use this implementation as the source for the later `omnipath_postgres` and
`omnipath_subsets` packages. The first migration milestone does not load a
PostgreSQL database or validate these legacy commands.

The retained Makefile, lockfile and setup instructions are historical. Their
original paths and environment assumptions are not the active workspace setup.
