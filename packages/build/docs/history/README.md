# Historical material

`notebooks/exploration.ipynb` retains the original exploratory notebook code with
saved outputs removed. It is not an entry point for the current pipeline.

`exploratory-artifacts-20261003.json` records the sizes and hashes of the previous
single-space root Parquet, original notebook and orphan root npm lockfile. Their
original bytes remain in the ignored local `local-artifacts/refactoring-20261003/`
directory, and are recoverable from the recorded baseline Git commit in any full
clone. This cleanup removes roughly 87 MB from the current source tree, not from
historical Git objects. No history rewrite was performed.

The dated migration reports and build log retain their historical paths and runtime
commits. The old numbered pipeline guide is labelled historical; current workflows
are in the root README and docs/architecture.md.
