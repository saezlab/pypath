# Golden products

What the exclusion and MetSigDB products publish over the small releases of the tests. They
were frozen from the PostgreSQL builds that the Parquet products replaced, and the Parquet
products equal them row for row: the exclusion files name each fixture (`release`,
`hierarchy`, `cofactor_release`) with the seed files the run used (`vendored`: the packaged
ones, `small`: the three-line set of the tests); `metsigdb_membership.json` is the
membership table of `metsigdb_fixture.py` (MACdb set ids carry their prefix).

The version of the seed files is masked in the comparison: it is a commit or a hash.

A deliberate change of a product changes these files. Review the difference, then rewrite
them with `UPDATE_GOLDEN=1 pytest packages/subsets/tests/test_exclusion_parquet.py
packages/subsets/tests/test_metsigdb_parquet.py`.
