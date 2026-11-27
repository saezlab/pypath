# Golden products

What the exclusion and MetSigDB products publish over the small releases of the tests. They
were frozen from the PostgreSQL builds that the Parquet products replaced, and the Parquet
products equal them row for row: the exclusion files name each fixture (`release`,
`hierarchy`, `cofactor_release`) with the seed files the run used (`vendored`: the packaged
ones, `small`: the three-line set of the tests); `metsigdb_membership.json` is the
membership table of `metsigdb_fixture.py` (MACdb set ids carry their prefix).

The version of the seed files is masked in the comparison: it is a commit or a hash.

`metabolite_class.json` holds the two tables of the metabolite class product over the
release of `metabolite_class_fixture.py`, built with all three signals of the rule file.
It was frozen from the first build of the product, not from an older one. The release
digest (`build_id`) is left out of the comparison.

A deliberate change of a product changes these files. Review the difference, then rewrite
them with `UPDATE_GOLDEN=1 pytest packages/subsets/tests/test_exclusion_parquet.py
packages/subsets/tests/test_metsigdb_parquet.py
packages/subsets/tests/test_metabolite_class.py`.
