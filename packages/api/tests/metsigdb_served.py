"""A served release with a MetSigDB subset, for the route's tests.

The subset's Parquet is written through the layout code of the build
(``omnipath_core.fixtures.write_metsigdb_subsets``), then the release's database is built
from it, so the tests read what a served release holds: the tables, their manifest and the
digest that ties them.
"""

from omnipath_api.release_database import build_served_database, publish_newest
from omnipath_core.fixtures import (
    METSIGDB_SETS as SETS,
    metsigdb_key as key,
    metsigdb_membership_rows as membership_rows,
    write_dataset,
    write_metsigdb_subsets,
)

__all__ = ["SETS", "key", "membership_rows", "serve"]


def serve(root, *, subsets=True, membership=None):
    """A data root with release 1 published and its database built."""
    write_dataset(root)
    publish_newest(root, "1")
    if subsets:
        write_metsigdb_subsets(root, membership=membership)
    build_served_database(root, "1", keep=None, memory_limit="1GB", threads=1)
    return root
