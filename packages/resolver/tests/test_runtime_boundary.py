"""The resolver wheel owns its complete public facade without a builder."""

import os
import subprocess
import sys


def test_public_matching_does_not_import_builder(tmp_path):
    from identity_snapshot import build_snapshot

    library = build_snapshot(tmp_path / "identity")
    script = r"""
import importlib.abc
import sys

class BlockBuild(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "omnipath_build" or fullname.startswith("omnipath_build."):
            raise AssertionError("Runtime imported build: " + fullname)

sys.meta_path.insert(0, BlockBuild())
from omnipath_resolver import EntityResolver, RawEntityObservation
from omnipath_resolver.canonical import votes_for, get_policy
from omnipath_resolver.observations import observation_bundle

observation = RawEntityObservation(
    entity_key="protein", entity_type="protein", namespace="uniprot",
    identifier="P04637", taxon="9606",
)
resolver = EntityResolver()
try:
    results = resolver.resolve_entity_targets({"protein": observation})
    target, = results["protein"]
    assert target.canonical_namespace == "uniprot"
    assert target.canonical_identifier == "P04637"
    assert target.taxon == "9606"
    assert target.node_id is None
finally:
    resolver.close()
votes, observed = votes_for(observation, get_policy("protein"))
query, rows = observation_bundle("protein", observation, votes, observed, "gene_protein")
assert query["target"] == 2
assert rows and rows[0]["lookup_key"].startswith(bytes([1, 2, 1]))
# The identity library was written by the parent process; reading it exercises kv
# lookup, native decisions and record assembly under the blocker.
resolver = EntityResolver(sys.argv[1])
try:
    target, = resolver.resolve_entity_targets({"protein": observation})["protein"]
    assert target.matched
    assert target.protein_identifier == "P04637"
    assert target.label == "TP53"
    assert target.entity_type == "protein"
finally:
    resolver.close()
assert not any(name.startswith("omnipath_build") for name in sys.modules)
"""
    environment = dict(os.environ)
    environment.pop("OMNIPATH_LIBRARY_DIR", None)
    environment["OMNIPATH_IDENTITY_CACHE"] = str(tmp_path / "cache")
    subprocess.run([sys.executable, "-c", script, str(library)], check=True, env=environment)
