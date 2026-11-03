"""Frozen policy fixtures and invariants for the sole Rust resolver."""

from dataclasses import asdict
from contextlib import closing
import random

import pytest

from library_fixture import build_fixture_library, WATER, ASPIRIN
from test_canonical import obs
from omnipath_resolver.canonical.match import LibraryMatcher


def outputs(matcher, batch):
    return {k: asdict(v) for k, v in matcher.match(batch).items()}


@pytest.mark.parametrize("defer", [True, False])
def test_resolution_golden_and_batch_boundaries(tmp_path, defer):
    library = build_fixture_library(tmp_path)
    rng = random.Random(803)
    chemicals = [
        ("bigg_metabolite", "h2o"),
        ("metanetx", "MNXM2"),
        ("inchikey", WATER),
        ("inchikey", ASPIRIN),
        ("chebi", "15377"),
        ("pubchem", "missing"),
        ("name", "water"),
    ]
    proteins = [
        ("uniprot", "P04637"),
        ("entrez", "7157"),
        ("genesymbol", "TP53"),
        ("genesymbol-syn", "TP53"),
        ("uniprot", "absent"),
    ]
    batch = {}
    for i in range(50):
        choices = chemicals if i % 2 else proteins
        ns, ident = rng.choice(choices)
        aliases = [rng.choice(choices) for _ in range(rng.randrange(5))]
        entity = obs(
            "chemical_entity" if i % 2 else "protein",
            ns,
            ident,
            aliases,
            taxon=rng.choice(["9606", "10090", ""]),
        )
        batch[str(i)] = entity
    with closing(LibraryMatcher(library, defer_aliases=defer)) as fast:
        expected = outputs(fast, batch)
        assert outputs(fast, batch) == expected  # answered from the matcher's memo
        assert fast.metrics["memo_hits"] == len(batch)
        # The memo is cleared so each batch below is decided again by the runtime.
        fast._memo.clear()
        assert outputs(fast, dict(reversed(list(batch.items())))) == expected
        # Repeated keys, cold misses, and cache persistence across smaller batches.
        for start in range(0, 50, 17):
            fast._memo.clear()
            subset = dict(list(batch.items())[start : start + 17])
            assert outputs(fast, subset) == {k: expected[k] for k in subset}
        assert fast.metrics["batches"] >= 3
