from omnipath_resolver.goslin import normalize


def test_goslin_identity_precision_and_chain_order():
    assert normalize("PC 18:1_16:0")["goslin"] == normalize("PC 16:0_18:1")["goslin"]
    assert normalize("PC(16:0/18:1)")["goslin"] == normalize("PC 16:0/18:1")["goslin"]
    assert normalize("PC 16:0/18:1")["goslin"] != normalize("PC 18:1/16:0")["goslin"]
    assert normalize("PC 34:1")["goslin"] != normalize("PC 16:0_18:1")["goslin"]
    assert normalize("PC 16:0_18:1")["goslin"] != normalize("PC 16:0/18:1")["goslin"]
    assert normalize("5-methyl-octadecanoic acid")["goslin"] == "full_structure:FA 18:0;5Me"
    assert normalize("water")["goslin"] is None
    assert normalize("PC")["goslin"] is None
    assert normalize("PC 16:0/18:1[M+H]1+")["status"] == "adduct_or_isotope"


def test_native_abbreviation_lists(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from omnipath_build.identity.goslin import build

    hub = tmp_path / "swiss.parquet"
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "hub_id": "SLM:1",
                    "source_type": "lipid_shorthand",
                    "source_id": "Cer(d18:0/24:0(2OH)) | Cer(d18:0/24:0-2OH)",
                },
            ]
        ),
        hub,
    )
    out = tmp_path / "result"
    build({"swisslipids": hub}, out, workers=1, memory="256MB")
    claims = pq.read_table(out / "claims.parquet").to_pylist()
    normalized = pq.read_table(out / "normalized.parquet").to_pylist()
    assert len(normalized) == 2
    assert all(r["status"] == "parsed" for r in normalized)
    # Goslin reports different precision for these native alternatives.
    # Retain the more informative descriptor rather than forcing equivalence.
    assert len(claims) == 1
    assert claims[0]["goslin"] == normalize("Cer(d18:0/24:0(2OH))")["goslin"]
    assert all(r["record_id"] == "swisslipids:SLM:1" for r in claims)


def test_cache_reused_across_builds_and_only_new_names_parsed(tmp_path, monkeypatch):
    import json
    import pyarrow as pa
    import pyarrow.parquet as pq
    from omnipath_build.identity import goslin

    hub = tmp_path / "hub.parquet"
    cache = tmp_path / "shared-cache"

    def write_names(names):
        pq.write_table(
            pa.Table.from_pylist(
                [
                    {"hub_id": str(i), "source_type": "name", "source_id": name}
                    for i, name in enumerate(names)
                ]
            ),
            hub,
        )

    def summary(out):
        return json.loads((out / "summary.json").read_text())

    write_names(["PC 16:0/18:1", "water"])
    first = tmp_path / "first"
    goslin.build({"hmdb": hub}, first, workers=1, memory="256MB", cache_dir=cache)
    assert summary(first)["parsed_names"] == 2
    assert summary(first)["cached_names"] == 0

    second = tmp_path / "second"
    with monkeypatch.context() as patch:

        def forbidden(*_):
            raise AssertionError("A cached name was sent to the parser")

        patch.setattr(goslin, "parse_batch", forbidden)
        goslin.build({"hmdb": hub}, second, workers=1, memory="256MB", cache_dir=cache)
    assert summary(second)["parsed_names"] == 0
    assert summary(second)["cached_names"] == 2  # failed parses are cached too
    assert (
        pq.read_table(first / "claims.parquet").to_pylist()
        == pq.read_table(second / "claims.parquet").to_pylist()
    )

    write_names(["PC 16:0/18:1", "water", "PE 36:2"])
    third = tmp_path / "third"
    goslin.build({"hmdb": hub}, third, workers=1, memory="256MB", cache_dir=cache)
    assert summary(third)["parsed_names"] == 1
    assert summary(third)["cached_names"] == 2


def test_cache_persists_batches_and_invalidates_changed_normalizer(tmp_path):
    from omnipath_resolver.goslin_cache import NormalizationCache

    result = normalize("PC 34:1")
    with NormalizationCache(tmp_path, normalize) as cache:
        first_path = cache.path
        cache.store([result])
    with NormalizationCache(tmp_path, normalize) as cache:
        assert cache.lookup(["PC 34:1"]) == [result]

    def changed_normalizer(name):
        return normalize(name)

    with NormalizationCache(tmp_path, changed_normalizer) as cache:
        assert cache.path != first_path
        assert cache.lookup(["PC 34:1"]) == []


def test_general_database_parse_must_match_the_record_structure(tmp_path):
    """Water's synonym HOH parses as a hydroxy fatty acid; its formula says otherwise."""
    import json
    import pyarrow as pa
    import pyarrow.parquet as pq
    from omnipath_build.identity.goslin import build

    chebi = tmp_path / "chebi.parquet"
    rows = [
        ("CHEBI:15377", "name", "water"),
        ("CHEBI:15377", "synonym", "HOH"),
        ("CHEBI:15377", "formula", "H2O"),
        ("CHEBI:16196", "name", "oleic acid"),
        ("CHEBI:16196", "synonym", "FA 18:1(9Z)"),
        ("CHEBI:16196", "formula", "C18H33O2"),  # the anion: hydrogen does not count
        ("CHEBI:9", "synonym", "PC 34:1"),  # no structure stated: nothing to contradict
        # Lipid notation is the source's own statement, even against its formula.
        ("CHEBI:89504", "name", "PC(36:2)"),
        ("CHEBI:89504", "formula", "C44H84NO7P"),
        ("CHEBI:17120", "name", "hexanoate"),  # read as WE 6:0/0:0: never a claim
    ]
    pq.write_table(
        pa.Table.from_pylist([dict(hub_id=h, source_type=t, source_id=v) for h, t, v in rows]),
        chebi,
    )
    out = tmp_path / "result"
    build({"chebi": chebi}, out, workers=1, memory="256MB")
    claims = {r["record_id"]: r["goslin"] for r in pq.read_table(out / "claims.parquet").to_pylist()}
    assert "chebi:CHEBI:15377" not in claims
    assert "chebi:CHEBI:17120" not in claims
    assert claims["chebi:CHEBI:16196"] == normalize("FA 18:1(9Z)")["goslin"]
    assert claims["chebi:CHEBI:9"] == normalize("PC 34:1")["goslin"]
    assert claims["chebi:CHEBI:89504"] == normalize("PC(36:2)")["goslin"]
    summary = json.loads((out / "summary.json").read_text())
    assert summary["structure_rejected"] == 1
