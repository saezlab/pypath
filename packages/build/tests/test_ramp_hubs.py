from omnipath_build.hubs.sources.ramp import assertions as ramp


def test_ramp_separates_domains_and_retains_conflicting_assertions():
    a = list(
        ramp("source", {"rampId": "RAMP_C_1", "IDtype": "hmdb", "sourceId": "hmdb:HMDB0000001"})
    )
    assert a[0][0:3] == ("ramp", "RAMP_C_1", "0")
    assert a[0][3]["hmdb"] == "HMDB0000001"
    assert (
        list(ramp("source", {"rampId": "RAMP_G_1", "IDtype": "entrez", "sourceId": "entrez:123"}))[
            0
        ][0]
        == "ramp_gene"
    )
    assert not list(ramp("pathway", {"pathwayRampId": "RAMP_P_1"}))


def test_registered_ramp_export_uses_native_cache(tmp_path, monkeypatch):
    import sqlite3
    import pyarrow.parquet as pq
    from omnipath_build.hubs.export import export_hubs
    from omnipath_build import discovery

    cache = tmp_path / "cache"
    (cache / "ramp").mkdir(parents=True)
    with sqlite3.connect(cache / "ramp/RaMP_SQLite_v3.0.7.sqlite") as con:
        con.execute("CREATE TABLE source(rampId TEXT, IDtype TEXT, sourceId TEXT)")
        con.executemany(
            "INSERT INTO source VALUES (?,?,?)",
            [
                ("RAMP_C_1", "chebi", "chebi:15377"),
                ("RAMP_G_1", "entrez", "entrez:7157"),
            ],
        )
        for name in ("analyte", "metabolite_class", "chem_props", "analytehaspathway"):
            con.execute(f"CREATE TABLE {name}(rampId TEXT)")
    monkeypatch.setattr(discovery, "setup_pypath_cache", lambda: cache)
    out = tmp_path / "hubs"
    export_hubs(out, hubs=["ramp", "ramp_gene"], include_dictionaries=False)
    chemical = pq.read_table(out / "ramp.parquet").to_pylist()
    genes = pq.read_table(out / "ramp_gene.parquet").to_pylist()
    assert {r["hub_id"] for r in chemical} == {"RAMP_C_1"}
    assert {r["hub_id"] for r in genes} == {"RAMP_G_1"}
    assert any(r["source_type"] == "entrez" and r["source_id"] == "7157" for r in genes)


def test_download_and_internal_cache_share_one_directory(tmp_path, monkeypatch):
    from omnipath_build.discovery import setup_pypath_cache
    from pypath.share.cache import get_cachedir
    from pypath.share import settings

    previous = settings.get("cachedir")
    monkeypatch.setenv("PYPATH_DOWNLOAD_DATADIR", str(tmp_path / "cache"))
    try:
        assert setup_pypath_cache() == get_cachedir() == tmp_path / "cache"
    finally:
        settings.setup(cachedir=previous)
