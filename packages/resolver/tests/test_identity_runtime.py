"""IdentityRuntime over a synthetic omnipath-identity-v2 layout (spec sections 3a, 3b, 4)."""

from __future__ import annotations

import hashlib
import json
import shutil
import sqlite3

import pytest

from identity_snapshot import (
    ASPIRIN,
    ETHANOL,
    FINGERPRINT,
    LIPID,
    LIPID3,
    AMINO,
    ANION,
    FATTY,
    NADH,
    NADH_2,
    NAMELESS,
    ONLY_PC,
    OTHERS,
    SYSTEMATIC,
    TWO_A,
    TWO_C,
    TWO_D,
    WATER,
    build_snapshot,
    part,
)
from omnipath_resolver import RawEntityObservation
from omnipath_resolver import identity_kv, identity_runtime
from omnipath_resolver.canonical import votes_for
from omnipath_resolver.canonical.match import LibraryMatcher
from omnipath_resolver.canonical.policy import CHEMICAL_POLICY
from omnipath_resolver.identity_runtime import (
    IdentityRuntime,
    decode_key,
    entity_num,
    open_runtime,
)
from omnipath_resolver.identity_kv import KvMissing
from omnipath_resolver.observations import CODES, key, observation_bundle

CODE_NAMES = {code: ns for ns, code in CODES.items()}
HUMAN, MOUSE = "9606", "10090"


LAYOUTS = {  # id: (shards, hubs with a kv store in the older format v2)
    "one-shard": (1, ()),
    "sixteen-shards": (16, ()),
    "format-v2": (1, "all"),
    "mixed-formats": (1, ("entrez",)),  # entrez v2, uniprot v3: a library part way to v3
}


@pytest.fixture(scope="module", params=list(LAYOUTS))
def snapshot(tmp_path_factory, request):
    """The same layout with one kv environment per store, with 16 shards per large store, and
    with hub kv stores of the older format (symbol rows not keyed by taxon)."""
    shards, v2_hubs = LAYOUTS[request.param]
    return build_snapshot(tmp_path_factory.mktemp("identity"), shards=shards, v2_hubs=v2_hubs)


@pytest.fixture(scope="module")
def snapshot_v2(tmp_path_factory):
    return build_snapshot(tmp_path_factory.mktemp("identity-v2"), v2_hubs="all")


@pytest.fixture
def runtime_v2(snapshot_v2, tmp_path):
    rt = IdentityRuntime(snapshot_v2, cache_dir=tmp_path / "cache")
    yield rt
    rt.close()


@pytest.fixture
def runtime(snapshot, tmp_path):
    rt = IdentityRuntime(snapshot, cache_dir=tmp_path / "cache")
    yield rt
    rt.close()


def ids(value):
    return [c[1] for c in value["candidates"]]


def look(rt, target, route, ns, identifier, scope=""):
    k = key(target, route, ns, scope, identifier)
    return rt.lookup_many([k])[k]


# ---------------------------------------------------------------------------- keys


def test_key_decoding_unscoped_and_scoped():
    assert decode_key(key(1, 1, "chebi", "", "CHEBI:1"), CODE_NAMES) == (
        1, 1, "chebi", None, "CHEBI:1"
    )  # fmt: skip
    assert decode_key(key(2, 2, "entrez", "9606", "7157"), CODE_NAMES) == (
        2, 2, "entrez", "9606", "7157"
    )  # fmt: skip
    assert decode_key(key(2, 1, "genesymbol", "10090", "Trp53"), CODE_NAMES)[3:] == (
        "10090",
        "Trp53",
    )
    for bad in (b"", b"\x02\x01\x01\x00\x01\x00x", key(1, 1, "chebi", "", "x")[:6]):
        with pytest.raises(ValueError):
            decode_key(bad, CODE_NAMES)
    with pytest.raises(ValueError):
        decode_key(bytes([1, 1, 1, 0xFF, 0xFF, 0]) + b"x", CODE_NAMES)


def test_only_identity_directories_open(snapshot, tmp_path):
    rt = open_runtime(snapshot, cache_dir=tmp_path)
    try:
        assert isinstance(rt, IdentityRuntime)
    finally:
        rt.close()
    with pytest.raises(OSError):
        IdentityRuntime(tmp_path)  # no manifest
    # the first-generation single-snapshot format is no longer served
    (tmp_path / "manifest.json").write_text(json.dumps({"format": "omnipath-identity-v1"}))
    with pytest.raises(ValueError):
        IdentityRuntime(tmp_path)


def test_hub_index_paths_relative_and_absolute(snapshot):
    rt = IdentityRuntime(snapshot, cache_dir=snapshot.parent / "c")
    try:
        assert rt.hub_dirs["uniprot"].is_absolute()  # written absolute
        assert (rt.hub_dirs["chebi"] / "kv" / "manifest.json").is_file()  # written relative
    finally:
        rt.close()


def test_runtime_reads_only_the_kv_stores(snapshot, tmp_path):
    # Delete every Parquet file of the layout: the runtime must not need any of them.
    root = tmp_path / "tree"
    shutil.copytree(snapshot.parent.parent, root)
    copy = root / "identity" / snapshot.name
    for path in root.rglob("*.parquet"):
        path.unlink()
    assert not list(root.rglob("*.parquet"))
    rt = IdentityRuntime(copy, cache_dir=tmp_path / "c")
    try:
        assert not hasattr(rt, "_db")
        k = key(1, 1, "chebi", "", "CHEBI:15377")
        assert ids(rt.lookup_many([k])[k]) == [f"inchikey:{WATER}"]
        assert rt.record(f"inchikey:{WATER}")["label"] == "water"
        assert rt.record("entrez:7157")["label"] == "TP53"
    finally:
        rt.close()


def test_missing_kv_fails_clearly(snapshot, tmp_path):
    root = tmp_path / "tree"
    shutil.copytree(snapshot.parent.parent, root)
    copy = root / "identity" / snapshot.name
    shutil.rmtree(root / "hubindex" / "hmdb" / "0123456789ab" / "kv")
    with pytest.raises(KvMissing, match="hmdb.*build-hub-kv"):
        IdentityRuntime(copy, cache_dir=tmp_path / "c")
    shutil.rmtree(copy / "kv")
    with pytest.raises(KvMissing, match="build-identity-kv"):
        IdentityRuntime(copy, cache_dir=tmp_path / "c")


def test_stale_kv_is_refused(snapshot, tmp_path):
    root = tmp_path / "tree"
    shutil.copytree(snapshot.parent.parent, root)
    copy = root / "identity" / snapshot.name
    manifest = root / "hubindex" / "chebi" / "0123456789ab" / "manifest.json"
    manifest.write_text(manifest.read_text() + " ")  # the index changed after its kv was built
    with pytest.raises(KvMissing, match="different index"):
        IdentityRuntime(copy, cache_dir=tmp_path / "c")


def test_only_the_hashed_shards_are_opened(snapshot, tmp_path):
    rt = IdentityRuntime(snapshot, cache_dir=tmp_path / "c")
    try:
        k1, k2 = key(1, 1, "chebi", "", "CHEBI:15377"), key(1, 1, "hmdb", "", "HMDB0002111")
        out = rt.lookup_many([k1, k2])
        assert ids(out[k1]) == [f"inchikey:{WATER}"] and ids(out[k2]) == [f"inchikey:{WATER}"]
        shards = rt.hub_kv["chebi"].id_shards
        want = 0 if shards == 1 else int(part("CHEBI:15377")[0], 16)
        assert rt.hub_kv["chebi"].opened() >= {("id", want)}
        assert {s for kind, s in rt.hub_kv["chebi"].opened() if kind == "id"} == {want}
        assert rt.hub_kv["pubchem"].opened() == set()  # hubs outside the key's domain stay closed
    finally:
        rt.close()


def test_miss_and_batch_shape(runtime):
    k = [key(1, 1, "chebi", "", "CHEBI:404"), key(2, 1, "hgnc", "", "HGNC:0")]
    out = runtime.lookup_many(k)
    assert set(out) == set(k)
    assert out[k[0]] == dict(candidates=[], gene=False, products=False)
    assert out[k[1]] == dict(candidates=[], gene=True, products=False)
    assert runtime.lookup_many([]) == {}


# ------------------------------------------------------------------------ admission


def test_native_claim_and_fallback_rows(runtime):
    water = f"inchikey:{WATER}"
    assert ids(look(runtime, 1, 1, "pubchem", "962")) == [water]  # native, own hub only
    assert ids(look(runtime, 1, 1, "inchikey", WATER)) == [water]
    assert ids(look(runtime, 1, 1, "cas", "7732-18-5")) == [water]  # claim, any chemical hub
    assert ids(look(runtime, 1, 1, "kegg", "C00001")) == [water]  # lone fallback survives
    # an xref of another hub is not a lookup key: HMDB ids are native only in the hmdb hub
    assert look(runtime, 1, 1, "inchi", "InChI=1S/x")["candidates"] == []


def test_native_beats_fallback(runtime):
    # bigg:x1 owns the id; bigg:y1 merely claims it (fallback)
    assert ids(look(runtime, 1, 1, "bigg", "x1")) == [f"inchikey:{ETHANOL}"]


def test_primary_accession_beats_secondary(runtime):
    assert ids(look(runtime, 2, 1, "uniprot", "Q99999")) == ["uniprot:Q99999"]
    assert ids(look(runtime, 2, 1, "uniprot", "Q15086")) == ["uniprot:P04637"]
    # a secondary accession of two entries stays plural (the kernel abstains)
    assert sorted(ids(look(runtime, 2, 1, "uniprot", "Q88888"))) == [
        "uniprot:P04637",
        "uniprot:P0DP23",
    ]
    assert ids(look(runtime, 2, 1, "uniprot-sec", "Q15086")) == ["uniprot:P04637"]
    assert ids(look(runtime, 2, 1, "uniprot_entry", "P53_HUMAN")) == ["uniprot:P04637"]


def test_exact_symbol_beats_synonym_within_taxon(runtime):
    for ns in ("genesymbol", "genesymbol-syn"):
        assert ids(look(runtime, 2, 1, ns, "ABC1", HUMAN)) == ["entrez:901"]
        assert ids(look(runtime, 2, 1, ns, "XYZ", HUMAN)) == ["entrez:902"]  # synonym only
        # TP53 is exact for human and only a synonym for mouse: each taxon keeps its own
        assert ids(look(runtime, 2, 1, ns, "TP53", HUMAN)) == ["entrez:7157"]
        assert ids(look(runtime, 2, 1, ns, "TP53", MOUSE)) == ["entrez:22059"]


def test_version_stripped_rows(runtime):
    for identifier in ("NP_000537", "NP_000537.3"):
        assert ids(look(runtime, 2, 1, "refseq_protein", identifier)) == ["uniprot:P04637"]
    # the gene hub's own refseq_protein rows are product ids and never keys
    assert ids(look(runtime, 2, 1, "refseq", "NM_000546")) == ["entrez:7157"]


# ---------------------------------------------------------------------------- scope


def test_scope_filters_on_entity_taxon(runtime):
    assert ids(look(runtime, 2, 1, "ensg", "ENSG00000141510")) == ["entrez:7157"]
    assert ids(look(runtime, 2, 1, "ensg", "ENSG00000141510", HUMAN)) == ["entrez:7157"]
    assert ids(look(runtime, 2, 1, "ensg", "ENSG00000141510", MOUSE)) == []
    # entities without a taxon are visible only unscoped
    assert ids(look(runtime, 2, 1, "ensg", "ENSG00000999999")) == ["entrez:555"]
    assert ids(look(runtime, 2, 1, "ensg", "ENSG00000999999", HUMAN)) == []
    assert ids(look(runtime, 2, 2, "entrez", "7157", MOUSE)) == []
    # products are scoped by their protein's taxon
    assert ids(look(runtime, 2, 1, "uniprot", "P02340", MOUSE)) == ["uniprot:P02340"]
    assert ids(look(runtime, 2, 1, "uniprot", "P02340", HUMAN)) == []


def test_symbols_exist_only_scoped(runtime):
    assert look(runtime, 2, 1, "genesymbol", "TP53")["candidates"] == []
    assert look(runtime, 2, 1, "genesymbol-syn", "TP53")["candidates"] == []
    assert ids(look(runtime, 2, 1, "genesymbol", "TP53", HUMAN)) == ["entrez:7157"]
    assert look(runtime, 2, 1, "genesymbol", "TP53", "7227")["candidates"] == []


def test_a_scoped_symbol_reads_only_the_records_of_its_taxon(runtime, monkeypatch):
    # TP53 is the symbol of the human gene 7157 and a synonym of the mouse gene 22059. The
    # records of the other taxon are left out when the rows are read, not after their
    # entities were looked up: a symbol is shared by hundreds of species in the real hubs.
    seen = []
    meta = runtime._meta
    monkeypatch.setattr(runtime, "_meta", lambda ids: seen.extend(ids) or meta(ids))
    assert ids(look(runtime, 2, 1, "genesymbol", "TP53", HUMAN)) == ["entrez:7157"]
    assert "entrez:7157" in seen and "entrez:22059" not in seen
    seen.clear()
    assert ids(look(runtime, 2, 1, "genesymbol-syn", "TP53", MOUSE)) == ["entrez:22059"]
    assert "entrez:22059" in seen and "entrez:7157" not in seen
    # The same symbol under two taxa in one call: each key gets the records of its own.
    human = key(2, 1, "genesymbol-syn", HUMAN, "TP53")
    mouse = key(2, 1, "genesymbol-syn", MOUSE, "TP53")
    both = runtime.lookup_many([human, mouse])
    assert ids(both[human]) == ids(look(runtime, 2, 1, "genesymbol-syn", "TP53", HUMAN))
    assert ids(both[mouse]) == ["entrez:22059"]


SYMBOLS = (
    "TP53", "P53", "TPX", "Trp53", "CALM1", "CALM2", "CALM3", "ABC1", "XYZ", "NOREC", "LINKED", "NONE",
)  # fmt: skip


def symbol_keys(taxa=(HUMAN, MOUSE, "7227")):
    return [
        key(2, 1, ns, taxon, symbol)
        for ns in ("genesymbol", "genesymbol-syn")
        for taxon in taxa
        for symbol in SYMBOLS
    ]


def test_a_scoped_symbol_reads_no_head_of_another_taxon(runtime, monkeypatch):
    # The kv stores keyed by taxon (format v3) read the rows of the scope and of the records
    # without a taxon; a store of the older format reads every species, so only the answer
    # is compared for it.
    heads = []
    for hub, kv in runtime.hub_kv.items():
        read = kv.heads
        monkeypatch.setattr(
            kv, "heads", lambda ids, hub=hub, read=read: heads.append((hub, set(ids))) or read(ids)
        )
    other = {HUMAN: {"22059", "P02340"}, MOUSE: {"7157", "P04637", "A0A0U1RQF1", "Q77777"}}
    for taxon, foreign in other.items():
        for symbol in ("TP53", "Trp53", "P53", "TPX", "NOREC"):
            for ns in ("genesymbol", "genesymbol-syn"):
                heads.clear()
                look(runtime, 2, 1, ns, symbol, taxon)
                for hub, local_ids in heads:
                    if runtime.hub_kv[hub].by_taxon:
                        assert not local_ids & foreign, (hub, taxon, symbol, local_ids & foreign)
    # A record without a taxon belongs to every scope: its row is read for each of them.
    for taxon, expected in ((HUMAN, ["entrez:4001"]), (MOUSE, [])):
        heads.clear()
        assert ids(look(runtime, 2, 1, "genesymbol", "LINKED", taxon)) == expected
        assert "Q55555" in {i for hub, local_ids in heads if hub == "uniprot" for i in local_ids}


def test_symbol_answers_are_those_of_format_v2(runtime, runtime_v2):
    keys = symbol_keys() + [
        key(2, 2, "entrez", HUMAN, "7157"),
        key(2, 1, "ensg", HUMAN, "ENSG00000141510"),
        key(2, 1, "uniprot", "", "P04637"),
        key(2, 1, "uniprot-sec", "", "Q88888"),
        key(1, 1, "chebi", "", "CHEBI:15377"),
    ]
    assert runtime.lookup_many(keys) == runtime_v2.lookup_many(keys)
    assert any(runtime.lookup_many(keys)[k]["candidates"] for k in symbol_keys())


def test_the_runtime_reads_both_kv_formats(snapshot):
    rt = IdentityRuntime(snapshot, cache_dir=snapshot.parent / "formats")
    try:
        formats = {hub: kv.manifest["format"] for hub, kv in rt.hub_kv.items()}
        assert set(formats.values()) <= set(identity_kv.HUB_KV_FORMATS)
        assert all(kv.by_taxon == (kv.manifest["format"] == identity_kv.HUB_KV_FORMAT) for kv in rt.hub_kv.values())
    finally:
        rt.close()


def test_an_unknown_kv_format_is_refused(snapshot, tmp_path):
    root = tmp_path / "tree"
    shutil.copytree(snapshot.parent.parent, root)
    manifest = root / "hubindex" / "chebi" / "0123456789ab" / "kv" / "manifest.json"
    data = json.loads(manifest.read_text())
    data["format"] = "omnipath-hub-kv-v9"
    manifest.write_text(json.dumps(data))
    with pytest.raises(KvMissing, match="omnipath-hub-kv-v9"):
        IdentityRuntime(root / "identity" / snapshot.name, cache_dir=tmp_path / "c")


def test_scoped_symbol_rows_are_read_in_rounds(runtime_v2, monkeypatch):
    # In the older format the answer does not depend on how many pairs a round reads.
    keys = symbol_keys((HUMAN,))
    whole = runtime_v2.lookup_many(keys)
    monkeypatch.setattr(type(runtime_v2), "_SCOPED_PAIRS", 1)
    assert runtime_v2.lookup_many(keys) == whole


def test_scoped_symbol_rows_are_read_once_in_format_v2(runtime_v2, monkeypatch):
    # The rows of a symbol for a taxon are kept: the second call reads no record again,
    # and answers as the first; a full memo is emptied and filled again.
    runtime = runtime_v2
    keys = symbol_keys((HUMAN, MOUSE))
    runtime._scoped_rows.clear()
    first = runtime.lookup_many(keys)
    assert runtime._scoped_rows
    reads = []
    rows = type(runtime)._rows
    monkeypatch.setattr(
        type(runtime), "_rows", staticmethod(lambda *a: reads.append(a[2]) or rows(*a))
    )
    assert runtime.lookup_many(keys) == first
    assert reads == []
    # A memo that holds three entries, read one pair at a time: it is emptied as it fills.
    monkeypatch.setattr(type(runtime), "_SCOPED_MEMO", 3)
    monkeypatch.setattr(type(runtime), "_SCOPED_PAIRS", 1)
    runtime._scoped_rows.clear()
    assert runtime.lookup_many(keys) == first
    # At most the limit and what one round adds (one pair under two taxa).
    assert reads and len(runtime._scoped_rows) <= 3 + 2
    # A key of another kind in the same call is not affected by the memo.
    mixed = keys + [key(2, 1, "ensg", HUMAN, "ENSG00000141510")]
    assert ids(runtime.lookup_many(mixed)[mixed[-1]]) == ids(
        look(runtime, 2, 1, "ensg", "ENSG00000141510", HUMAN)
    )


def test_format_v3_keeps_no_memo(runtime):
    runtime.lookup_many(symbol_keys())
    if all(kv.by_taxon for kv in runtime.hub_kv.values()):
        assert runtime._scoped_rows == {}


def test_gene_flag(runtime):
    expect = {
        (1, "ensg", "ENSG00000141510"): True,
        (1, "hgnc", "HGNC:11998"): True,
        (1, "refseq", "NM_000546"): True,
        (1, "genesymbol", "ABC1"): True,
        (1, "uniprot", "P04637"): False,
        (1, "uniprot_entry", "P53_HUMAN"): False,
        (2, "entrez", "7157"): True,  # route 2 is always a gene vote
    }
    for (route, ns, identifier), flag in expect.items():
        value = look(runtime, 2, route, ns, identifier, HUMAN if ns.startswith("gene") else "")
        assert value["candidates"] and value["gene"] is flag and value["products"] is False, ns
    assert look(runtime, 1, 1, "chebi", "CHEBI:15377")["gene"] is False


# --------------------------------------------------------- gene-level postings (rule 10)


def test_gene_level_postings_are_genes_from_genes_and_single_gene_proteins(runtime):
    # entrez record and linked protein both speak: one gene, not two candidates
    assert ids(look(runtime, 2, 1, "hgnc", "HGNC:11998")) == ["entrez:7157"]
    assert ids(look(runtime, 2, 1, "genesymbol", "TP53", HUMAN)) == ["entrez:7157"]
    # a protein linking three genes is evidence for none of them
    assert look(runtime, 2, 1, "hgnc", "HGNC:1442")["candidates"] == []
    assert ids(look(runtime, 2, 1, "genesymbol", "CALM1", HUMAN)) == ["entrez:801"]
    # a gene known only from gene_products is a valid target
    assert ids(look(runtime, 2, 1, "genesymbol", "NOREC", HUMAN)) == ["entrez:999"]
    gene = look(runtime, 2, 1, "genesymbol", "NOREC", HUMAN)["candidates"][0]
    assert gene[2:] == [3, None, False, False, ["999"]]


def test_more_than_ten_genes_are_all_kept(runtime):
    value = look(runtime, 2, 1, "ensg", "ENSG00000000001")
    assert sorted(ids(value)) == sorted(f"entrez:{3000 + n}" for n in range(12))
    nums = [c[0] for c in value["candidates"]]
    assert nums == sorted(nums)


def test_product_postings_are_proteins_with_their_gene_ids(runtime):
    (fact,) = look(runtime, 2, 1, "uniprot", "P0DP23")["candidates"]
    assert fact[1:4] == ["uniprot:P0DP23", 2, "uniprot:P0DP23"]
    assert fact[6] == ["801", "805", "808"]
    (fact,) = look(runtime, 2, 1, "uniprot", "Q99999")["candidates"]
    assert fact[6] == []


def test_isoform_projects_to_its_primary_parent(runtime):
    for ns, identifier in (("uniprot", "P04637-2"), ("ensp", "ENSP00000999999")):
        assert ids(look(runtime, 2, 1, ns, identifier)) == ["uniprot:P04637"]
    # no primary parent: the isoform entity stays as it is
    assert ids(look(runtime, 2, 1, "ensp", "ENSP00000888888")) == ["uniprot:Q11111-2"]
    (fact,) = look(runtime, 2, 1, "ensp", "ENSP00000888888")["candidates"]
    assert fact[3] == "uniprot:Q11111-2"  # the dash keeps it from being a product


def test_route_two_gene_identity(runtime):
    assert ids(look(runtime, 2, 2, "entrez", "7157")) == ["entrez:7157"]
    assert ids(look(runtime, 2, 2, "entrez", "999")) == ["entrez:999"]  # via gene_products only
    assert look(runtime, 2, 2, "entrez", "424242")["candidates"] == []
    # RaMP source genes map through the exceptions
    assert ids(look(runtime, 2, 2, "ramp_gene", "RAMP_G_1")) == ["entrez:7157"]
    assert ids(look(runtime, 2, 2, "ramp_gene", "RAMP_G_2")) == ["ramp_gene:RAMP_G_2"]
    assert look(runtime, 2, 2, "kegg_gene", "hsa:1")["candidates"] == []


# -------------------------------------------------- exceptions, quarantine, lipids


def test_exceptions_override_anchors_and_add_members(runtime):
    water = f"inchikey:{WATER}"
    # the record carries WATER's key but its exception makes it its own entity
    assert ids(look(runtime, 1, 1, "chembl", "CHEMBL999")) == ["chembl:CHEMBL999"]
    assert ids(look(runtime, 1, 1, "inchikey", WATER)) == [water]  # not polluted by the override
    assert ids(look(runtime, 1, 1, "pubchem", "7777")) == [water]  # attached by decision
    assert ids(look(runtime, 1, 1, "hmdb", "HMDB0070001")) == ["chebi:CHEBI:70001"]  # grouped
    assert ids(look(runtime, 1, 1, "chebi", "CHEBI:70001")) == ["chebi:CHEBI:70001"]


def test_quarantine(runtime):
    (bad,) = look(runtime, 1, 1, "chebi", "CHEBI:99999")["candidates"]
    assert bad[1:] == ["chebi:CHEBI:99999", 1, None, True, False, []]
    # a record claiming two keys creates neither key's entity
    assert look(runtime, 1, 1, "inchikey", TWO_A)["candidates"] == []


def test_lipid_names(runtime):
    (lipid,) = look(runtime, 1, 1, "goslin", "species:PE 36:2")["candidates"]
    assert lipid[1:4] == ["goslin:species:PE 36:2", 1, None]  # unanchored for the kernel
    # species level: every record carrying the name speaks
    assert sorted(ids(look(runtime, 1, 1, "goslin", "species:PC 34:1"))) == sorted(
        [f"inchikey:{LIPID}", "inchikey:AAAAAAAAAAAAAA-MMMMMMMMMM-N"]
    )
    # full structure: two records disagree on the key, so the name keys no structure ...
    assert look(runtime, 1, 1, "goslin", "full_structure:PC 16:0/18:1")["candidates"] == []
    # ... and a name with exactly one key resolves to it
    assert ids(look(runtime, 1, 1, "goslin", "full_structure:PE 18:0/18:1")) == [
        f"inchikey:{LIPID3}"
    ]


# --------------------------------------------------------------------- candidate facts


def test_candidate_facts_and_num(runtime):
    (fact,) = look(runtime, 2, 1, "uniprot", "P04637")["candidates"]
    digest = hashlib.sha256(b"uniprot:P04637").digest()
    assert fact == [
        int.from_bytes(digest[:8], "big") >> 1,
        "uniprot:P04637",
        2,
        "uniprot:P04637",
        False,
        True,
        ["7157"],
    ]
    (gene,) = look(runtime, 2, 2, "entrez", "7157")["candidates"]
    assert gene[1:] == ["entrez:7157", 3, None, False, False, ["7157"]]
    (chem,) = look(runtime, 1, 1, "chebi", "CHEBI:15377")["candidates"]
    assert chem[1:4] == [f"inchikey:{WATER}", 1, f"inchikey:{WATER}"]


def test_num_is_stable_unique_and_sorts_candidates(runtime):
    value = look(runtime, 2, 1, "uniprot", "Q88888")
    nums = [c[0] for c in value["candidates"]]
    assert nums == sorted(nums) and len(set(nums)) == len(nums) == 2
    k = key(2, 1, "uniprot", "", "Q88888")
    alone = runtime.lookup_many([k])[k]
    many = runtime.lookup_many([k, key(1, 1, "chebi", "", "CHEBI:15377")])[k]
    assert alone == many == value
    assert {c[1]: c[0] for c in value["candidates"]} == {
        eid: entity_num(eid) for eid in ("uniprot:P04637", "uniprot:P0DP23")
    }


def test_num_collision_in_a_batch_raises(runtime, monkeypatch):
    monkeypatch.setattr(identity_runtime, "entity_num", lambda eid: 7)
    ks = [key(1, 1, "chebi", "", "CHEBI:15377"), key(1, 1, "chebi", "", "CHEBI:15365")]
    with pytest.raises(ValueError, match="collision"):
        runtime.lookup_many(ks)


def test_no_candidate_cutoff(runtime):
    value = look(runtime, 1, 1, "cas", "0-0-0")
    assert len(value["candidates"]) == 12
    nums = [c[0] for c in value["candidates"]]
    assert nums == sorted(nums)


def test_lookup_and_record_wrappers(runtime):
    k = key(1, 1, "chebi", "", "CHEBI:15377")
    assert runtime.lookup(k) == runtime.lookup_many([k])[k]
    assert (
        runtime.record(f"inchikey:{WATER}")
        == runtime.record_many([f"inchikey:{WATER}"])[f"inchikey:{WATER}"]
    )


# ------------------------------------------------------------------------------ records


def rec(runtime, eid):
    return runtime.record(eid)


def pairs(record):
    return {tuple(p) for p in record["identifiers"]}


def test_record_shape_and_members(runtime):
    r = rec(runtime, "uniprot:P04637")
    assert set(r) == {"entity_id", "kind", "anchor", "taxon", "label", "identifiers", "gene_ids"}
    assert (r["entity_id"], r["kind"], r["anchor"], r["taxon"]) == (
        "uniprot:P04637", 2, "uniprot:P04637", "9606"
    )  # fmt: skip
    assert r["gene_ids"] == ["7157"]
    listed = [tuple(p) for p in r["identifiers"]]
    assert listed == sorted(listed)
    assert {
        ("uniprot", "P04637"),
        ("uniprot_entry", "P53_HUMAN"),
        ("genesymbol-syn", "P53"),
        ("uniprot-sec", "Q15086"),
        ("hgnc", "HGNC:11998"),
    } <= set(listed)
    chem = rec(runtime, f"inchikey:{WATER}")
    found = pairs(chem)
    assert ("inchikey", WATER) in found and ("chebi", "CHEBI:15377") in found
    assert ("name", "water") in found and ("hmdb", "HMDB0002111") in found
    assert ("pubchem", "7777") in found  # attached through exception_members
    assert ("chembl", "CHEMBL999") not in found  # overridden by its exception
    assert not any(ns in ("smiles", "synonym") for ns, _ in found)
    assert chem["kind"] == 1 and chem["anchor"] == f"inchikey:{WATER}" and chem["taxon"] is None
    grouped = pairs(rec(runtime, "chebi:CHEBI:70001"))
    assert {("chebi", "CHEBI:70001"), ("hmdb", "HMDB0070001")} <= grouped
    assert pairs(rec(runtime, "chembl:CHEMBL999")) >= {
        ("chembl", "CHEMBL999"),
        ("name", "overridden"),
    }


def test_gene_record_gets_aliases_only_from_single_gene_proteins(runtime):
    p53 = rec(runtime, "entrez:7157")
    assert (p53["kind"], p53["anchor"], p53["taxon"], p53["gene_ids"]) == (3, None, "9606", [])
    assert {
        ("ensg", "ENSG00000141510"),
        ("hgnc", "HGNC:11998"),
        ("genesymbol-syn", "P53"),
    } <= pairs(p53)
    calm = pairs(rec(runtime, "entrez:801"))
    assert ("genesymbol", "CALM1") in calm
    assert ("hgnc", "HGNC:1442") not in calm  # that protein links three genes
    # a gene known only from gene_products: no record rows, but linked protein aliases
    only = rec(runtime, "entrez:999")
    assert only["taxon"] == "9606" and ("genesymbol", "NOREC") in pairs(only)
    assert only["label"] == "999"


def test_record_taxon_extras_and_missing(runtime):
    assert rec(runtime, "ramp_gene:RAMP_G_2")["taxon"] == "9606"
    assert rec(runtime, "ramp_gene:RAMP_G_2")["kind"] == 3
    assert rec(runtime, "entrez:555")["taxon"] is None
    assert rec(runtime, "chebi:CHEBI:99999")["kind"] == 1
    for missing in ("inchikey:" + "Z" * 14 + "-UHFFFAOYSA-N", "entrez:424242", "uniprot:ZZZZZZ"):
        with pytest.raises(ValueError, match="absent entity"):
            runtime.record_many([missing])


@pytest.mark.parametrize(
    "entity, label",
    [
        # chemical: ChEBI > HMDB > ChEMBL > PubChem > other hub > systematic
        (f"inchikey:{WATER}", "water"),  # chebi: shortest of two, ties alphabetical
        (f"inchikey:{ASPIRIN}", "ASPIRIN"),  # chebi name is an InChIKey, hmdb's is over 80
        (f"inchikey:{ONLY_PC}", "pc-name"),  # pubchem before an unlisted hub
        (f"inchikey:{OTHERS}", "a"),  # other hubs: shorter then alphabetical
        (f"inchikey:{SYSTEMATIC}", "sys-name"),  # systematic name last
        (f"inchikey:{NAMELESS}", NAMELESS),  # id local part, no CHEBI:/CID: formatting
        (f"inchikey:{LIPID}", "PC 16:0/18:1"),  # most specific Goslin name beats a chemical name
        ("goslin:species:PE 36:2", "PE 36:2"),
        # a fatty-acyl shorthand only where a lipid database holds the compound
        (f"inchikey:{FATTY}", "FA 18:1(9Z)"),
        (f"inchikey:{AMINO}", "L-glutamine"),
        # an "-oate" misread as a wax ester without its alcohol never labels
        (f"inchikey:{ANION}", "FA 6:0"),
        ("chebi:CHEBI:99999", "multi"),
        ("chebi:CHEBI:70001", "grouped-chebi"),
        # gene: NCBI symbol, else the id
        ("entrez:7157", "TP53"),
        ("entrez:555", "555"),
        # protein: primary gene name > entry name > accession
        ("uniprot:P04637", "TP53"),
        ("uniprot:A0A0U1RQF1", "A0A0U1RQF1_HUMAN"),  # only a genesymbol-syn: never a label
        ("uniprot:Q99999", "Q99999"),
        ("uniprot:P04637-2", "P04637-2"),
    ],
)
def test_labels(runtime, entity, label):
    assert rec(runtime, entity)["label"] == label


def test_cache_hit_path(snapshot, tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    first = IdentityRuntime(snapshot, cache_dir=cache)
    try:
        wanted = [f"inchikey:{WATER}", "uniprot:P04637", "entrez:7157"]
        built = first.record_many(wanted)
    finally:
        first.close()
    database = cache / FINGERPRINT / "entities.sqlite"
    assert database.is_file()
    con = sqlite3.connect(database)
    try:
        assert con.execute("PRAGMA journal_mode").fetchone()[0] == "wal"
        assert con.execute("SELECT count(*) FROM records").fetchone()[0] == 3
    finally:
        con.close()

    def boom(self, entity_ids):
        raise AssertionError("cache miss: " + ",".join(entity_ids))

    monkeypatch.setattr(IdentityRuntime, "_build_records", boom)
    second = IdentityRuntime(snapshot, cache_dir=cache)
    try:
        assert second.record_many(wanted) == built  # served from sqlite
        with pytest.raises(AssertionError):
            second.record_many([f"inchikey:{ASPIRIN}"])  # a genuine miss still builds
    finally:
        second.close()


def test_cache_location_default_and_env(snapshot, tmp_path, monkeypatch):
    monkeypatch.delenv("OMNIPATH_IDENTITY_CACHE", raising=False)
    assert identity_runtime.default_cache_dir(snapshot) == snapshot.resolve().parent / "cache"
    monkeypatch.setenv("OMNIPATH_IDENTITY_CACHE", str(tmp_path / "elsewhere"))
    assert identity_runtime.default_cache_dir(snapshot) == tmp_path / "elsewhere"


# ----------------------------------------------------------------------------- rule 11


def chemical(identifiers, **kw):
    return RawEntityObservation(
        "e",
        "chemical_entity",
        kw.pop("namespace", "name"),
        kw.pop("identifier", "x"),
        identifiers=[dict(ns=ns, id=i) for ns, i in identifiers],
        **kw,
    )


@pytest.fixture
def structures(tmp_path, monkeypatch):
    from omnipath_resolver.canonical import structures as module

    monkeypatch.setenv("OMNIPATH_STRUCTURE_CACHE", str(tmp_path / "structures"))
    monkeypatch.setenv("OMNIPATH_GOSLIN_CACHE", str(tmp_path / "goslin"))
    module._cached.cache_clear()
    yield module
    module._cached.cache_clear()


def anchors(obs):
    votes, observed = votes_for(obs, CHEMICAL_POLICY)
    _, rows = observation_bundle("e", obs, votes, observed, "chemical")
    return {r["identifier"]: r["anchor"] for r in rows if r["ns"] == "inchikey"}


def test_derived_inchikey_is_not_an_anchor_but_a_stated_one_is(structures):
    assert anchors(chemical([("smiles", "CCO")])) == {ETHANOL: ""}
    assert anchors(chemical([("inchikey", ETHANOL)])) == {ETHANOL: ETHANOL}
    # a stated key suppresses derivation entirely
    assert anchors(chemical([("smiles", "CCO"), ("inchikey", WATER)])) == {WATER: WATER}


# ------------------------------------------------------------------------ end to end


@pytest.fixture
def matcher(snapshot, tmp_path, monkeypatch, structures):
    monkeypatch.setenv("OMNIPATH_IDENTITY_CACHE", str(tmp_path / "cache"))
    m = LibraryMatcher(snapshot)
    yield m
    m.close()


def run(matcher, **observations):
    return matcher.targets(observations)


def test_matcher_chooses_identity_runtime(matcher):
    assert isinstance(matcher.runtime, IdentityRuntime)
    assert matcher.libraries == ["gene_protein", "chemical", "reaction"]


def test_compact_references_are_rejected(tmp_path):
    (tmp_path / "manifest.json").write_text('{"format": "omnipath-full-two-index-v1"}')
    with pytest.raises(ValueError, match="identity directory is required"):
        open_runtime(tmp_path)


def test_chemical_stated_vs_derived_inchikey(matcher):
    water = f"inchikey:{WATER}"
    out = run(
        matcher,
        # a stated InChIKey decides alone, overriding the conflicting ChEBI id
        stated=chemical([("inchikey", ETHANOL)], namespace="chebi", identifier="CHEBI:15377"),
        # the same structure derived from SMILES is an ordinary vote: it conflicts with ChEBI,
        # and the source's own structure then decides
        derived=chemical([("smiles", "CCO")], namespace="chebi", identifier="CHEBI:15377"),
        # alone, a derived key still resolves
        alone=chemical([("smiles", "CCO")]),
        agreeing=chemical([("smiles", "O")], namespace="chebi", identifier="CHEBI:15377"),
        plain=chemical([], namespace="chebi", identifier="CHEBI:15377"),
        attached=chemical([], namespace="pubchem", identifier="7777"),
        override=chemical([], namespace="chembl", identifier="CHEMBL999"),
        quarantined=chemical([], namespace="chebi", identifier="CHEBI:99999"),
    )
    (stated,) = out["stated"]
    assert stated.matched and stated.node_id == f"inchikey:{ETHANOL}" and stated.label == "ethanol"
    (derived,) = out["derived"]
    assert derived.node_id == f"inchikey:{ETHANOL}"
    (alone,) = out["alone"]
    assert alone.node_id == f"inchikey:{ETHANOL}"
    (agreeing,) = out["agreeing"]
    assert agreeing.node_id == water and agreeing.label == "water"
    assert out["plain"][0].node_id == water
    assert "CHEBI:15377" in out["plain"][0].aliases["chebi"]
    assert out["attached"][0].node_id == water
    assert out["override"][0].node_id == "chembl:CHEMBL999"
    assert not out["quarantined"][0].matched  # a quarantined entity is never accepted


def test_chemical_ambiguity_has_no_cutoff(matcher):
    (match,) = run(matcher, many=chemical([], namespace="cas", identifier="0-0-0"))["many"]
    assert not match.matched
    k = key(1, 1, "cas", "", "0-0-0")
    assert len(matcher.runtime.lookup(k)["candidates"]) == 12
    resolved, metrics = matcher.runtime.resolve(
        [dict(input_id="q", target=1)],
        [
            dict(
                input_id="q", ns="cas", identifier="0-0-0", scope="", anchor="", target=1,
                route=1, ordinal=0, lookup_key=k,
            )
        ],
    )  # fmt: skip
    assert resolved["results"][0]["outcome"] == "Ambiguous"
    assert resolved["results"][0]["candidate_count"] == 12
    assert resolved["results"][0]["entities"] == []
    assert set(metrics) == {
        "lookup_seconds", "decision_seconds", "entity_fetch_seconds", "total_seconds",
        "unique_keys", "candidate_records", "entity_records",
    }  # fmt: skip


def protein(identifier, identifiers=(), taxon=None, entity_type="protein", namespace="uniprot"):
    return RawEntityObservation(
        "p",
        entity_type,
        namespace,
        identifier,
        taxon=taxon,
        identifiers=[dict(ns=ns, id=i) for ns, i in identifiers],
    )


def test_protein_with_geneid_resolves_to_gene_with_product(matcher):
    out = run(
        matcher,
        both=protein("P04637", [("entrez", "7157")], taxon=HUMAN),
        secondary=protein("Q15086"),
        primary_wins=protein("Q99999"),
        shared=protein("Q88888"),
        conflict=protein("P04637", [("entrez", "801")]),
        isoform=protein("P04637-2"),
    )
    (both,) = out["both"]
    assert both.node_id == "entrez:7157" and both.gene_mapping_status == "resolved"
    assert both.protein_identifier == "P04637" and both.protein_label == "TP53"
    assert both.protein_gene_candidates == ("entrez:7157",)
    assert both.reference_library == "gene_protein" and both.taxon == HUMAN
    (secondary,) = out["secondary"]
    assert secondary.node_id == "entrez:7157" and secondary.protein_identifier == "P04637"
    (primary_wins,) = out["primary_wins"]
    assert primary_wins.protein_identifier == "Q99999" and primary_wins.matched
    (shared,) = out["shared"]
    assert not shared.matched  # two entries claim Q88888 as secondary: abstain
    (conflict,) = out["conflict"]
    assert conflict.gene_mapping_status == "conflict" and conflict.node_id == "uniprot:P04637"
    assert out["isoform"][0].protein_node_id == "uniprot:P04637"


def test_gene_symbol_with_taxon(matcher):
    out = run(
        matcher,
        human=protein("TP53", taxon=HUMAN, entity_type="gene", namespace="genesymbol"),
        mouse=protein("TP53", taxon=MOUSE, entity_type="gene", namespace="genesymbol"),
        nothing=protein("TP53", taxon="7227", entity_type="gene", namespace="genesymbol"),
        no_taxon=protein("TP53", entity_type="gene", namespace="genesymbol"),
        scoped=protein("P04637", [("genesymbol", "TP53")], taxon=HUMAN),
        wrong_scope=protein("P04637", [("genesymbol", "TP53")], taxon=MOUSE),
        by_ensg=protein("ENSG00000141510", entity_type="gene", namespace="ensembl"),
        many=protein("ENSG00000000001", entity_type="gene", namespace="ensembl"),
        ramp=protein("RAMP_G_1", entity_type="gene", namespace="ramp"),
    )
    assert [m.node_id for m in out["human"]] == ["entrez:7157"]
    assert out["human"][0].label == "TP53" and out["human"][0].taxon == HUMAN
    assert [m.node_id for m in out["mouse"]] == ["entrez:22059"]
    assert out["mouse"][0].label == "Trp53"
    assert not out["nothing"][0].matched
    assert not out["no_taxon"][0].matched  # symbols never match unscoped
    assert out["scoped"][0].node_id == "entrez:7157"
    assert out["scoped"][0].protein_identifier == "P04637"
    # the symbol scopes the whole observation: the human accession is invisible under mouse
    # (a miss never vetoes), so only the mouse synonym speaks
    assert out["wrong_scope"][0].node_id == "entrez:22059"
    assert out["by_ensg"][0].node_id == "entrez:7157"
    assert not out["many"][0].matched  # twelve genes share the id: ambiguous, no cutoff
    assert out["ramp"][0].node_id == "entrez:7157"


# ------------------------------------------------------ records pointing to several anchors


def test_a_record_with_candidates_votes_for_each_of_them(runtime):
    # pubchem:P2 is quarantined and claims two InChIKeys: it contributes both, not itself
    out = look(runtime, 1, 1, "pubchem", "P2")
    assert set(ids(out)) == {f"inchikey:{TWO_C}", f"inchikey:{TWO_D}"} and len(ids(out)) == 2
    assert not any(c[4] for c in out["candidates"])  # the anchors are clean entities
    assert all(c[3] == c[1] for c in out["candidates"])
    # ...but it is never a member of those entities
    members = runtime.record(f"inchikey:{TWO_C}")["identifiers"]
    assert ["pubchem", "P2"] not in members and ["hmdb", "HMDB0088888"] in members
    # a record without candidates is unchanged
    assert ids(look(runtime, 1, 1, "chebi", "CHEBI:99999")) == ["chebi:CHEBI:99999"]


def test_candidates_intersect_with_a_second_identifier(matcher):
    alone = chemical([], namespace="pubchem", identifier="P2")
    agreeing = chemical([("hmdb", "HMDB0088888")], namespace="pubchem", identifier="P2")
    out = run(matcher, alone=alone, agreeing=agreeing)
    assert not out["alone"][0].matched  # two anchors: Ambiguous
    resolved, _ = matcher.runtime.resolve(
        [dict(input_id="q", target=1)],
        [
            dict(
                input_id="q", ns="pubchem", identifier="P2", scope="", anchor="", target=1,
                route=1, ordinal=0, lookup_key=key(1, 1, "pubchem", "", "P2"),
            )
        ],
    )  # fmt: skip
    assert resolved["results"][0]["outcome"] == "Ambiguous"
    assert resolved["results"][0]["candidate_count"] == 2
    (agreeing,) = out["agreeing"]
    assert agreeing.matched and agreeing.node_id == f"inchikey:{TWO_C}"


def test_coarse_cross_references_do_not_veto_specific_identifiers(matcher):
    water, ethanol = f"inchikey:{WATER}", f"inchikey:{ETHANOL}"
    out = run(
        matcher,
        # bigg:x1 is ethanol: a secondary BiGG id cannot veto the ChEBI id
        vetoed=chemical([("bigg", "x1")], namespace="chebi", identifier="CHEBI:15377"),
        agreeing=chemical([("kegg", "C00001")], namespace="chebi", identifier="CHEBI:15377"),
        # alone it still resolves
        alone=chemical([("bigg", "x1")]),
        # the source's own identifier is no cross-reference: the conflict stands
        primary=chemical([("chebi", "CHEBI:15377")], namespace="bigg", identifier="x1"),
        # specific identifiers that disagree among themselves still abstain
        specific=chemical([("smiles", "CCO"), ("bigg", "x1")], namespace="chebi", identifier="CHEBI:15377"),
    )  # fmt: skip
    assert out["vetoed"][0].node_id == water
    assert out["agreeing"][0].node_id == water
    assert out["alone"][0].node_id == ethanol
    assert not out["primary"][0].matched
    # the conflict among specific identifiers is decided by the source's structure
    assert out["specific"][0].node_id == ethanol


def test_conflicting_structures_and_protonation_states(matcher):
    nadh, nadh_2 = f"inchikey:{NADH}", f"inchikey:{NADH_2}"
    out = run(
        matcher,
        # two SMILES that disagree are no structure to decide by
        structures=chemical([("smiles", "CCO"), ("smiles", "O")], namespace="chebi", identifier="CHEBI:15365"),
        # one molecule in two protonation states: the neutral form decides, whichever is primary
        neutral=chemical([("chebi", "CHEBI:80002")], namespace="chebi", identifier="CHEBI:80001"),
        charged=chemical([("chebi", "CHEBI:80001")], namespace="chebi", identifier="CHEBI:80002"),
        # without the neutral form among them, the primary id decides
        no_neutral=chemical([("chebi", "CHEBI:80003")], namespace="chebi", identifier="CHEBI:80002"),
        # a coarse cross-reference naming another molecule does not block the neutral form
        coarse=chemical([("chebi", "CHEBI:80002"), ("bigg", "x1")], namespace="chebi", identifier="CHEBI:80001"),
        # different molecules are no protonation states: still a conflict
        different=chemical([("chebi", "CHEBI:16236")], namespace="chebi", identifier="CHEBI:80001"),
    )  # fmt: skip
    assert not out["structures"][0].matched
    assert out["neutral"][0].node_id == nadh
    assert out["charged"][0].node_id == nadh
    assert out["no_neutral"][0].node_id == nadh_2
    assert out["coarse"][0].node_id == nadh
    assert not out["different"][0].matched


def test_gene_label_falls_back_to_gene_info_name():
    from omnipath_resolver.identity_runtime import choose_label

    # A tRNA gene has no UniProt product, so only gene_info's name row carries its symbol.
    assert choose_label(3, "entrez:100189199", [("entrez", "name", "TRL-CAA6-1")]) == "TRL-CAA6-1"
    # A product's symbol still wins.
    rows = [("entrez", "genesymbol", "TP53"), ("entrez", "name", "tumor protein")]
    assert choose_label(3, "entrez:7157", rows) == "TP53"
    assert choose_label(3, "entrez:42", []) == "42"


def test_workers_opening_a_new_cache_together(runtime):
    import threading

    # Switching a new cache to WAL needs an exclusive lock, and SQLite answers "locked"
    # at once instead of waiting when several workers try it together.
    path = runtime.cache_dir / str(runtime.fingerprint) / "entities.sqlite"
    errors = []
    for _ in range(10):
        runtime._cache.close()
        for suffix in ("", "-wal", "-shm"):
            path.with_name(path.name + suffix).unlink(missing_ok=True)
        start = threading.Barrier(8)

        def open_cache():
            start.wait()
            try:
                runtime._open_cache().close()
            except sqlite3.OperationalError as exc:
                errors.append(exc)

        threads = [threading.Thread(target=open_cache) for _ in range(8)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        runtime._cache = runtime._open_cache()
    assert errors == []
