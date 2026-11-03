"""Hub indexes and identity decisions on tiny fixture hubs: one case per rule, exact rows."""

from __future__ import annotations

import hashlib
import json

import duckdb
import pyarrow.parquet as pq
import pytest

from identity_fixtures import key, write_hubs
from omnipath_build.identity import build_hub_index, build_identity
from omnipath_build.identity.common import HUBS

KA, KB, KC, KD, KE, KF, KG, KH = (key(x) for x in "ABCDEFGH")
ARGS = dict(memory="1GB", threads=2, min_free_gib=0.01)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("identity")
    hubs = write_hubs(root / "hubs")
    manifests = {
        h: build_hub_index(h, hubs, root / "hubindex", goslin_cache=root / "goslin", **ARGS)
        for h in HUBS
        if (hubs / f"{h}.parquet").exists()
    }
    snapshot = build_identity(root / "hubindex", root / "identity", **ARGS)
    return root, manifests, snapshot


@pytest.fixture(scope="module")
def q(built):
    root, manifests, snapshot = built
    out = root / "identity" / snapshot["fingerprint"]
    c = duckdb.connect()

    def index(hub):
        return next((root / "hubindex" / hub).glob("*/manifest.json")).parent

    def query(sql):
        for hub in HUBS:
            if (root / "hubindex" / hub).exists():
                d = index(hub)
                sql = (
                    sql.replace(f"${hub}.records", f"'{d}/records.parquet'")
                    .replace(f"${hub}.by_id", f"read_parquet('{d}/by_id/*/*.parquet')")
                    .replace(f"${hub}.by_record", f"read_parquet('{d}/by_record/*/*.parquet')")
                    .replace(f"${hub}.xrefs", f"'{d}/xrefs.parquet'")
                )
        for name in ("exceptions", "exception_members", "entities_extra", "lipid_structures"):
            sql = sql.replace(f"${name}", f"'{out}/{name}.parquet'")
        sql = sql.replace("$products_by_gene", f"'{out}/gene_products_by_gene.parquet'")
        sql = sql.replace("$products_by_protein", f"'{out}/gene_products_by_protein.parquet'")
        return c.execute(sql).fetchall()

    query.out = out
    return query


# ------------------------------------------------------------------- hub index
def test_hub_index_records(q):
    assert q("SELECT * FROM $chebi.records WHERE local_id IN ('CHEBI:1','CHEBI:5')") == [
        ("CHEBI:1", "chebi:CHEBI:1", None, "inchikey:" + KA, 1, False),
        ("CHEBI:5", "chebi:CHEBI:5", None, "inchikey:" + KE, 1, False),
    ]
    # Two InChIKeys: no anchor, anchor_count 2. The hub id HMDB2 is normalized to HMDB0000002.
    assert q("SELECT * FROM $hmdb.records ORDER BY local_id") == [
        ("HMDB0000002", "hmdb:HMDB0000002", None, None, 2, False),
        ("HMDB0000004", "hmdb:HMDB0000004", None, "inchikey:" + KD, 1, False),
    ]
    # The excluded "empty" InChIKey is no anchor; pubchem 100 shares chebi 1's key.
    assert q("SELECT local_id,anchor,anchor_count FROM $pubchem.records ORDER BY local_id") == [
        ("100", "inchikey:" + KA, 1),
        ("101", None, 0),
    ]
    assert q("SELECT * FROM $uniprot.records WHERE local_id IN ('P04637','P04637-2','A0A000')") == [
        ("A0A000", "uniprot:A0A000", "9606", "uniprot:A0A000", 1, False),
        ("P04637", "uniprot:P04637", "9606", "uniprot:P04637", 1, True),
        ("P04637-2", "uniprot:P04637-2", "9606", "uniprot:P04637-2", 1, False),
    ]
    # An NCBI Gene record is not anchored in the index (its entity is `entrez:<local_id>`).
    assert q("SELECT * FROM $entrez.records WHERE local_id='9999'") == [
        ("9999", "entrez:9999", "10090", None, 0, False)
    ]
    assert q("SELECT local_id FROM $chebi.records") == sorted(
        q("SELECT local_id FROM $chebi.records")
    )


def test_hub_index_by_id_tags_and_expansions(q):
    cols = "ns,identifier,tag,taxon,anchor,anchor_count"
    p = "uniprot:P04637"
    rows = q(
        f"SELECT {cols} FROM $uniprot.by_id WHERE local_id='P04637' ORDER BY ns,identifier,tag"
    )
    assert rows == [
        ("ensg", "ENSG00000141510", "claim", "9606", p, 1),
        ("ensp", "ENSP00000269305", "claim", "9606", p, 1),
        ("entrez", "7157", "claim", "9606", p, 1),
        ("genbank", "AAA12345", "version_stripped", "9606", p, 1),
        ("genbank", "AAA12345.1", "claim", "9606", p, 1),
        ("genesymbol", "P53", "symbol_synonym", "9606", p, 1),
        ("genesymbol", "TP53", "claim", "9606", p, 1),
        ("genesymbol-syn", "P53", "claim", "9606", p, 1),
        ("hgnc", "HGNC:11998", "claim", "9606", p, 1),
        ("refseq_protein", "NP_000537", "version_stripped", "9606", p, 1),
        ("refseq_protein", "NP_000537.3", "claim", "9606", p, 1),
        ("uniprot", "P04637", "native", "9606", p, 1),
        ("uniprot", "Q15086", "secondary", "9606", p, 1),
        ("uniprot", "Q16535", "secondary", "9606", p, 1),
        ("uniprot-sec", "Q15086", "claim", "9606", p, 1),
        ("uniprot-sec", "Q16535", "claim", "9606", p, 1),
        ("uniprot_entry", "P53_HUMAN", "claim", "9606", p, 1),
    ]
    # Chemical: native id, cross-references and aliases; names are not lookup keys.
    assert q(
        f"SELECT {cols} FROM $chebi.by_id WHERE local_id='CHEBI:1' ORDER BY ns,identifier"
    ) == [
        ("cas", "50-00-0", "claim", None, "inchikey:" + KA, 1),
        ("chebi", "CHEBI:1", "native", None, "inchikey:" + KA, 1),
        ("chebi", "CHEBI:99", "claim", None, "inchikey:" + KA, 1),
        ("inchikey", KA, "claim", None, "inchikey:" + KA, 1),
        ("kegg", "C00001", "claim", None, "inchikey:" + KA, 1),
    ]
    assert q(
        "SELECT ns,tag,anchor,anchor_count FROM $hmdb.by_id WHERE ns='inchikey' AND local_id='HMDB0000002' ORDER BY identifier"
    ) == [
        ("inchikey", "claim", None, 2),
        ("inchikey", "claim", None, 2),
    ]


def test_hub_index_goslin_rows(q):
    assert q(
        "SELECT ns,identifier,local_id,tag,anchor_count FROM $swisslipids.by_id WHERE ns='goslin' ORDER BY identifier"
    ) == [
        ("goslin", "full_structure:FA 18:0;5Me", "SLM:5", "claim", 0),
        ("goslin", "sn_position:PC 16:0/18:1", "SLM:1", "claim", 0),
    ]
    # Only the most specific level of a record is kept (PC 34:1 is species level).
    assert q(
        "SELECT source_type,value FROM $swisslipids.by_record WHERE local_id='SLM:1' ORDER BY 1,2"
    ) == [
        ("goslin", "sn_position:PC 16:0/18:1"),
        ("lipid_shorthand", "PC 16:0/18:1"),
        ("name", "PC 34:1"),
        ("swisslipids", "SLM:1"),
    ]
    assert q("SELECT identifier FROM $chebi.by_id WHERE ns='goslin' AND local_id='CHEBI:7'") == [
        ("species:PC 34:1",)
    ]
    assert q(
        "SELECT identifier FROM $lipidmaps.by_id WHERE ns='goslin' AND local_id='LMQ' ORDER BY 1"
    ) == [
        ("sn_position:PE 16:0/18:1",),
        ("sn_position:PE 18:1/16:0",),
    ]


def test_hub_index_by_record_and_xrefs(q):
    assert q(
        "SELECT source_type,value FROM $chebi.by_record WHERE local_id='CHEBI:1' ORDER BY 1,2"
    ) == [
        ("cas", "50-00-0"),
        ("chebi", "CHEBI:1"),
        ("chebi", "CHEBI:99"),
        ("formula", "C2H6"),
        ("inchikey", KA),
        ("kegg", "C00001"),
        ("name", "Alpha"),
        ("smiles", "CC"),
        ("synonym", "alpha synonym"),
    ]
    assert q("SELECT * FROM $chembl.xrefs") == [("CHEMBL1", "chebi", "CHEBI:1")]
    assert q("SELECT * FROM $chebi.xrefs ORDER BY 1,2,3") == [
        ("CHEBI:1", "chebi", "CHEBI:99"),
        ("CHEBI:3", "metanetx", "MNXM3"),
    ]
    assert q("SELECT ns,identifier FROM $uniprot.xrefs WHERE local_id='P04637' ORDER BY 1") == [
        ("entrez", "7157")
    ]


def test_hub_index_layout_manifest_and_sorting(built):
    root, manifests, _ = built
    m = manifests["uniprot"]
    assert set(m) >= {"hub", "input", "code_sha256", "counts", "timings"} and m["hub"] == "uniprot"
    assert len(m["input"]["sha256"]) == 64 and m["input"]["bytes"] > 0
    d = root / "hubindex" / "uniprot" / m["input"]["sha256"][:12]
    assert (d / "manifest.json").is_file() and not list(
        (root / "hubindex" / "uniprot").glob("*.building")
    )
    for path in (d / "by_id").glob("part=*/data.parquet"):
        t = pq.read_table(path)
        keys = list(
            zip(
                t.column("ns").to_pylist(),
                t.column("identifier").to_pylist(),
                t.column("local_id").to_pylist(),
            )
        )
        assert keys == sorted(keys)
        assert {
            hashlib.md5(i.encode()).hexdigest()[:2] for i in t.column("identifier").to_pylist()
        } == {path.parent.name[5:]}
        assert pq.ParquetFile(path).metadata.row_group(0).column(0).statistics is not None
    for path in (d / "by_record").glob("part=*/data.parquet"):
        t = pq.read_table(path)
        ids = t.column("local_id").to_pylist()
        assert ids == sorted(ids)
        assert {hashlib.md5(f"uniprot:{i}".encode()).hexdigest()[:2] for i in ids} == {
            path.parent.name[5:]
        }


def test_hub_index_is_reused_unchanged(built):
    root, manifests, _ = built
    again = build_hub_index(
        "chebi", root / "hubs", root / "hubindex", goslin_cache=root / "goslin", **ARGS
    )
    assert again == manifests["chebi"]


# ------------------------------------------------------------------- decisions
EXPECTED_ENTITY_DEFINING = [
    # (record_id, entity_id, decision)
    ("bigg:b1", "bigg:b1", "ambiguous_native"),  # reaches two different anchors
    ("bigg:b2", "bigg:b2", "ambiguous_native"),  # reaches a quarantined record
    ("bigg:b3", "bigg:b3", "structureless"),
    ("bigg:g10", "metanetx:MNXM10", "grouped"),  # preferred hub of the group: metanetx
    ("chembl:CHEMBL1", "inchikey:" + KA, "attached"),  # rule 3, outgoing edge
    ("hmdb:HMDB0000002", "hmdb:HMDB0000002", "quarantined"),  # rule 2: two InChIKeys
    ("lipidmaps:LM1", "goslin:sn_position:PC 16:0/18:1", "lipid_name"),  # rule 1 lipid name
    ("lipidmaps:LMQ", "lipidmaps:LMQ", "quarantined"),  # rule 2: two lipid names
    ("metanetx:MNXM10", "metanetx:MNXM10", "grouped"),
    ("metanetx:MNXM20", "metanetx:MNXM20", "ambiguous_native"),  # rule 4: two metanetx in a group
    ("metanetx:MNXM21", "metanetx:MNXM21", "ambiguous_native"),
    ("metanetx:MNXM3", "inchikey:" + KD, "attached"),  # rule 3, incoming edge
    ("metanetx:MNXM4", "metanetx:MNXM4", "structureless"),  # chains are not followed
    ("pubchem:101", "pubchem:101", "structureless"),
    ("ramp:RAMP_C_1", "ramp:RAMP_C_1", "ambiguous_native"),
    ("swisslipids:SLM:1", "goslin:sn_position:PC 16:0/18:1", "lipid_name"),
    ("swisslipids:SLM:5", "goslin:full_structure:FA 18:0;5Me", "lipid_name"),
    # Reactions: a MetaNetX reaction attaches to the one Rhea master it names; without one it is
    # its own entity; naming two masters it is ambiguous (and evidence for both, record_candidates).
    ("metanetx_reaction:MNXR1", "rhea:10000", "attached"),
    ("metanetx_reaction:MNXR2", "metanetx_reaction:MNXR2", "structureless"),
    ("metanetx_reaction:MNXR3", "metanetx_reaction:MNXR3", "structureless"),
    ("metanetx_reaction:MNXR4", "metanetx_reaction:MNXR4", "ambiguous_native"),
]
ROUTE2 = ("explicit_source_gene", "source_gene_only")


def test_exceptions_are_the_non_trivial_decisions(q):
    rows = q(
        "SELECT record_id,entity_id,decision FROM $exceptions WHERE NOT starts_with(record_id,'ramp_gene:') ORDER BY 1,2,3"
    )
    assert rows == sorted(EXPECTED_ENTITY_DEFINING)
    assert q(
        "SELECT record_id,entity_id,decision,quarantined FROM $exceptions WHERE decision='quarantined' ORDER BY 1"
    ) == [
        ("hmdb:HMDB0000002", "hmdb:HMDB0000002", "quarantined", True),
        ("lipidmaps:LMQ", "lipidmaps:LMQ", "quarantined", True),
    ]
    assert q("SELECT count(*) FROM $exceptions WHERE quarantined") == [(2,)]


def test_anchored_records_are_never_exceptions_unless_quarantined(q):
    anchored = q(
        "SELECT record_id FROM $chebi.records WHERE anchor_count=1 UNION SELECT record_id FROM $pubchem.records WHERE anchor_count=1 "
        "UNION SELECT record_id FROM $hmdb.records WHERE anchor_count=1 UNION SELECT record_id FROM $lipidmaps.records WHERE anchor_count=1 "
        "UNION SELECT record_id FROM $uniprot.records WHERE anchor_count=1"
    )
    assert len(anchored) > 20
    exceptions = {r[0] for r in q("SELECT record_id FROM $exceptions")}
    assert not exceptions & {r[0] for r in anchored}
    # No entrez record is an exception either (its entity is entrez:<local_id>).
    assert not {r for r in exceptions if r.startswith("entrez:")}
    assert q(
        "SELECT count(*) FROM $exceptions WHERE decision='quarantined' AND record_id IN (SELECT record_id FROM $hmdb.records WHERE anchor_count=2)"
    ) == [(1,)]


def test_exception_members_and_entities_extra(q):
    assert q("SELECT entity_id,record_id FROM $exception_members ORDER BY 1,2") == sorted(
        [(e, r) for r, e, _ in EXPECTED_ENTITY_DEFINING]
        + [(e, r) for r, e, _ in EXPECTED_RAMP_GENES]
    )
    assert q("SELECT * FROM $entities_extra ORDER BY entity_id") == [
        ("bigg:b1", "chemical", None, False, "bigg:b1"),
        ("bigg:b2", "chemical", None, False, "bigg:b2"),
        ("bigg:b3", "chemical", None, False, "bigg:b3"),
        ("goslin:full_structure:FA 18:0;5Me", "chemical", None, False, "swisslipids:SLM:5"),
        ("goslin:sn_position:PC 16:0/18:1", "chemical", None, False, "lipidmaps:LM1"),
        # Lipid names only a quarantined record claims: entities so they can be candidates.
        ("goslin:sn_position:PE 16:0/18:1", "chemical", None, False, None),
        ("goslin:sn_position:PE 18:1/16:0", "chemical", None, False, None),
        ("hmdb:HMDB0000002", "chemical", None, True, "hmdb:HMDB0000002"),
        # The keys of a quarantined record still get entities (no record is anchored on them).
        ("inchikey:" + KB, "chemical", None, False, None),
        ("inchikey:" + KC, "chemical", None, False, None),
        ("lipidmaps:LMQ", "chemical", None, True, "lipidmaps:LMQ"),
        ("metanetx:MNXM10", "chemical", None, False, "metanetx:MNXM10"),
        ("metanetx:MNXM20", "chemical", None, False, "metanetx:MNXM20"),
        ("metanetx:MNXM21", "chemical", None, False, "metanetx:MNXM21"),
        ("metanetx:MNXM4", "chemical", None, False, "metanetx:MNXM4"),
        ("metanetx_reaction:MNXR2", "reaction", None, False, "metanetx_reaction:MNXR2"),
        ("metanetx_reaction:MNXR3", "reaction", None, False, "metanetx_reaction:MNXR3"),
        ("metanetx_reaction:MNXR4", "reaction", None, False, "metanetx_reaction:MNXR4"),
        ("pubchem:101", "chemical", None, False, "pubchem:101"),
        ("ramp:RAMP_C_1", "chemical", None, False, "ramp:RAMP_C_1"),
        ("ramp_gene:RAMP_G_3", "gene", None, False, "ramp_gene:RAMP_G_3"),
        ("ramp_gene:RAMP_G_4", "gene", None, False, "ramp_gene:RAMP_G_4"),
    ]


EXPECTED_RAMP_GENES = [
    ("ramp_gene:RAMP_G_1", "entrez:7157", "explicit_source_gene"),  # via its UniProt source id
    ("ramp_gene:RAMP_G_2", "entrez:7157", "explicit_source_gene"),  # via its NCBI Gene source id
    ("ramp_gene:RAMP_G_3", "ramp_gene:RAMP_G_3", "source_gene_only"),
    ("ramp_gene:RAMP_G_4", "ramp_gene:RAMP_G_4", "ambiguous_native"),  # P99999 links two genes
    ("ramp_gene:RAMP_G_5", "entrez:7157", "explicit_source_gene"),  # uniprot and entrez agree
]


def test_ramp_gene_source_gene_mappings(q):
    # Exactly one row per ramp_gene record.
    assert (
        q(
            "SELECT record_id,entity_id,decision FROM $exceptions WHERE starts_with(record_id,'ramp_gene:') ORDER BY 1"
        )
        == EXPECTED_RAMP_GENES
    )
    assert q(
        "SELECT count(*) FROM (SELECT record_id FROM $exceptions GROUP BY 1 HAVING count(*)>1)"
    ) == [(0,)]


def test_gene_products(q):
    rows = q("SELECT protein_entity_id,entrez_id,taxon FROM $products_by_protein")
    assert rows == sorted(rows) and len(rows) == 20
    assert ("uniprot:P04637", "7157", "9606") in rows
    # Neither gene's symbol is the gene name CYCS: both stay.
    assert ("uniprot:P99999", "7158", "9606") in rows and ("uniprot:P99999", "7159", "9606") in rows
    assert [r[1] for r in rows if r[0] == "uniprot:P0C0L4"] == ["720"]
    assert [r[1] for r in rows if r[0] == "uniprot:P0DN86"] == ["1082", "93659", "94115"]
    assert ("uniprot:Q00003", "5555", "9606") in rows  # target record absent: kept
    assert not [r for r in rows if r[0] == "uniprot:Q00002"]  # known, different organisms: dropped
    by_gene = q("SELECT entrez_id,protein_entity_id,taxon FROM $products_by_gene")
    assert by_gene == sorted(by_gene) and len(by_gene) == 20


def test_lipid_structures_rule5(q):
    # Both records with the name agree on one InChIKey. FA 20:0;5Me is carried by two different
    # keys and the species-level PC 34:1 never resolves to a structure.
    assert q("SELECT * FROM $lipid_structures") == [("full_structure:FA 18:0;5Me", KE)]


def test_manifest_and_resume(built):
    root, manifests, snapshot = built
    assert snapshot["format"] == "omnipath-identity-v2"
    assert set(snapshot["hub_indexes"]) == set(manifests)
    assert snapshot["counts"]["decisions"]["groups_accepted"] == 1
    assert snapshot["counts"]["decisions"]["groups_rejected"] == 1
    assert {"gene_products", "structures", "decisions", "total"} <= set(snapshot["timings"])
    again = build_identity(root / "hubindex", root / "identity", **ARGS)
    assert again["fingerprint"] == snapshot["fingerprint"]
    assert (root / "identity" / snapshot["fingerprint"] / "manifest.json").is_file()
    assert json.loads((root / "identity" / snapshot["fingerprint"] / "manifest.json").read_text())[
        "rules_sha256"
    ]


def test_record_candidates_for_unmerged_records(q):
    # Records that point to several anchors are not merged but stay evidence for each anchor.
    rows = q(
        "SELECT record_id,entity_id FROM read_parquet('"
        + str(q.out / "record_candidates.parquet")
        + "') ORDER BY 1,2"
    )
    assert rows == sorted(rows)
    by_record = {}
    for record_id, entity_id in rows:
        by_record.setdefault(record_id, []).append(entity_id)
    # A record claiming two InChIKeys: both keys.
    assert by_record["hmdb:HMDB0000002"] == ["inchikey:" + KB, "inchikey:" + KC]
    # A lipid record claiming two names at its most specific level: both names.
    assert len(by_record["lipidmaps:LMQ"]) == 2 and all(
        e.startswith("goslin:") for e in by_record["lipidmaps:LMQ"]
    )
    # A record cross-referencing two anchored records: both anchors.
    assert by_record["bigg:b1"] == ["inchikey:" + KA, "inchikey:" + KD]
    # A record reaching a quarantined record: that record's keys.
    assert by_record["bigg:b2"] == ["inchikey:" + KB, "inchikey:" + KC]
    # A ramp_gene record whose protein links two genes: both genes.
    assert by_record["ramp_gene:RAMP_G_4"] == ["entrez:7158", "entrez:7159"]
    # Rejected groups touch no anchor: no candidates.
    assert "metanetx:MNXM20" not in by_record
    # Every goslin candidate is an entity.
    extra = {r[0] for r in q("SELECT entity_id FROM $entities_extra")}
    lipid = {r[0] for r in q("SELECT entity_id FROM $exceptions WHERE decision='lipid_name'")}
    assert {e for e in by_record["lipidmaps:LMQ"]} <= extra | lipid
