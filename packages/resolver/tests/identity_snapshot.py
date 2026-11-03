"""A small synthetic ``omnipath-identity-v2`` layout written straight from spec 3a and 3b.

Per-hub indexes (``records``, ``by_id``, ``by_record``) plus the identity decisions
(``exceptions``, ``entities_extra``, ``gene_products_*``, ``lipid_structures``). The index and
decision builders (``build_hub_index``, ``build_identity``) are not involved: the helpers below
derive ``by_id`` rows from hub rows with the expansion rules of spec 3a, so each fixture record
is one case for a rule. The runtime reads LMDB only, so the real kv writers
(``build_hub_kv_dir``, ``build_identity_kv``) turn those Parquet files into the kv stores.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

import lmdb
import pyarrow as pa
import pyarrow.parquet as pq

from omnipath_resolver import identity_kv

WATER = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
ASPIRIN = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
ETHANOL = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
LIPID = "AAAAAAAAAAAAAA-LLLLLLLLLL-N"
LIPID2 = "AAAAAAAAAAAAAA-MMMMMMMMMM-N"
LIPID3 = "AAAAAAAAAAAAAA-NNNNNNNNNN-N"
FATTY = "AAAAAAAAAAAAAA-FFFFFFFFFF-N"
AMINO = "AAAAAAAAAAAAAA-GGGGGGGGGG-N"
ANION = "AAAAAAAAAAAAAA-HHHHHHHHHH-N"
NAMELESS = "BBBBBBBBBBBBBB-CCCCCCCCCC-N"
ONLY_PC = "CCCCCCCCCCCCCC-DDDDDDDDDD-N"
OTHERS = "DDDDDDDDDDDDDD-EEEEEEEEEE-N"
SYSTEMATIC = "EEEEEEEEEEEEEE-FFFFFFFFFF-N"
TWO_A, TWO_B = "FFFFFFFFFFFFFF-GGGGGGGGGG-N", "FFFFFFFFFFFFFF-HHHHHHHHHH-N"
TWO_C, TWO_D = "IIIIIIIIIIIIII-JJJJJJJJJJ-N", "IIIIIIIIIIIIII-KKKKKKKKKK-N"
NADH, NADH_2 = "PPPPPPPPPPPPPP-QQQQQQQQQQ-N", "PPPPPPPPPPPPPP-QQQQQQQQQQ-L"  # protonation states
NADH_1 = "PPPPPPPPPPPPPP-QQQQQQQQQQ-M"
FINGERPRINT = "synthetic-v2-0001"
LONG_NAME = "x" * 81

CHEMICAL_HUBS = ("chebi", "hmdb", "chembl", "pubchem", "lipidmaps", "swisslipids", "bigg", "kegg")
GENE_HUBS = ("uniprot", "entrez", "ramp_gene")
# source types that are never lookup rows
NOT_KEYS = {"name", "synonym", "lipid_shorthand", "systematic_name", "smiles", "inchi", "formula"}


def part(value: str) -> str:
    return hashlib.md5(value.encode(), usedforsecurity=False).hexdigest()[:2]


def hub_records():
    """hub -> {local_id: (taxon, reviewed, [(source_type, value)])}."""
    h = {hub: {} for hub in CHEMICAL_HUBS + GENE_HUBS}

    def rec(hub, local, rows, taxon=None, reviewed=False):
        h[hub][local] = (taxon, reviewed, rows)

    # ---- chemicals -----------------------------------------------------------------
    rec("chebi", "CHEBI:15377", [
        ("chebi", "CHEBI:15377"), ("name", "water"), ("name", "Water molecule"),
        ("synonym", "H2O"), ("inchikey", WATER), ("smiles", "O"), ("kegg", "C00001"),
    ])  # fmt: skip
    rec(
        "hmdb", "HMDB0002111", [("name", "Water (hmdb)"), ("inchikey", WATER), ("cas", "7732-18-5")]
    )
    rec("pubchem", "962", [("name", "PUBCHEM WATER"), ("inchikey", WATER)])
    rec("pubchem", "7777", [("name", "attached to water")])  # no key: attached by decision
    rec("chembl", "CHEMBL999", [("inchikey", WATER), ("name", "overridden")])  # exception wins
    rec("chebi", "CHEBI:15365", [("name", ASPIRIN), ("inchikey", ASPIRIN)])  # name is a key
    rec("hmdb", "HMDB0001879", [("name", LONG_NAME), ("inchikey", ASPIRIN)])  # over 80 chars
    rec("chembl", "CHEMBL25", [("name", "ASPIRIN"), ("name", "aspirin2"), ("inchikey", ASPIRIN)])
    rec("pubchem", "2244", [("name", "aspirin-pc"), ("inchikey", ASPIRIN)])
    rec("chebi", "CHEBI:16236", [("name", "ethanol"), ("inchikey", ETHANOL)])
    rec("chebi", "CHEBI:50000", [
        ("name", "phosphatidylcholine"), ("goslin", "species:PC 34:1"),
        ("goslin", "full_structure:PC 16:0/18:1"), ("inchikey", LIPID),
    ])  # fmt: skip
    rec(
        "lipidmaps",
        "LM9",
        [
            ("goslin", "species:PC 34:1"),
            ("goslin", "full_structure:PC 16:0/18:1"),
            ("inchikey", LIPID2),
        ],
    )  # fmt: skip  (disagrees with CHEBI:50000 on the full-structure name)
    rec("chebi", "CHEBI:50001", [("goslin", "full_structure:PE 18:0/18:1"), ("inchikey", LIPID3)])
    # Goslin reads a systematically named acid as a fatty acyl: a real fatty acid,
    # which a lipid database also holds, and an amino acid, which none does.
    rec("chebi", "CHEBI:16196", [
        ("name", "oleic acid"), ("goslin", "full_structure:FA 18:1(9Z)"), ("inchikey", FATTY),
    ])  # fmt: skip
    rec("lipidmaps", "LMFA01030002", [("goslin", "full_structure:FA 18:1(9Z)"), ("inchikey", FATTY)])
    rec("chebi", "CHEBI:18050", [
        ("name", "L-glutamine"), ("goslin", "full_structure:FA 5:0;2NH2,5NH2;5oxo"), ("inchikey", AMINO),
    ])  # fmt: skip
    rec("chebi", "CHEBI:17120", [
        ("name", "hexanoate"), ("goslin", "full_structure:WE 6:0/0:0"), ("inchikey", ANION),
    ])  # fmt: skip
    rec("lipidmaps", "LMFA01010006", [
        ("goslin", "full_structure:WE 6:0/0:0"), ("goslin", "full_structure:FA 6:0"), ("inchikey", ANION),
    ])  # fmt: skip
    rec("swisslipids", "SLM:1", [("name", "a long lipid name"), ("goslin", "species:PE 36:2")])
    rec("hmdb", "HMDB0000001", [("inchikey", NAMELESS)])
    rec("pubchem", "5", [("name", "pc-name"), ("inchikey", ONLY_PC)])
    rec("lipidmaps", "LM1", [("name", "lm"), ("inchikey", ONLY_PC)])
    rec(
        "swisslipids", "SLM:2", [("name", "zz"), ("name", "b"), ("name", "a"), ("inchikey", OTHERS)]
    )
    rec("lipidmaps", "LM2", [("systematic_name", "sys-name"), ("inchikey", SYSTEMATIC)])
    rec("chebi", "CHEBI:99999", [("name", "multi"), ("inchikey", TWO_A), ("inchikey", TWO_B)])
    rec("pubchem", "P2", [("name", "points to two anchors"), ("inchikey", TWO_C), ("inchikey", TWO_D)])
    rec("hmdb", "HMDB0088888", [("name", "second anchor holder"), ("inchikey", TWO_C)])
    rec("chebi", "CHEBI:80001", [("name", "NADH"), ("inchikey", NADH)])
    rec("chebi", "CHEBI:80002", [("name", "NADH(2-)"), ("inchikey", NADH_2)])
    rec("chebi", "CHEBI:80003", [("name", "NADH(1-)"), ("inchikey", NADH_1)])
    rec("chebi", "CHEBI:70001", [("name", "grouped-chebi")])
    rec("hmdb", "HMDB0070001", [("name", "grouped-hmdb")])
    for n in range(12):  # twelve structure-less records behind one cas number
        rec("pubchem", str(1000 + n), [("name", f"cas-twin-{n}"), ("cas", "0-0-0")])
    rec("bigg", "x1", [("inchikey", ETHANOL)])
    rec("bigg", "y1", [("bigg", "x1"), ("inchikey", ASPIRIN)])  # a claim: fallback behind x1

    # ---- genes and proteins --------------------------------------------------------------
    def gene(n, taxon, rows):
        rec("entrez", str(n), [("entrez", str(n)), *rows], taxon)

    gene(7157, "9606", [
        ("genesymbol", "TP53"), ("genesymbol-syn", "P53"), ("ensg", "ENSG00000141510"),
        ("refseq", "NM_000546"), ("refseq_protein", "NP_000537.3"), ("hgnc", "HGNC:11998"),
    ])  # fmt: skip
    gene(801, "9606", [("genesymbol", "CALM1")])
    gene(805, "9606", [("genesymbol", "CALM2")])
    gene(808, "9606", [("genesymbol", "CALM3")])
    gene(22059, "10090", [("genesymbol", "Trp53"), ("genesymbol-syn", "TP53")])
    gene(555, None, [("ensg", "ENSG00000999999")])  # no taxon: visible only unscoped
    gene(901, "9606", [("genesymbol", "ABC1")])
    gene(4001, "9606", [("genesymbol", "LINKED")])
    gene(902, "9606", [("genesymbol-syn", "ABC1"), ("genesymbol-syn", "XYZ")])
    for n in range(12):  # one ENSG id shared by twelve genes: no cap, so no resolution
        gene(3000 + n, "9606", [("ensg", "ENSG00000000001")])

    def protein(acc, taxon, rows, reviewed=False):
        rec("uniprot", acc, [("uniprot", acc), *rows], taxon, reviewed)

    protein("P04637", "9606", [
        ("uniprot_entry", "P53_HUMAN"), ("genesymbol", "TP53"), ("genesymbol-syn", "P53"),
        ("entrez", "7157"), ("ensg", "ENSG00000141510"), ("hgnc", "HGNC:11998"),
        ("uniprot-sec", "Q15086"), ("uniprot-sec", "Q99999"), ("uniprot-sec", "Q88888"),
        ("refseq_protein", "NP_000537.3"),
    ], True)  # fmt: skip
    protein("P04637-2", "9606", [("ensp", "ENSP00000999999")])  # isoform, parent exists
    protein("Q11111-2", "9606", [("ensp", "ENSP00000888888")])  # isoform, parent absent
    protein("A0A0U1RQF1", "9606", [
        ("genesymbol-syn", "TPX"), ("uniprot_entry", "A0A0U1RQF1_HUMAN"), ("entrez", "7157"),
    ])  # fmt: skip
    protein("Q99999", "9606", [])
    protein("P0DP23", "9606", [
        ("genesymbol", "CALM1"), ("hgnc", "HGNC:1442"), ("uniprot-sec", "Q88888"),
        ("entrez", "801"), ("entrez", "805"), ("entrez", "808"),
    ], True)  # fmt: skip
    protein("P02340", "10090", [("genesymbol", "Trp53"), ("entrez", "22059")], True)
    protein("Q77777", "9606", [("entrez", "999"), ("genesymbol", "NOREC")])  # gene without record
    # a protein record without a taxon: its symbol rows belong to every scope
    protein("Q55555", None, [("genesymbol", "LINKED"), ("entrez", "4001")])
    rec("ramp_gene", "RAMP_G_1", [("ramp_gene", "RAMP_G_1")])
    rec("ramp_gene", "RAMP_G_2", [("ramp_gene", "RAMP_G_2")], "9606")
    return h


EXCEPTIONS = [  # record_id, entity_id, decision, quarantined
    ("swisslipids:SLM:1", "goslin:species:PE 36:2", "lipid_name", False),
    ("chebi:CHEBI:99999", "chebi:CHEBI:99999", "quarantined", True),
    ("pubchem:P2", "pubchem:P2", "quarantined", True),
    ("pubchem:7777", f"inchikey:{WATER}", "attached", False),
    ("chebi:CHEBI:70001", "chebi:CHEBI:70001", "grouped", False),
    ("hmdb:HMDB0070001", "chebi:CHEBI:70001", "grouped", False),
    ("chembl:CHEMBL999", "chembl:CHEMBL999", "ambiguous_native", False),
    ("ramp_gene:RAMP_G_1", "entrez:7157", "explicit_source_gene", False),
    ("ramp_gene:RAMP_G_2", "ramp_gene:RAMP_G_2", "source_gene_only", False),
]
EXTRAS = [  # entity_id, kind, taxon, quarantined, preferred_record
    ("goslin:species:PE 36:2", "chemical", None, False, "swisslipids:SLM:1"),
    ("chebi:CHEBI:99999", "chemical", None, True, "chebi:CHEBI:99999"),
    ("pubchem:P2", "chemical", None, True, "pubchem:P2"),
    ("chebi:CHEBI:70001", "chemical", None, False, "chebi:CHEBI:70001"),
    ("chembl:CHEMBL999", "chemical", None, False, "chembl:CHEMBL999"),
    ("ramp_gene:RAMP_G_2", "gene", "9606", False, "ramp_gene:RAMP_G_2"),
]
GENE_PRODUCTS = [  # protein_entity_id, entrez_id, taxon
    ("uniprot:P04637", "7157", "9606"),
    ("uniprot:A0A0U1RQF1", "7157", "9606"),
    ("uniprot:P0DP23", "801", "9606"),
    ("uniprot:P0DP23", "805", "9606"),
    ("uniprot:P0DP23", "808", "9606"),
    ("uniprot:P02340", "22059", "10090"),
    ("uniprot:Q77777", "999", "9606"),
]
CANDIDATES = [  # record_id, entity_id: a quarantined record points to the anchors it claims
    ("pubchem:P2", f"inchikey:{TWO_C}"),
    ("pubchem:P2", f"inchikey:{TWO_D}"),
]
LIPID_STRUCTURES = [("full_structure:PE 18:0/18:1", f"inchikey:{LIPID3}")]  # rule 5: unique names


def anchor_of(hub, rows, local):
    """(anchor, anchor_count) of a record, as the hub index would compute them."""
    if hub in CHEMICAL_HUBS:
        keys = sorted({v for st, v in rows if st == "inchikey"})
        return (f"inchikey:{keys[0]}" if len(keys) == 1 else None), len(keys)
    if hub == "uniprot":
        return f"uniprot:{local}", 1
    return None, 0


def by_id_rows(hub, local, rows):
    """Lookup rows of one record with the data-level expansions of spec 3a."""
    out = [(hub, local, "native")]
    for st, value in rows:
        if st in NOT_KEYS or (st == hub and value == local):
            continue
        out.append((st, value, "claim"))
        if st == "uniprot-sec":
            out.append(("uniprot", value, "secondary"))
        if st == "genesymbol-syn":
            out.append(("genesymbol", value, "symbol_synonym"))
        if st in ("refseq_protein", "genbank") and "." in value:
            out.append((st, value.rsplit(".", 1)[0], "version_stripped"))
    return out


def _write(path: Path, table: pa.Table) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path, row_group_size=64 * 1024)


S = pa.string()
RECORDS = pa.schema(
    [("local_id", S), ("record_id", S), ("taxon", S), ("anchor", S), ("anchor_count", pa.int64()),
     ("reviewed", pa.bool_())]
)  # fmt: skip
BY_ID = pa.schema(
    [("ns", S), ("identifier", S), ("local_id", S), ("tag", S), ("taxon", S), ("anchor", S),
     ("anchor_count", pa.int64())]
)  # fmt: skip
BY_RECORD = pa.schema([("local_id", S), ("source_type", S), ("value", S)])


def write_hub_index(directory: Path, hub: str, records) -> None:
    rows, by_id, by_record = [], defaultdict(list), defaultdict(list)
    for local in sorted(records):
        taxon, reviewed, items = records[local]
        anchor, count = anchor_of(hub, items, local)
        rows.append(dict(local_id=local, record_id=f"{hub}:{local}", taxon=taxon, anchor=anchor,
                         anchor_count=count, reviewed=reviewed))  # fmt: skip
        for ns, identifier, tag in by_id_rows(hub, local, items):
            by_id[part(identifier)].append(dict(ns=ns, identifier=identifier, local_id=local,
                tag=tag, taxon=taxon, anchor=anchor, anchor_count=count))  # fmt: skip
        for st, value in items:
            by_record[part(f"{hub}:{local}")].append(
                dict(local_id=local, source_type=st, value=value)
            )
    _write(directory / "records.parquet", pa.Table.from_pylist(rows, schema=RECORDS))
    for p, items in by_id.items():
        items.sort(key=lambda r: (r["ns"], r["identifier"], r["local_id"]))
        _write(
            directory / f"by_id/part={p}/data.parquet", pa.Table.from_pylist(items, schema=BY_ID)
        )
    for p, items in by_record.items():
        items.sort(key=lambda r: r["local_id"])
        _write(directory / f"by_record/part={p}/data.parquet",
               pa.Table.from_pylist(items, schema=BY_RECORD))  # fmt: skip
    (directory / "manifest.json").write_text(
        json.dumps(dict(hub=hub, counts=dict(records=len(rows))))
    )


def downgrade_to_v2(directory: Path) -> None:
    """Rewrite the id store of a built hub kv in format v2: no taxon in the symbol keys."""
    kv = Path(directory) / identity_kv.KV_DIR
    manifest = json.loads((kv / identity_kv.MANIFEST).read_text())
    rows = defaultdict(set)
    if (Path(directory) / "by_id").is_dir():  # a hub without records has none
        for row in pq.read_table(Path(directory) / "by_id").to_pylist():
            rows[row["ns"], row["identifier"]].add((row["local_id"], row["tag"]))
    shards = defaultdict(dict)
    for (ns, identifier), items in rows.items():
        key = identity_kv.id_key(ns, identifier)
        if key is not None:
            value = identity_kv.ROW.join(identity_kv.FIELD.join(i) for i in sorted(items))
            shards[identity_kv.shard_of(identifier, manifest["id_shards"])][key] = value.encode()
    for shard in range(manifest["id_shards"]):
        shutil.rmtree(kv / f"id.{shard}")
        env = lmdb.open(str(kv / f"id.{shard}"), map_size=1 << 30, max_dbs=2)
        db = env.open_db(b"id")
        with env.begin(write=True, db=db) as txn:
            for key, value in sorted(shards[shard].items()):
                txn.put(key, value, db=db)
        env.close()
    manifest["format"] = "omnipath-hub-kv-v2"
    del manifest["taxon_keyed_namespaces"]
    (kv / identity_kv.MANIFEST).write_text(json.dumps(manifest))


def build_snapshot(
    root: Path, *, fingerprint: str = FINGERPRINT, shards: int = 1, v2_hubs=()
) -> Path:
    """Write hub indexes under ``root/hubindex`` and the identity directory; return the latter.

    ``shards`` (1 or 16) is the id/rec shard count of every hub kv store. The kv stores of
    ``v2_hubs`` ("all" for every hub) are written in the older format v2.
    """
    from omnipath_build.identity import build_hub_kv_dir, build_identity_kv

    root = Path(root)
    hubs = hub_records()
    identity = root / "identity" / fingerprint
    identity.mkdir(parents=True)
    indexes = {}
    for hub, records in hubs.items():
        directory = root / "hubindex" / hub / "0123456789ab"
        write_hub_index(directory, hub, records)
        build_hub_kv_dir(directory, id_shards=shards, rec_shards=shards, memory="256MB", threads=1, min_free_gib=0.01)
        if v2_hubs == "all" or hub in v2_hubs:
            downgrade_to_v2(directory)
        # relative to the identity directory, except uniprot (absolute): both must resolve
        indexes[hub] = (
            str(directory) if hub == "uniprot" else str(Path("../..") / directory.relative_to(root))
        )

    def table(name, schema, rows, key):
        rows = sorted(rows, key=lambda r: r[key])
        _write(identity / name, pa.Table.from_pylist(rows, schema=pa.schema(schema)))

    exc = [dict(record_id=r, entity_id=e, decision=d, quarantined=q) for r, e, d, q in EXCEPTIONS]
    cols = [("record_id", S), ("entity_id", S), ("decision", S), ("quarantined", pa.bool_())]
    table("exceptions.parquet", cols, exc, "record_id")
    table("exception_members.parquet", [cols[1], cols[0], cols[2], cols[3]], exc, "entity_id")
    table(
        "entities_extra.parquet",
        [
            ("entity_id", S),
            ("kind", S),
            ("taxon", S),
            ("quarantined", pa.bool_()),
            ("preferred_record", S),
        ],  # fmt: skip
        [
            dict(entity_id=e, kind=k, taxon=t, quarantined=q, preferred_record=p)
            for e, k, t, q, p in EXTRAS
        ],  # fmt: skip
        "entity_id",
    )
    gp = [dict(protein_entity_id=p, entrez_id=g, taxon=t) for p, g, t in GENE_PRODUCTS]
    gp_cols = [("protein_entity_id", S), ("entrez_id", S), ("taxon", S)]
    table("gene_products_by_protein.parquet", gp_cols, gp, "protein_entity_id")
    table("gene_products_by_gene.parquet", gp_cols, gp, "entrez_id")
    table(
        "record_candidates.parquet",
        [("record_id", S), ("entity_id", S)],
        [dict(record_id=r, entity_id=e) for r, e in CANDIDATES],
        "record_id",
    )
    table(
        "lipid_structures.parquet",
        [("goslin", S), ("inchikey", S)],
        [dict(goslin=g, inchikey=k) for g, k in LIPID_STRUCTURES],
        "goslin",
    )
    (identity / "manifest.json").write_text(
        json.dumps(
            dict(
                format="omnipath-identity-v2",
                fingerprint=fingerprint,
                hub_indexes=indexes,
                rules_sha256="0" * 64,
                counts={},
                timings={},
            )
        )
    )
    build_identity_kv(identity, min_free_gib=0.001)
    return identity
