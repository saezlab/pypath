"""The MetSigDB Parquet product over a small release, against the rows frozen in ``golden/``
(see its README).

The rows carry more than the frozen ones (exclusion classes, structure flag); the columns
below are compared one by one. MACdb set ids are prefixed in the frozen rows.
"""

import json
import os
from pathlib import Path

import duckdb

from omnipath_subsets.parquet_products import build_release_subsets
import metsigdb_fixture as fixture

GOLDEN = Path(__file__).parent / "golden" / "metsigdb_membership.json"


def parquet_rows(root, version):
    path = root / "subsets" / version / "metsigdb" / "metsigdb_membership.parquet"
    rows = (
        duckdb.connect()
        .execute(
            "SELECT resource, set_source_id, metabolite_entity_key, metabolite_label, "
            "metabolite_entity_type, inchikey, metabolite_structure_key, smiles, hmdb, pubchem, "
            "chebi, kegg, set_entity_key, set_label, set_type, set_sub_type, organism, set_size, "
            f"set_context, provenance_source, provenance_record FROM '{path}'"
        )
        .fetchall()
    )
    normalised = []
    for row in rows:
        row = list(row)
        record = json.loads(row[20])
        row[20] = json.dumps(
            {"row_id": str(record["row_id"]), "via_reaction": record.get("via_reaction")},
            sort_keys=True,
        )
        row[18] = json.dumps(json.loads(row[18]), sort_keys=True) if row[18] is not None else None
        normalised.append(tuple(row))
    return sorted(normalised)


def test_the_product_equals_the_frozen_rows(tmp_path):
    version = fixture.write_parquet_release(tmp_path)
    build_release_subsets(tmp_path, version)
    actual = parquet_rows(tmp_path, version)
    if os.environ.get("UPDATE_GOLDEN"):
        GOLDEN.write_text(json.dumps([list(row) for row in actual], indent=1) + "\n")
    expected = sorted(tuple(row) for row in json.loads(GOLDEN.read_text()))
    assert expected, "the fixture must give rows"
    assert {row[0] for row in actual} == {"Reactome", "WikiPathways", "KEGG", "MACdb", "ClassyFire"}
    assert actual == expected
