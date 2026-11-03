import pyarrow.parquet as pq

from omnipath_build.hubs.sources.refmet import emit
from omnipath_build.hubs.writer import HubParquetWriter


def test_refmet_preserves_precise_and_unspecified_records(tmp_path):
    path = tmp_path / "refmet.parquet"
    writer = HubParquetWriter(path)
    emit(
        writer,
        [
            {
                "refmet_id": "RM1",
                "refmet_name": "precise compound",
                "inchi_key": "AAAAAAAAAAAAAA-BBBBBBBBSA-N",
                "pubchem_cid": "12",
            },
            {"refmet_id": "RM2", "refmet_name": "PC 34:1", "inchi_key": "", "pubchem_cid": ""},
        ],
    )
    writer.close()
    rows = pq.read_table(path).to_pylist()
    fields = {(r["hub_id"], r["source_type"], r["source_id"]) for r in rows}
    assert ("RM1", "inchikey", "AAAAAAAAAAAAAA-BBBBBBBBSA-N") in fields
    assert ("RM1", "pubchem", "12") in fields
    assert ("RM2", "refmet", "RM2") in fields
    assert ("RM2", "name", "PC 34:1") in fields
    assert not any(r["hub_id"] == "RM2" and r["source_type"] == "inchikey" for r in rows)
