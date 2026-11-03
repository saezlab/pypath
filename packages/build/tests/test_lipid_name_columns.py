from unittest.mock import patch

import pyarrow.parquet as pq

from omnipath_build.hubs.sources import lipidmaps, swisslipids
from omnipath_build.hubs.writer import HubParquetWriter


def test_swisslipids_abbreviation_is_preserved(tmp_path):
    path = tmp_path / "swiss.parquet"
    writer = HubParquetWriter(path)
    lines = [
        "Lipid ID\tName\tAbbreviation*\n",
        "SLM:1\tPhosphatidylinositol (16:0/18:1)\tPI(16:0/18:1)\n",
    ]
    with patch.object(swisslipids, "iter_url_lines", return_value=iter(lines)):
        swisslipids.emit(writer)
    writer.close()
    values = {(r["source_type"], r["source_id"]) for r in pq.read_table(path).to_pylist()}
    assert ("lipid_shorthand", "PI(16:0/18:1)") in values
    assert ("name", "Phosphatidylinositol (16:0/18:1)") in values


def test_lipidmaps_current_name_and_goslin_columns(tmp_path):
    path = tmp_path / "lipidmaps.parquet"
    writer = HubParquetWriter(path)
    row = {
        "LM_ID": "LMFA01020216",
        "NAME": "fatty acid label",
        "ABBREVIATION": "FA 18:0;5Me",
        "SYSTEMATIC_NAME": "5-methyl-octadecanoic acid",
    }
    lines = [line for k, v in row.items() for line in [f"> <{k}>\n", v + "\n", "\n"]] + ["$$$$\n"]
    with patch.object(lipidmaps, "iter_url_lines", return_value=iter(lines)):
        lipidmaps.emit(writer)
    writer.close()
    values = {(r["source_type"], r["source_id"]) for r in pq.read_table(path).to_pylist()}
    assert ("name", "fatty acid label") in values
    assert ("lipid_shorthand", "FA 18:0;5Me") in values
    assert ("systematic_name", "5-methyl-octadecanoic acid") in values
