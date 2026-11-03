"""End-to-end integration test for the full-parquet build pipeline."""

import tempfile
import unittest
from pathlib import Path
import pyarrow.parquet as pq

from omnipath_build.silver import SilverExtractor
from omnipath_resolver import EntityResolver
from omnipath_build.writer import ParquetWriter


class TestBuildPipelineEndToEnd(unittest.TestCase):
    def test_pipeline_transformation_flow(self):
        # 1. Simulate silver extraction from raw items
        extractor = SilverExtractor(source="signor", dataset="test_ds")
        raw_items = [
            {
                "subject": {
                    "type": "protein",
                    "identifiers": [
                        {"type": "uniprot", "value": "P04637"},
                        {"type": "genesymbol", "value": "TP53"},
                    ],
                    "annotations": [{"term": "description", "value": "nucleus"}],
                },
                "predicate": "affects",
                "object": {
                    "type": "protein",
                    "identifiers": [
                        {"type": "uniprot", "value": "P38398"},
                        {"type": "genesymbol", "value": "BRCA1"},
                    ],
                    "annotations": [{"term": "description", "value": "nucleus"}],
                },
                "annotations": [
                    {"term": "description", "value": "phosphorylation"},
                    {"term": "object_direction_qualifier", "value": "increased"},
                ],
            }
        ]

        for idx, item in enumerate(raw_items):
            extractor.process_record(item, raw_payload=item, row_id=str(idx), row_number=idx)

        self.assertEqual(len(extractor.entities), 2)
        self.assertEqual(len(extractor.relations), 1)

        # 2. Local entity resolution
        resolver = EntityResolver(library_dir=None)
        with tempfile.TemporaryDirectory() as tmpdir:
            target = Path(tmpdir) / "resources" / "signor" / "vtest"
            writer = ParquetWriter(target)
            writer.append_observations(extractor, resolver)
            resolver.close()
            result = writer.close()
            (relation,) = pq.read_table(result["files"]["relation"]).to_pylist()
            self.assertEqual(relation["predicate"], "affects")
            self.assertEqual(relation["category"], "interaction")
            self.assertEqual(relation["sign"], 1)
            self.assertTrue(all(path.exists() for path in result["files"].values()))
            self.assertEqual(result["rows"]["entity"], 2)
            self.assertEqual(result["rows"]["relation"], 1)
            self.assertEqual(result["rows"]["evidence_payloads"], 1)


if __name__ == "__main__":
    unittest.main()


def test_record_limit_reaches_parser_without_reading_an_extra_row(
    tmp_path, monkeypatch, mapper_module
):
    from omnipath_build.discovery import DiscoveredDataset
    from omnipath_build import pipeline

    class RawDataset:
        def raw(self, **kwargs):
            assert kwargs["max_records"] == 10
            for index in range(10):
                yield {"id": str(index)}
            raise AssertionError("The pipeline requested an eleventh record")

    dataset = DiscoveredDataset(
        source="fixture",
        dataset_name="items",
        qualified_module=mapper_module.__name__,
        call=lambda: (),
        raw_dataset=RawDataset(),
        mapper=lambda row: {
            "type": "named_thing",
            "identifiers": [{"type": "fixture", "value": row["id"]}],
        },
    )
    monkeypatch.setattr(pipeline, "discover_datasets", lambda **kwargs: ("fixture", [dataset], {}))
    pipeline.build_resource(
        "fixture", batch_workers=1, version="1", output_dir=tmp_path, max_records=10, progress=False
    )
    assert pq.read_metadata(tmp_path / "resources/fixture/1/entity.parquet").num_rows == 10
