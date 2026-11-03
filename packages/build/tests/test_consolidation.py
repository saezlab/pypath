"""Unit tests for consolidation, deduplication, and nested structure generation."""

import tempfile
import unittest
from pathlib import Path

from writer_fixture import write_observations
from omnipath_resolver import EntityResolver
from omnipath_build.silver import RawEntityObservation, RawRelationObservation

from library_fixture import build_fixture_library


class TestConsolidation(unittest.TestCase):
    def test_consolidation_merges_duplicate_relations(self):
        resolver = EntityResolver(library_dir=None)

        raw_entities = {
            "e1": RawEntityObservation(
                entity_key="e1",
                entity_type="protein",
                namespace="uniprot",
                identifier="P04637",
                taxon="9606",
                identifiers=[
                    {"ns": "uniprot", "id": "P04637", "is_canonical": True, "source": "s1"}
                ],
            ),
            "e2": RawEntityObservation(
                entity_key="e2",
                entity_type="protein",
                namespace="uniprot",
                identifier="P38398",
                taxon="9606",
                identifiers=[
                    {"ns": "uniprot", "id": "P38398", "is_canonical": True, "source": "s1"}
                ],
            ),
        }

        raw_relations = [
            RawRelationObservation(
                relation_key="r1",
                subject_entity_key="e1",
                predicate="affects",
                object_entity_key="e2",
                source="signor",
                dataset="complex",
                row_id="1",
                upstream_id="up1",
            ),
            RawRelationObservation(
                relation_key="r2",
                subject_entity_key="e1",
                predicate="affects",
                object_entity_key="e2",
                source="reactome",
                dataset="pathways",
                row_id="2",
                upstream_id="up2",
            ),
        ]

        raw_payloads = [
            {"relation_key": "r1", "source": "signor", "row_id": "1", "payload_json": '{"id": 1}'},
            {
                "relation_key": "r2",
                "source": "reactome",
                "row_id": "2",
                "payload_json": '{"id": 2}',
            },
        ]

        tables = write_observations(resolver, raw_entities, raw_relations, raw_payloads)
        entities, relations, payloads = (
            tables["entity"],
            tables["relation"],
            tables["evidence_payloads"],
        )

        resolver.close()
        self.assertEqual(len(entities), 2)
        self.assertEqual(
            len(relations), 1, "Duplicate edge should be merged into single canonical relation"
        )
        rel = relations[0]
        self.assertEqual(rel["evidence_count"], 2)
        self.assertIn("signor", rel["sources"])
        self.assertIn("reactome", rel["sources"])
        self.assertEqual(len(payloads), 2)

    def test_resolver_aliases_are_saved_on_entities(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            library = build_fixture_library(Path(tmpdir))
            resolver = EntityResolver(library_dir=library)
            tables = write_observations(
                resolver,
                {
                    "e1": RawEntityObservation(
                        entity_key="e1",
                        entity_type="protein",
                        namespace="uniprot",
                        identifier="P04637",
                        taxon="9606",
                        identifiers=[
                            {"ns": "entrez", "id": "7157", "is_canonical": True, "source": "s1"}
                        ],
                    )
                },
                [],
                [],
            )
            resolver.close()
            entities = tables["entity"]
            self.assertEqual(len(entities), 2)
            ent = next(e for e in entities if e["namespace"] == "entrez")
            product = next(e for e in entities if e["namespace"] == "uniprot")
            self.assertEqual(ent["identifier"], "7157")
            self.assertEqual(ent["entity_type"], "protein")
            self.assertEqual(product["identifier"], "P04637")
            self.assertEqual(ent["reference_entity_key"], "entrez:7157")
            self.assertEqual(product["reference_entity_key"], ent["reference_entity_key"])
            self.assertNotEqual(product["entity_key"], ent["entity_key"])
            self.assertEqual(ent["label"], "TP53")
            self.assertEqual(ent["taxon"], "9606")
            identifiers = tables.children("entity_identifier", ent)
            pairs = {(item["ns"], item["id"]) for item in identifiers}
            self.assertIn(("entrez", "7157"), pairs)
            self.assertIn(("uniprot", "P04637"), pairs)
            self.assertIn(("genesymbol", "TP53"), pairs)
            self.assertIn(("ensg", "ENSG00000141510"), pairs)
            self.assertTrue(
                any(
                    item["ns"] == "genesymbol" and item["source"] == "resolver"
                    for item in identifiers
                )
            )


if __name__ == "__main__":
    unittest.main()
