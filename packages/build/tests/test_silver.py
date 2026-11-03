"""Unit tests for silver extraction."""

import unittest
from omnipath_build.silver import SilverExtractor, entity_key


class TestSilverExtractor(unittest.TestCase):
    def test_entity_extraction(self):
        extractor = SilverExtractor(source="test_source", dataset="test_ds")
        record = {
            "type": "protein",
            "identifiers": [
                {"type": "uniprot", "value": "P04637"},
                {"type": "genesymbol", "value": "TP53"},
            ],
            "annotations": [
                {"term": "description", "value": "tumor suppressor"},
            ],
        }
        extractor.process_record(record, raw_payload={"raw": 1}, row_id="0", row_number=0)
        self.assertEqual(len(extractor.entities), 1)
        ent = list(extractor.entities.values())[0]
        self.assertEqual(ent.entity_type, "protein")
        self.assertEqual(ent.identifier, "P04637")
        self.assertEqual(len(ent.identifiers), 2)
        self.assertEqual(len(ent.annotations), 1)

    def test_relation_extraction(self):
        extractor = SilverExtractor(source="signor", dataset="test_ds")
        record = {
            "subject": {
                "type": "protein",
                "identifiers": [{"type": "uniprot", "value": "P04637"}],
            },
            "predicate": "affects",
            "object": {
                "type": "protein",
                "identifiers": [{"type": "uniprot", "value": "P38398"}],
            },
            "annotations": [
                {"term": "description", "value": "phosphorylation"},
            ],
        }
        extractor.process_record(record, raw_payload={"raw": "rel"}, row_id="42", row_number=42)
        self.assertEqual(len(extractor.entities), 2)
        self.assertEqual(len(extractor.relations), 1)
        self.assertEqual(len(extractor.payloads), 1)

        rel = extractor.relations[0]
        self.assertEqual(rel.predicate, "affects")
        self.assertEqual(rel.source, "signor")
        self.assertEqual(rel.row_id, "42")
        self.assertIn("phosphorylation", [a["value"] for a in rel.annotations])

    def test_annotation_scope_separation(self):
        extractor = SilverExtractor(source="signor", dataset="interactions")
        record = {
            "subject": {
                "type": "protein",
                "identifiers": [{"type": "uniprot", "value": "P04637"}],
                "annotations": [{"term": "description", "value": "Ser-15"}],
            },
            "predicate": "affects",
            "object": {
                "type": "protein",
                "identifiers": [{"type": "uniprot", "value": "P38398"}],
                "annotations": [{"term": "name", "value": "BRCT"}],
            },
            "annotations": [
                {"term": "publications", "value": "12345678"},
            ],
        }
        extractor.process_record(
            record, raw_payload={"raw": "payload"}, row_id="101", row_number=101
        )

        # Global entities should NOT contain participant annotations
        for ent in extractor.entities.values():
            self.assertEqual(
                len(ent.annotations),
                0,
                "Participant annotations must not leak into global entity table",
            )

        # Relation should contain all three scoped annotations
        rel = extractor.relations[0]
        scopes = {a["term"]: a["scope"] for a in rel.annotations}
        self.assertEqual(scopes.get("description"), "subject")
        self.assertEqual(scopes.get("name"), "object")
        self.assertEqual(scopes.get("publications"), "relation")

    def test_entity_associations_become_associated_with_edges(self):
        extractor = SilverExtractor(source="uniprot", dataset="proteins")
        record = {
            "type": "protein",
            "identifiers": [{"type": "uniprot", "value": "P04637"}],
            "associations": [
                {
                    "object": {
                        "type": "ontology_class",
                        "identifier_type": "cv_term_accession",
                        "identifier": "KW-0002",
                    },
                    "predicate": None,
                },
                {
                    "object": {
                        "type": "ontology_class",
                        "identifier_type": "cv_term_accession",
                        "identifier": "GO:0005515",
                    },
                },
            ],
        }
        extractor.process_record(record, raw_payload={"id": "P04637"}, row_id="0", row_number=0)

        self.assertEqual(len(extractor.entities), 3)
        assoc_rels = [r for r in extractor.relations if r.predicate == "associated_with"]
        self.assertEqual(len(assoc_rels), 2)

        protein = extractor.entities[entity_key("protein", "uniprot", "P04637")]
        keyword = extractor.entities[entity_key("ontology_class", "cv_term_accession", "KW-0002")]
        go_term = extractor.entities[
            entity_key("ontology_class", "cv_term_accession", "GO:0005515")
        ]
        self.assertEqual(protein.entity_type, "protein")
        self.assertEqual(keyword.namespace, "cv_term_accession")
        self.assertEqual(keyword.identifier, "KW-0002")
        self.assertEqual(
            {r.object_entity_key for r in assoc_rels},
            {keyword.entity_key, go_term.entity_key},
        )
        self.assertTrue(all(r.subject_entity_key == protein.entity_key for r in assoc_rels))

    def test_ontology_relation_entity_ref_keeps_cv_term_type(self):
        extractor = SilverExtractor(source="uniprot", dataset="keywords")
        record = {
            "type": "ontology_class",
            "identifiers": [{"type": "uniprot_keyword", "value": "KW-0002"}],
            "ontology_relations": [
                {
                    "predicate": "subclass_of",
                    "object": {
                        "type": "ontology_class",
                        "identifier_type": "uniprot_keyword",
                        "identifier": "KW-9993",
                    },
                    "ontology_id": "uniprot_keywords",
                }
            ],
        }
        extractor.process_record(record, raw_payload={}, row_id="1", row_number=1)
        self.assertEqual(len(extractor.relations), 1)
        rel = extractor.relations[0]
        self.assertEqual(rel.predicate, "subclass_of")
        parent = extractor.entities[entity_key("ontology_class", "uniprot_keyword", "KW-9993")]
        self.assertEqual(parent.entity_type, "ontology_class")
        self.assertEqual(rel.object_entity_key, parent.entity_key)

    def test_signor_curie_stays_signor_not_entrez(self):
        extractor = SilverExtractor(source="signor", dataset="complexes")

        from omnipath_core import Namespace

        extractor.process_record(
            {
                "type": "macromolecular_complex",
                "identifiers": [
                    {"type": Namespace.SIGNOR, "value": "SIGNOR-C1"},
                    {"type": "name", "value": "AP-1"},
                ],
            },
            raw_payload={},
            row_id="0",
            row_number=0,
        )
        extractor.process_record(
            {
                "type": "phenotypic_feature",
                "identifiers": [{"type": "signor", "value": "SIGNOR-PH1"}],
            },
            raw_payload={},
            row_id="1",
            row_number=1,
        )
        namespaces = {ent.namespace for ent in extractor.entities.values()}
        self.assertEqual(namespaces, {"signor"})
        self.assertTrue(
            all(ent.identifier.startswith("SIGNOR-") for ent in extractor.entities.values())
        )


if __name__ == "__main__":
    unittest.main()
