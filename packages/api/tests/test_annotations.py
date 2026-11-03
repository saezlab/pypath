"""Unit tests for presentation annotation and CURIE transformations."""

import unittest
from omnipath_api.queries.annotations import (
    normalize_identifier_type,
    split_compiled_annotations,
    unwrap_attribute_value,
)


class TestAnnotations(unittest.TestCase):
    def test_unwrap_attribute_value(self):
        self.assertEqual(unwrap_attribute_value('{"value": "123", "canonical": true}'), "123")
        self.assertEqual(unwrap_attribute_value("plain_val"), "plain_val")
        self.assertEqual(unwrap_attribute_value(None), "")

    def test_normalize_identifier_type(self):
        self.assertEqual(normalize_identifier_type("uniprot"), "uniprot")
        self.assertEqual(normalize_identifier_type("entrez"), "entrez")
        self.assertEqual(normalize_identifier_type("genesymbol"), "genesymbol")
        self.assertEqual(normalize_identifier_type("chebi"), "chebi")

    def test_compiled_annotations_preserve_semantics(self):
        annotations = [
            {"term": "object_direction_qualifier", "value": "increased", "scope": "relation"},
            {"term": "description", "value": "NA", "scope": "subject"},
            {"term": "original_predicate", "value": "MI:2236", "scope": "relation"},
        ]
        result = split_compiled_annotations(annotations)
        self.assertEqual(result["subject"], [{"term": "description", "value": "NA"}])
        self.assertEqual(result["relation"][1]["value"], "MI:2236")
        self.assertEqual(normalize_identifier_type("unknown", "7157"), "unknown")


if __name__ == "__main__":
    unittest.main()
