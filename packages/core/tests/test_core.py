"""Unit tests for omnipath_core components."""

import unittest
from biolink_model.datamodel.model import (
    slots,
    Protein,
    Gene,
    SmallMolecule,
    MacromolecularComplex,
)
from omnipath_core import (
    ENTITY_TABLE,
    RELATION_TABLE,
    PAYLOAD_SCHEMA,
    SILVER_ENTITY_SCHEMA,
    SILVER_RELATION_SCHEMA,
    Entity,
    Relation,
    Identifier,
    Annotation,
    Namespace,
    format_term,
    entity_key,
    relation_key,
    validate_version,
    validate_source,
    normalize_namespace,
)


class TestOmnipathCore(unittest.TestCase):
    def test_biolink_model_and_types(self):
        self.assertEqual(Protein.class_name, "protein")
        self.assertEqual(Gene.class_name, "gene")
        self.assertEqual(SmallMolecule.class_name, "small molecule")
        self.assertEqual(MacromolecularComplex.class_name, "macromolecular complex")
        self.assertEqual(slots.regulates.curie, "biolink:regulates")
        self.assertEqual(slots.physically_interacts_with.curie, "biolink:physically_interacts_with")

    def test_format_term(self):
        self.assertEqual(format_term(Protein), "protein")
        self.assertEqual(format_term(MacromolecularComplex), "macromolecular_complex")
        self.assertEqual(format_term(slots.regulates), "regulates")
        self.assertEqual(format_term(slots.physically_interacts_with), "physically_interacts_with")
        self.assertEqual(format_term(Namespace.UNIPROT), "uniprot")

    def test_silver_schema_models(self):
        ident = Identifier(type=Namespace.UNIPROT, value="P04637")
        self.assertEqual(repr(ident), "uniprot:P04637")
        self.assertEqual(ident.type, Namespace.UNIPROT)
        self.assertEqual(ident.type, "uniprot")

        ann = Annotation(term=slots.direction_qualifier, value="upregulated")
        self.assertEqual(repr(ann), "direction_qualifier=upregulated")

        ent = Entity(type=Protein, identifiers=[ident], annotations=[ann])
        self.assertIn("protein", str(ent))
        self.assertEqual(len(ent.identifiers), 1)

        rel = Relation(subject=ent, predicate=slots.regulates, object=ent)
        self.assertIn("protein", str(rel))
        self.assertIn("regulates", str(rel))

    def test_normalize_namespace(self):
        self.assertEqual(normalize_namespace("pubchem"), Namespace.PUBCHEM)
        self.assertEqual(normalize_namespace("ensg"), Namespace.ENSG)
        self.assertEqual(normalize_namespace("uniprot"), Namespace.UNIPROT)
        self.assertEqual(normalize_namespace("complexportal"), Namespace.COMPLEXPORTAL)
        self.assertEqual(normalize_namespace("genesymbol"), Namespace.GENESYMBOL)
        self.assertEqual(normalize_namespace("inchikey"), Namespace.INCHIKEY)
        self.assertEqual(normalize_namespace(Namespace.CHEBI), Namespace.CHEBI)
        self.assertEqual(normalize_namespace("custom_db"), "custom_db")

    def test_serving_schemas(self):
        self.assertIn("entity_key", ENTITY_TABLE.names)
        self.assertIn("label", ENTITY_TABLE.names)
        self.assertIn("relation_key", RELATION_TABLE.names)
        self.assertIn("payload_json", PAYLOAD_SCHEMA.names)
        self.assertIn("type", SILVER_ENTITY_SCHEMA.names)
        self.assertIn("predicate", SILVER_RELATION_SCHEMA.names)

    def test_keys(self):
        k1 = entity_key("protein", "uniprot", "P04637")
        k2 = entity_key("protein", "uniprot", "P04637")
        self.assertEqual(k1, k2)
        rk = relation_key(k1, "regulates", k2)
        self.assertEqual(len(rk), 64)

    def test_versioning(self):
        self.assertEqual(validate_version("1.0.0"), "1.0.0")
        with self.assertRaises(ValueError):
            validate_version("invalid-version")
        self.assertEqual(validate_source("signor"), "signor")
        with self.assertRaises(ValueError):
            validate_source("Bad Source!")


if __name__ == "__main__":
    unittest.main()
