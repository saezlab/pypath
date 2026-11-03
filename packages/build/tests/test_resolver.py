"""EntityResolver: the pipeline-facing wrapper around the library matcher."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from library_fixture import WATER, build_fixture_library
from omnipath_resolver import EntityResolver, locate_library_dir
from omnipath_build.silver import RawEntityObservation


def obs(key, entity_type, ns, ident, extra=(), taxon="9606"):
    return RawEntityObservation(
        entity_key=key,
        entity_type=entity_type,
        namespace=ns,
        identifier=ident,
        taxon=taxon,
        identifiers=[{"ns": a, "id": b} for a, b in extra],
    )


class TestEntityResolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.library = build_fixture_library(Path(cls.tmp.name))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_resolves_genes_and_chemicals_in_one_batch(self):
        resolver = EntityResolver(library_dir=self.library)
        try:
            out = resolver.resolve_entities(
                {
                    "p": obs("p", "protein", "entrez", "7157"),
                    "c": obs("c", "chemical_entity", "hmdb", "HMDB0002111"),
                    "kw": obs("kw", "ontology_class", "cv_term_accession", "KW-9993", taxon=None),
                    "x": obs("x", "protein", "uniprot", "Q00000"),
                },
                progress=False,
            )
            self.assertEqual(
                (out["p"].canonical_namespace, out["p"].canonical_identifier), ("entrez", "7157")
            )
            self.assertEqual(out["p"].label, "TP53")
            self.assertTrue(out["p"].matched)
            self.assertEqual(out["c"].canonical_identifier, WATER)
            self.assertEqual(out["c"].label.lower(), "water")
            self.assertEqual(
                (out["kw"].canonical_namespace, out["kw"].label), ("uniprot_keyword", "Ligand")
            )
            self.assertFalse(out["x"].matched)
            stats = resolver.resolution_stats()
            self.assertEqual(stats["input_entities"], 4)
            self.assertEqual(stats["resolved_entities"], 2)
            self.assertEqual(stats["unresolved_entities"], 1)
            self.assertEqual(stats["not_applicable_entities"], 1)
            self.assertEqual(stats["by_entity_type"]["protein"]["resolved"], 1)
            self.assertEqual(sorted(stats["libraries"]), ["chemical", "gene_protein", "reaction"])
        finally:
            resolver.close()

    def test_observation_aliases_do_not_leak_into_later_requests(self):
        resolver = EntityResolver(library_dir=None)
        try:
            first = obs("same", "protein", "uniprot", "X00001")
            later = obs("same", "protein", "uniprot", "X00001", [("genesymbol", "LATE")])
            resolver.resolve_entities({"same": first}, progress=False)
            info = resolver.resolve_entities({"same": later}, progress=False)["same"]
            self.assertEqual(info.label, "LATE")
            self.assertEqual(info.aliases["genesymbol"], ["LATE"])
            repeated = resolver.resolve_entities({"same": first}, progress=False)["same"]
            self.assertNotIn("LATE", repeated.aliases.get("genesymbol", []))
            stats = resolver.resolution_stats()
            self.assertEqual(stats["input_entities"], 1)
            self.assertEqual(stats["unresolved_entities"], 1)
            self.assertEqual(stats["libraries"], [])
        finally:
            resolver.close()

    def test_without_library_everything_keeps_native_ids(self):
        resolver = EntityResolver(library_dir=None)
        try:
            info = resolver.resolve_entities(
                {"p": obs("p", "protein", "entrez", "7157")}, progress=False
            )["p"]
            self.assertEqual(
                (info.canonical_namespace, info.canonical_identifier), ("entrez", "7157")
            )
            self.assertEqual(info.resolved_by, "unmatched")
        finally:
            resolver.close()

    def test_locate_library_dir(self):
        with mock.patch.dict(os.environ, {"OMNIPATH_LIBRARY_DIR": "/env/library"}):
            self.assertEqual(locate_library_dir(None, "/data"), Path("/env/library"))
            self.assertEqual(locate_library_dir("/explicit", "/data"), Path("/explicit"))
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(locate_library_dir(None, "/data"), Path("/data/reference/library"))
            self.assertIsNone(locate_library_dir(None, None))


if __name__ == "__main__":
    unittest.main()
