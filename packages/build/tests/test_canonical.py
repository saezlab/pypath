"""Library builder + matcher: the canonicalization rules, one case per rule."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


from omnipath_resolver.canonical import (
    CHEMICAL_POLICY,
    GENE_PROTEIN_POLICY,
    LibraryMatcher,
    ResolutionTracker,
    assign_label,
    get_policy,
    normalize_identifier,
)
from omnipath_build.silver import RawEntityObservation

from library_fixture import ASPIRIN, LINKED, MIX_A, WATER, build_fixture_library


def obs(entity_type, ns, ident, extra=(), taxon="9606", key=None):
    return RawEntityObservation(
        entity_key=key or f"{entity_type}|{ns}|{ident}|{taxon}|{len(extra)}",
        entity_type=entity_type,
        namespace=ns,
        identifier=ident,
        taxon=taxon,
        identifiers=[{"ns": a, "id": b} for a, b in extra],
    )


class LibraryTestCase(unittest.TestCase):
    tmp: tempfile.TemporaryDirectory
    library: Path
    matcher: LibraryMatcher

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.library = build_fixture_library(Path(cls.tmp.name))
        cls.matcher = LibraryMatcher(cls.library)

    @classmethod
    def tearDownClass(cls):
        cls.matcher.close()
        cls.tmp.cleanup()

    def one(self, observation):
        return self.matcher.match({observation.entity_key: observation})[observation.entity_key]


class TestPolicy(unittest.TestCase):
    def test_policy_lookup(self):
        self.assertIs(get_policy("protein"), GENE_PROTEIN_POLICY)
        self.assertIs(get_policy("gene"), GENE_PROTEIN_POLICY)
        self.assertIs(get_policy("small_molecule"), CHEMICAL_POLICY)
        self.assertEqual(get_policy("ontology_class").entity_class, "cv_term")
        self.assertFalse(get_policy("ontology_class").matched)
        with self.assertRaises(ValueError):
            get_policy("some_random_type")

    def test_normalize_identifier(self):
        self.assertEqual(normalize_identifier("ChEBI", "15377"), [("chebi", "CHEBI:15377")])
        self.assertEqual(normalize_identifier("hgnc", "11998"), [("hgnc", "HGNC:11998")])
        self.assertEqual(
            normalize_identifier("inchikey", "InChIKey=" + WATER.lower()), [("inchikey", WATER)]
        )
        self.assertEqual(
            normalize_identifier("ensg", "ENSG00000141510.15"), [("ensg", "ENSG00000141510")]
        )
        self.assertEqual(
            normalize_identifier("uniprot", "P04637-2"),
            [("uniprot", "P04637-2"), ("uniprot", "P04637")],
        )
        self.assertEqual(normalize_identifier("uniprot", ""), [])
        self.assertEqual(
            normalize_identifier("kegg_gene", "mmu:16348"),
            [("kegg_gene", "mmu:16348"), ("entrez", "16348")],
        )
        # Other organisms use their own gene names, not NCBI Gene IDs.
        self.assertEqual(
            normalize_identifier("kegg_gene", "dme:Dmel_CG1234"), [("kegg_gene", "dme:Dmel_CG1234")]
        )
        self.assertEqual(
            normalize_identifier("uniprot_trembl", "A0A0U1RQF1-2"),
            [("uniprot", "A0A0U1RQF1-2"), ("uniprot", "A0A0U1RQF1")],
        )

    def test_assign_label(self):
        self.assertEqual(
            assign_label("entrez", "7157", {"genesymbol": ["TP53"]}, get_policy("protein")), "TP53"
        )
        self.assertEqual(
            assign_label("inchikey", "ABC", {"name": ["Aspirin"]}, get_policy("chemical_entity")),
            "Aspirin",
        )
        self.assertEqual(assign_label("kegg", "C1", {}, get_policy("chemical_entity")), "C1")


class TestLibraryBuild(LibraryTestCase):
    def test_complete_identity_library(self):
        import json

        manifest = json.loads((self.library / "manifest.json").read_text())
        self.assertEqual(manifest["format"], "omnipath-identity-v2")
        self.assertTrue(manifest["fingerprint"])
        self.assertIn("uniprot", manifest["hub_indexes"])
        self.assertTrue((self.library / "kv" / "manifest.json").is_file())


class TestGeneProteinMatching(LibraryTestCase):
    def test_every_identifier_of_a_gene_hits_the_same_node(self):
        for ns, ident in [
            ("uniprot", "P04637"),
            ("uniprot", "P04637-2"),
            ("uniprot-sec", "Q15086"),
            ("entrez", "7157"),
            ("ensg", "ENSG00000141510.15"),
            ("hgnc", "11998"),
            ("kegg_gene", "hsa:7157"),
        ]:
            m = self.one(obs("protein", ns, ident))
            self.assertEqual(
                (m.canonical_namespace, m.canonical_identifier), ("entrez", "7157"), (ns, ident)
            )
            self.assertEqual(m.label, "TP53")
            self.assertEqual(m.resolved_by, "parquet")
            self.assertEqual(m.taxon, "9606")
            self.assertIn("7157", m.aliases["entrez"])

    def test_gene_reference_retains_source_type_without_product_fanout(self):
        for etype in ("gene", "rna_product"):
            m = self.one(obs(etype, "uniprot", "P04637"))
            self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("entrez", "7157"))
            self.assertEqual(m.label, "TP53")
        m = self.one(obs("protein", "entrez", "7157"))
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("entrez", "7157"))

    def test_symbol_needs_taxon_and_is_taxon_bound(self):
        human = self.one(obs("protein", "genesymbol", "TP53"))
        self.assertEqual(human.canonical_identifier, "7157")
        self.assertEqual(human.resolved_by, "parquet")
        mouse = self.one(obs("protein", "genesymbol", "TP53", taxon="10090"))
        self.assertEqual(mouse.canonical_identifier, "22059")  # via synonym
        self.assertEqual(mouse.label, "Trp53")
        synonym = self.one(obs("protein", "genesymbol", "P53"))
        self.assertEqual(synonym.canonical_identifier, "7157")
        no_taxon = self.one(obs("protein", "genesymbol", "TP53", taxon=None))
        self.assertFalse(no_taxon.matched)
        self.assertEqual(
            (no_taxon.canonical_namespace, no_taxon.canonical_identifier), ("genesymbol", "TP53")
        )

    def test_stale_nonmatching_evidence_is_ignored_but_conflicts_stay_unresolved(self):
        stale = self.one(obs("protein", "uniprot", "P04637", [("genesymbol", "WRONG")]))
        self.assertEqual(stale.canonical_identifier, "7157")
        self.assertEqual(stale.resolved_by, "parquet")
        conflict = self.one(obs("protein", "uniprot", "P04637", [("entrez", "55")]))
        self.assertEqual(conflict.canonical_identifier, "P04637")
        self.assertTrue(conflict.matched)
        self.assertEqual(conflict.gene_mapping_status, "conflict")
        self.assertIn("55", conflict.aliases["entrez"])  # observed ids are kept as aliases

    def test_calmodulin_and_multi_reviewed(self):
        calm = self.one(obs("protein", "uniprot", "P0DP23"))
        self.assertEqual(
            (calm.canonical_namespace, calm.canonical_identifier), ("uniprot", "P0DP23")
        )
        self.assertCountEqual(calm.aliases["entrez"], ["801", "805", "808"])
        gene = self.one(obs("gene", "entrez", "801"))
        self.assertEqual((gene.canonical_namespace, gene.canonical_identifier), ("entrez", "801"))
        self.assertEqual(gene.label, "801")
        two = self.one(obs("protein", "uniprot", "P11111"))
        self.assertEqual((two.canonical_namespace, two.canonical_identifier), ("entrez", "999"))
        self.assertEqual(two.protein_identifier, "P11111")
        self.assertEqual(two.protein_aliases["uniprot"], ["P11111"])

    def test_unmatched_keeps_native_identifier(self):
        m = self.one(obs("protein", "uniprot", "Q00000", [("name", "Some protein")]))
        self.assertFalse(m.matched)
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("uniprot", "Q00000"))
        self.assertEqual(m.label, "Some protein")
        self.assertEqual(m.resolved_by, "unmatched")
        yeast = self.one(obs("protein", "uniprot", "Q99999", taxon="4932"))
        self.assertTrue(yeast.matched)
        hiv = self.one(obs("protein", "uniprot", "P12345", taxon="11676"))
        self.assertTrue(hiv.matched)
        self.assertEqual((hiv.canonical_namespace, hiv.canonical_identifier), ("uniprot", "P12345"))

    def test_generic_types_are_not_matched(self):
        m = self.one(obs("microrna", "mirbase", "MI0000001", [("name", "hsa-mir-1")]))
        self.assertFalse(m.matched)
        self.assertEqual(m.canonical_namespace, "mirbase")
        self.assertEqual(m.label, "hsa-mir-1")


class TestResolutionMemo(LibraryTestCase):
    def test_recurring_observation_reuses_the_result(self):
        first = self.matcher.targets({"a": obs("protein", "uniprot", "P04637", key="a")})["a"]
        hits = self.matcher.metrics["memo_hits"]
        # The same observation under another key, in a later batch.
        again = self.matcher.targets({"b": obs("protein", "uniprot", "P04637", key="b")})["b"]
        self.assertEqual(self.matcher.metrics["memo_hits"], hits + 1)
        self.assertEqual(again, first)
        self.assertIsNot(again[0], first[0])
        # A different observation is resolved, not served from the memo.
        other = self.matcher.targets({"c": obs("protein", "entrez", "7157", key="c")})["c"]
        self.assertEqual(self.matcher.metrics["memo_hits"], hits + 1)
        self.assertEqual(other[0].canonical_identifier, "7157")


class TestChemicalMatching(LibraryTestCase):
    def test_every_identifier_of_water(self):
        for ns, ident in [
            ("chebi", "CHEBI:15377"),
            ("chebi", "15377"),
            ("hmdb", "HMDB0002111"),
            ("kegg", "C00001"),
            ("cas", "7732-18-5"),
            ("bigg", "h2o"),
            ("bigg_metabolite", "h2o"),
            ("bigg_metabolite", "h2o_c"),
            ("metanetx", "MNXM2"),
            ("inchikey", "InChIKey=" + WATER),
        ]:
            m = self.one(obs("chemical_entity", ns, ident))
            self.assertEqual(
                (m.canonical_namespace, m.canonical_identifier), ("inchikey", WATER), (ns, ident)
            )
            self.assertEqual(m.label.lower(), "water")
            self.assertIn("HMDB0002111", m.aliases["hmdb"])

    def test_model_identifiers_vote_as_secondary_ids(self):
        for ns, ident in [("bigg_metabolite", "h2o"), ("metanetx", "MNXM2")]:
            match = self.one(obs("chemical_entity", "imm1415", "local_water", [(ns, ident)]))
            self.assertEqual(match.canonical_identifier, WATER)
            self.assertEqual(match.resolved_by, "parquet")
            self.assertIn("h2o", match.aliases["bigg"])
            self.assertIn("MNXM2", match.aliases["metanetx"])

    def test_supplied_structure_overrides_model_ids(self):
        match = self.one(
            obs(
                "chemical_entity",
                "imm1415",
                "unknown",
                [("metanetx", "MNXM2"), ("inchikey", ASPIRIN)],
            )
        )
        self.assertTrue(match.matched)
        self.assertEqual(match.canonical_identifier, ASPIRIN)
        # Explicit structure also outranks a unique primary model ID.
        match = self.one(obs("chemical_entity", "bigg_metabolite", "h2o", [("inchikey", ASPIRIN)]))
        self.assertEqual(match.canonical_identifier, ASPIRIN)
        self.assertEqual(match.resolved_by, "parquet")

    def test_reaction_ids_are_not_chemical_votes(self):
        for ns in ("bigg_reaction", "metanetx_reaction"):
            self.assertFalse(self.one(obs("chemical_entity", ns, "h2o")).matched)
        self.assertEqual(normalize_identifier("bigg_metabolite", "h2o_c"), [("bigg", "h2o_c")])

    def test_agreement_and_conflict(self):
        agree = self.one(obs("chemical_entity", "chebi", "CHEBI:15365", [("chembl", "CHEMBL25")]))
        self.assertEqual(agree.canonical_identifier, ASPIRIN)
        self.assertEqual(agree.resolved_by, "parquet")
        conflict = self.one(
            obs("chemical_entity", "chebi", "CHEBI:15377", [("chembl", "CHEMBL25")])
        )
        self.assertEqual(conflict.canonical_identifier, "CHEBI:15377")
        self.assertFalse(conflict.matched)

    def test_chemical_canonical_namespaces(self):
        m = self.one(
            obs(
                "chemical_entity",
                "chembl",
                "CHEMBL25",
                [("inchikey", ASPIRIN)],
            )
        )
        self.assertTrue(m.matched)
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("inchikey", ASPIRIN))
        self.assertEqual(m.resolved_by, "parquet")

    def test_structureless_class_is_its_own_node(self):
        m = self.one(obs("chemical_entity", "chebi", "CHEBI:24431"))
        self.assertTrue(m.matched)
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("chebi", "CHEBI:24431"))
        self.assertEqual(m.label, "chemical entity")

    def test_class_folded_into_compound(self):
        m = self.one(obs("chemical_entity", "chebi", "16974"))
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("inchikey", LINKED))
        self.assertIn("DB09145", m.aliases["drugbank"])

    def test_mixture(self):
        mix = self.one(obs("chemical_entity", "chebi", "CHEBI:1"))
        self.assertEqual((mix.canonical_namespace, mix.canonical_identifier), ("chebi", "CHEBI:1"))
        key = self.one(obs("chemical_entity", "inchikey", MIX_A))
        self.assertEqual((key.canonical_namespace, key.canonical_identifier), ("inchikey", MIX_A))

    def test_unmatched_chemical(self):
        m = self.one(obs("chemical_entity", "kegg", "C99999", [("name", "Mystery")]))
        self.assertFalse(m.matched)
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("kegg", "C99999"))
        self.assertEqual(m.label, "Mystery")

    def test_chebi_ontology_terms_are_not_chemicals(self):
        m = self.one(obs("ontology_class", "chebi", "CHEBI:15377", [("name", "water")]))
        self.assertFalse(m.matched)
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("chebi", "CHEBI:15377"))
        self.assertEqual(m.label.lower(), "water")


class TestOtherClasses(LibraryTestCase):
    def test_cv_term_rewrite_and_glossary_label(self):
        kw = self.one(obs("ontology_class", "cv_term_accession", "KW-9993"))
        self.assertEqual(
            (kw.canonical_namespace, kw.canonical_identifier), ("uniprot_keyword", "KW-9993")
        )
        self.assertEqual(kw.label, "Ligand")
        go = self.one(obs("ontology_class", "go", "GO:0005515", [("name", "protein binding")]))
        self.assertEqual(go.canonical_namespace, "go")
        self.assertEqual(go.label, "protein binding")

    def test_complex(self):
        m = self.one(obs("macromolecular_complex", "signor", "SIGNOR-C1", [("name", "AP-1")]))
        self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("signor", "SIGNOR-C1"))
        self.assertEqual(m.label, "AP-1")
        unnamed = self.one(obs("macromolecular_complex", "signor", "SIGNOR-C2"))
        self.assertEqual(unnamed.label, "SIGNOR-C2")

    def test_no_library(self):
        matcher = LibraryMatcher(None)
        try:
            m = matcher.match({"k": obs("protein", "uniprot", "P04637", key="k")})["k"]
            self.assertFalse(m.matched)
            self.assertEqual((m.canonical_namespace, m.canonical_identifier), ("uniprot", "P04637"))
        finally:
            matcher.close()


class TestTracker(LibraryTestCase):
    def test_counts(self):
        tracker = ResolutionTracker()
        items = {
            "a": obs("protein", "uniprot", "P04637", key="a"),
            "b": obs("chemical_entity", "kegg", "C99999", key="b"),
            "c": obs("ontology_class", "go", "GO:1", key="c"),
        }
        matches = self.matcher.match(items)
        for key, o in items.items():
            tracker.record(key, o.entity_type, get_policy(o.entity_type), matches[key])
        tracker.record(
            "a", "protein", get_policy("protein"), matches["a"]
        )  # re-seen, no double count
        s = tracker.summary()
        self.assertEqual(s["input_entities"], 3)
        self.assertEqual(s["resolved_entities"], 1)
        self.assertEqual(s["unresolved_entities"], 1)
        self.assertEqual(s["not_applicable_entities"], 1)
        self.assertEqual(s["by_rule"]["gene_protein/parquet"], 1)

    def test_structure_only_chemicals_and_pending_types(self):
        tracker = ResolutionTracker()
        items = {
            # A full InChIKey the fixture library does not hold.
            "s": obs("chemical_entity", "inchikey", "AAAAAAAAAAAAAA-BBBBBBBBBB-N", key="s"),
            "m": obs("microrna", "mirbase", "MIMAT0000001", key="m"),
        }
        matches = self.matcher.match(items)
        for key, o in items.items():
            tracker.record(key, o.entity_type, get_policy(o.entity_type), matches[key])
        s = tracker.summary()
        self.assertEqual(s["structure_entities"], 1)
        self.assertEqual(s["unresolved_entities"], 0)
        self.assertEqual(s["not_applicable_entities"], 1)
        self.assertEqual(s["lookup_entities"], 1)
        self.assertEqual(s["by_entity_type"]["chemical_entity"]["structure"], 1)


if __name__ == "__main__":
    unittest.main()
