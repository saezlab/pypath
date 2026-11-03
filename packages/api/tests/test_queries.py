"""Unit tests for DuckDB ParquetQueries."""

import tempfile
import time
import unittest
from pathlib import Path

from omnipath_api.queries.parquet import ParquetQueries

from omnipath_core.fixtures import write_dataset
from omnipath_core.fixtures import write_resource_version as _write_resource


class TestParquetQueries(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.engine = ParquetQueries(data_root=write_dataset(cls.tmp.name))
        cls.addClassCleanup(cls.engine.close)

    def test_01_list_resources(self):
        resources = self.engine.list_resources()
        self.assertIsInstance(resources, list)
        self.assertGreater(len(resources), 0, "Expected at least one resource indexed")

    def test_02_search_entities(self):
        result = self.engine.search_entities_api(query="TP53", limit=10)
        self.assertEqual(result["entities"][0]["label"], "TP53")
        self.assertIn("elapsed_ms", result)

    def test_03_search_relations(self):
        result = self.engine.search_relations(filters={}, limit=10)
        self.assertIn("rows", result)
        self.assertIn("total", result)
        self.assertIn("elapsed_ms", result)
        self.assertGreater(result["total"], 0, "Expected relations in dataset")

        first = result["rows"][0]
        self.assertIn("subject_label", first)
        self.assertIn("predicate", first)
        self.assertIn("object_label", first)
        # Verify lightweight projection by default
        self.assertNotIn("evidence", first)

        # Verify details when explicitly requested
        detailed = self.engine.search_relations(filters={}, limit=1, include_details=True)
        self.assertIn("evidence", detailed["rows"][0])

    def test_04_facet_counts(self):
        facets = self.engine.get_facets(filters={})
        self.assertIn("categories", facets)
        self.assertIn("predicates", facets)
        self.assertIn("signs", facets)
        self.assertIn("elapsed_ms", facets)

    def test_05_export_parquet(self):
        data, mime = self.engine.export_slice(filters={"genes": ["TP53"]})
        self.assertEqual(mime, "application/vnd.apache.parquet")
        self.assertGreater(len(data), 0)

    def test_07_svelte_entity_search(self):
        sample = self.engine.search_entities_api(limit=1)
        self.assertTrue(sample["entities"])
        test_query = sample["entities"][0]["label"]
        result = self.engine.search_entities_api(query=test_query, limit=10)
        self.assertIn("entities", result)
        self.assertGreater(len(result["entities"]), 0)
        first = result["entities"][0]
        self.assertIn("entityPk", first)
        self.assertIn("canonicalIdentifier", first)

    def test_08_svelte_relation_search(self):
        result = self.engine.search_relations_api(filters={}, limit=10)
        self.assertIn("relations", result)
        self.assertIn("rows", result)
        self.assertNotIn("records", result)
        self.assertGreater(result["total"], 0)
        self.assertIn("relationPk", result["relations"][0])
        self.assertIn("subjectEntity", result["rows"][0])
        self.assertIn("entityPk", result["rows"][0]["subjectEntity"])

    def test_entity_search_is_scalar_first(self):
        sample = self.engine.search_entities_api(limit=1)
        self.assertTrue(sample["entities"])
        query = sample["entities"][0]["label"]
        result = self.engine.search_entities_api(query=query, limit=5)
        self.assertTrue(result["entities"])
        first = result["entities"][0]
        self.assertIsNone(first.get("entityAttributes"))
        self.assertLessEqual(len(first.get("identifiers") or []), 4)
        self.assertLess(result["elapsed_ms"], 1500)

    def test_entity_search_ranks_exact_label_first(self):
        result = self.engine.search_entities_api(query="ADA", limit=5)
        self.assertTrue(result["entities"])
        self.assertEqual(result["entities"][0].get("label"), "ADA")

    def test_entity_search_nested_identifier_fallback(self):
        result = self.engine.search_entities_api(query="P00813", limit=5)
        self.assertTrue(result["entities"])
        labels = {str(item.get("label") or "") for item in result["entities"]}
        identifiers = {str(item.get("canonicalIdentifier") or "") for item in result["entities"]}
        self.assertTrue("ADA" in labels or "100" in identifiers)

    def test_scoped_entity_facets_include_source_counts(self):
        facets = self.engine.get_scoped_entity_facets({"query": "ADA", "facetLimit": 5})
        sources = [item for item in facets if item["facetName"] == "source"]
        self.assertTrue(sources)
        self.assertTrue(any(item["scopedCount"] > 0 for item in sources))

    def test_entity_search_respects_source_filter(self):
        all_results = self.engine.search_entities_api(query="ADA", limit=5)
        self.assertTrue(all_results["entities"])
        for resource in ("signor", "uniprot", "chebi"):
            scoped = self.engine.search_entities_api(
                query="ADA", limit=5, filters={"sources": [resource]}
            )
            self.assertLessEqual(len(scoped["entities"]), len(all_results["entities"]))

    def test_a_source_filter_reads_only_that_resources_relations(self):
        # The same relations and counts as a filter over every file, from that file alone.
        resources = sorted(
            {info["resource"] for info in self.engine._selected_resource_infos(None)}
        )
        self.assertGreater(len(resources), 1)
        for resource in resources:
            filters = {"sources": [resource]}
            selection, _ = self.engine._relation_selection(filters)
            self.assertIn(f"/{resource}/", selection)
            self.assertEqual(selection.count("relation.parquet"), 1)
            scoped = self.engine.search_relations(filters=filters, limit=1000)
            full = self.engine.search_relations(filters=filters, resources=resources, limit=1000)
            def key(rows):
                return sorted((r["resource"], r["relation_id"]) for r in rows)

            self.assertEqual(key(scoped["rows"]), key(full["rows"]))
            self.assertEqual(scoped["total"], full["total"])
            self.assertEqual(
                self.engine._compute_relation_facets({"sources": [resource]}),
                self.engine._compute_relation_facets({"sources": [resource]}, resources),
            )
        # A source the release does not have matches nothing.
        self.assertEqual(self.engine.search_relations(filters={"sources": ["nosuch"]})["rows"], [])

    def test_partner_types_select_pairs_of_endpoint_types(self):
        rows = self.engine.search_relations(limit=1000)["rows"]
        pairs = {(r["subject_type"], r["object_type"]) for r in rows}
        self.assertTrue(pairs)
        options = {
            (f["facetName"], f["facetValue"]) for f in self.engine.get_scoped_relation_facets({})
        }
        for first, second in sorted(pairs):
            filters = {"interaction_types": [first], "partner_types": [second]}
            found = self.engine.search_relations(filters=filters, limit=1000)
            # Either direction, never another type at either end.
            self.assertTrue(found["rows"])
            for row in found["rows"]:
                self.assertIn(
                    (row["subject_type"], row["object_type"]), {(first, second), (second, first)}
                )
            # Both types are options, and so is their pair: what choosing both keeps.
            self.assertIn(("participant_pair", "|".join(sorted((first, second)))), options)
            self.assertIn(("participant_type", first), options)
        # The pair options under another group's filter are the pairs of its relations.
        some = sorted(pairs)[0]
        predicate = self.engine.search_relations(
            filters={"interaction_types": [some[0]], "partner_types": [some[1]]}, limit=1
        )["rows"][0]["predicate"]
        kept = {
            tuple(sorted((r["subject_type"], r["object_type"])))
            for r in self.engine.search_relations(filters={"predicates": [predicate]}, limit=1000)[
                "rows"
            ]
            if r["subject_type"] and r["object_type"]
        }
        offered = {
            tuple(f["facetValue"].split("|"))
            for f in self.engine.get_scoped_relation_facets(
                {
                    "predicates": [predicate],
                    "interactionTypes": [some[0]],
                    "partnerTypes": [some[1]],
                }
            )
            if f["facetName"] == "participant_pair"
        }
        self.assertEqual(offered, kept)

    def test_entity_key_lookup_skips_identifier_scan(self):
        fake = "a" * 64
        self.assertEqual(self.engine.resolve_entity_keys([fake]), [fake])
        sample = self.engine.search_relations(limit=1)
        self.assertTrue(sample.get("rows"))
        key = sample["rows"][0]["subject_entity_key"]
        found = self.engine.get_entities_by_pks([key], slim=True)["entities"]
        self.assertTrue(
            any(
                item["entityPk"] == key or key in (item.get("sourceEntityPks") or [])
                for item in found
            )
        )

    def test_empty_relation_facets_are_cached(self):
        self.engine._facet_cache.clear()
        first = self.engine.get_scoped_relation_facets({})
        self.assertTrue(self.engine._facet_cache.entries)
        second = self.engine.get_scoped_relation_facets({})
        self.assertEqual(first, second)

    def test_concurrent_relation_queries(self):
        from concurrent.futures import ThreadPoolExecutor

        def run(_):
            return self.engine.search_relations(limit=5)

        with ThreadPoolExecutor(4) as pool:
            results = list(pool.map(run, range(8)))
        self.assertTrue(all(item.get("rows") for item in results))

    def test_09_scoped_relation_facets(self):
        sample_rel = self.engine.search_relations(limit=1)
        self.assertTrue(sample_rel.get("rows"))
        entity_id = sample_rel["rows"][0]["subject_label"]
        counts = self.engine.get_scoped_relation_facets({"entityIds": [entity_id]})
        self.assertTrue(any(item["facetName"] == "predicate" for item in counts))

    def test_selection_options_follow_filters_from_kept_combinations(self):
        self.engine._facet_cache.clear()
        key = self.engine.search_relations(limit=1)["rows"][0]["subject_entity_key"]
        scope = {"entityIds": [key]}
        options = self.engine.get_scoped_relation_facets(scope)
        predicates = [o["facetValue"] for o in options if o["facetName"] == "predicate"]
        self.assertTrue(predicates)
        for predicate in predicates:
            filters = {"entity_pks": [key], "predicates": [predicate]}
            rows = self.engine.search_relations(filters=filters, limit=1000)["rows"]
            self.assertTrue(rows)
            narrowed = self.engine.get_scoped_relation_facets({**scope, "predicates": [predicate]})
            self.assertEqual(
                {o["facetValue"] for o in narrowed if o["facetName"] == "source"},
                {source for row in rows for source in row["sources"]},
            )
            # The relation types stay listed: choosing another one adds it.
            self.assertEqual(
                {o["facetValue"] for o in narrowed if o["facetName"] == "predicate"},
                set(predicates),
            )
        # The selection's combinations were read once and kept.
        kept = [k for k in self.engine._facet_cache.entries if k[0] == "relation-combinations"]
        self.assertEqual(len(kept), 1)

    def test_10_relation_evidence(self):
        result = self.engine.search_relations(limit=1, include_details=True)
        self.assertTrue(result.get("rows"))
        pk = result["rows"][0]["relation_key"]
        evidence = self.engine.get_relation_evidence(pk)
        self.assertIn("evidence", evidence)
        self.assertIn("annotations", evidence)

    def test_11_ontology_hierarchy(self):
        rows = self.engine.search_entities_api(query="KW-0001", limit=100)["entities"]
        matching = [row for row in rows if row["canonicalIdentifier"] == "KW-0001"]
        self.assertTrue(matching)
        details = self.engine.get_entity_details(matching[0]["entityPk"])
        self.assertIsNotNone(details)
        ent = details["entity"]
        self.assertEqual(ent["entityPk"], matching[0]["entityPk"])
        hierarchy = ent.get("ontologyHierarchy")
        self.assertIsNotNone(hierarchy)
        self.assertEqual(hierarchy["termId"], "KW-0001")
        # The local resource release can change; this term must still have ancestry.
        self.assertGreater(hierarchy["parentCount"], 0)
        self.assertGreaterEqual(hierarchy["childCount"], 0)


class TestInventoryHotReload(unittest.TestCase):
    def test_picks_up_new_version_without_reload_call(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_resource(root, "uniprot", "v1", n_relations=1)
            engine = ParquetQueries(data_root=root)
            self.assertEqual(engine._latest_by_resource["uniprot"]["version"], "v1")
            (path,) = engine._table_paths("relation", ["uniprot"]).values()
            self.assertTrue(path.endswith("/v1/relation.parquet"))

            time.sleep(0.05)
            _write_resource(root, "uniprot", "v2", n_relations=3)
            infos = engine._selected_resource_infos(["uniprot"])
            self.assertEqual(infos[0]["version"], "v2")
            self.assertTrue(str(infos[0]["tables"]["relation"]).endswith("/v2/relation.parquet"))
            listed = {item["key"]: item["relations_count"] for item in engine.list_resources()}
            self.assertEqual(listed["uniprot/v2"], 3)

    def test_catalog_is_kept_until_the_data_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_resource(root, "uniprot", "v1", n_relations=1)
            engine = ParquetQueries(data_root=root)
            builds = []
            build = engine._build_resource_catalog
            engine._build_resource_catalog = lambda: builds.append(1) or build()
            first = engine.list_resource_catalog()
            self.assertEqual(engine.list_resource_catalog(), first)
            self.assertEqual(len(builds), 1)
            # A caller's change does not reach the kept catalog.
            first[0]["version"] = "changed"
            self.assertEqual(engine.list_resource_catalog()[0]["version"], "v1")
            time.sleep(0.05)
            _write_resource(root, "uniprot", "v2", n_relations=3)
            catalog = engine.list_resource_catalog()
            self.assertEqual(len(builds), 2)
            self.assertEqual(catalog[0]["version"], "v2")
            self.assertEqual(catalog[0]["interaction_count"], 3)
            # So does a change of the query policy.
            (root / "query_policy.json").write_text('{"disabled_resources": ["uniprot"]}')
            self.assertFalse(engine.list_resource_catalog()[0]["queryable"])
            self.assertEqual(len(builds), 3)

    def test_topic_map_lists_each_topics_resources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("signor", "uniprot", "chebi"):
                _write_resource(root, name, "v1", n_relations=1)
            engine = ParquetQueries(data_root=root)
            catalog = engine.list_resource_catalog()
            expected = {}
            for resource in catalog:
                for tag in resource.get("tags") or []:
                    if tag["dimension"] == "topic":
                        expected.setdefault(tag["id"], set()).add(resource["resource_id"])
            self.assertTrue(expected, "the fixture resources have topics")
            topics = engine.list_resource_topics()
            self.assertEqual({t["id"]: set(t["resources"]) for t in topics}, expected)
            self.assertTrue(
                all(t["label"] and t["resources"] == sorted(t["resources"]) for t in topics)
            )
            sizes = [len(t["resources"]) for t in topics]
            self.assertEqual(sizes, sorted(sizes, reverse=True))

    def test_skips_version_with_a_missing_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_resource(root, "uniprot", "v1", n_relations=1)
            engine = ParquetQueries(data_root=root)
            time.sleep(0.05)
            pending = _write_resource(root, "uniprot", "v2", n_relations=9)
            serving = pending / "entity_term.parquet"
            serving.rename(pending / "entity_term.partial")
            infos = engine._selected_resource_infos(["uniprot"])
            self.assertEqual(infos[0]["version"], "v1")
            (pending / "entity_term.partial").rename(serving)
            infos = engine._selected_resource_infos(["uniprot"])
            self.assertEqual(infos[0]["version"], "v2")

    def test_picks_up_in_place_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_resource(root, "uniprot", "v1", n_relations=1)
            engine = ParquetQueries(data_root=root)
            self.assertEqual(engine.list_resources()[0]["relations_count"], 1)
            time.sleep(0.05)
            _write_resource(root, "uniprot", "v1", n_relations=4)
            self.assertEqual(engine.list_resources()[0]["relations_count"], 4)


if __name__ == "__main__":
    unittest.main()


def test_catalog_exposes_sample_limit_and_partial_build(tmp_path):
    import json

    folder = _write_resource(tmp_path, "sample", "1")
    (folder / "build_manifest.json").write_text(json.dumps({"max_records": 10, "datasets": None}))
    engine = ParquetQueries(data_root=tmp_path)
    item = engine.list_resource_catalog()[0]
    assert item["sample_build"] is True
    assert item["max_records"] == 10
    assert engine.get_stats_build_manifest()["partialBuild"] is True
