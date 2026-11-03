"""Tests for streaming identifier hub Parquet exports."""

from __future__ import annotations

import gzip
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

import pyarrow.parquet as pq

from omnipath_build.hubs.export import HUB_NAMES, export_hubs
from omnipath_build.hubs.sources import emit_bigg, emit_from_records, emit_metanetx
from omnipath_build.hubs.sources.metanetx_reaction import emit as emit_metanetx_reaction
from omnipath_build.hubs.sources.rhea import emit as emit_rhea
from omnipath_build.hubs.sources.uniprot import emit as emit_uniprot
from omnipath_build.hubs.stream import iter_decoded_lines
from omnipath_build.hubs.writer import HubParquetWriter


class TestHubWriter(unittest.TestCase):
    def test_caps_output_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hub.parquet"
            writer = HubParquetWriter(path, max_records=3)
            for index in range(10):
                writer.add("chebi", str(index), "15377", "0", "chebi")
            self.assertEqual(writer.close(), 3)
            table = pq.read_table(path)
            self.assertEqual(table.num_rows, 3)
            self.assertEqual(
                table.column_names,
                ["source_type", "source_id", "hub_id", "taxonomy_id", "backend"],
            )


class TestHubExplode(unittest.TestCase):
    def test_inputs_v2_style_record_maps_xrefs_to_own_hub(self):
        records = [
            {
                "chebi_id": "CHEBI:15377",
                "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                "pubchem_compound": ["962", "22247451"],
                "alt_ids": ["CHEBI:44819"],
            }
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chebi.parquet"
            writer = HubParquetWriter(path, max_records=20)
            emit_from_records(
                writer,
                records,
                hub_column="chebi_id",
                hub_type="chebi",
                backend="chebi",
                columns={
                    "inchikey": "inchikey",
                    "pubchem": "pubchem_compound",
                },
                extra=lambda row: {"chebi": row.get("alt_ids") or []},
            )
            writer.close()
            rows = pq.read_table(path).to_pylist()
            self.assertIn(
                {
                    "source_type": "chebi",
                    "source_id": "CHEBI:15377",
                    "hub_id": "CHEBI:15377",
                    "taxonomy_id": "0",
                    "backend": "chebi",
                },
                rows,
            )
            self.assertIn(
                {
                    "source_type": "inchikey",
                    "source_id": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                    "hub_id": "CHEBI:15377",
                    "taxonomy_id": "0",
                    "backend": "chebi",
                },
                rows,
            )
            self.assertIn(
                {
                    "source_type": "pubchem",
                    "source_id": "962",
                    "hub_id": "CHEBI:15377",
                    "taxonomy_id": "0",
                    "backend": "chebi",
                },
                rows,
            )
            self.assertIn(
                {
                    "source_type": "chebi",
                    "source_id": "CHEBI:44819",
                    "hub_id": "CHEBI:15377",
                    "taxonomy_id": "0",
                    "backend": "chebi",
                },
                rows,
            )


class TestIdCleaning(unittest.TestCase):
    def test_drops_placeholder_values(self):
        from omnipath_build.hubs.schema import as_values

        self.assertEqual(as_values("none"), [])
        self.assertEqual(as_values("InChI=none"), [])
        self.assertEqual(as_values(["962", "none", ""]), ["962"])


class TestMetanetxStream(unittest.TestCase):
    def test_property_export_adds_structures_and_deduplicates_identity(self):
        properties = [
            "#ID\tname\treference\tformula\tcharge\tmass\tInChI\tInChIKey\tSMILES",
            "MNXM2\twater\tchebi:15377\tH2O\t0\t18.0\tInChI=1S/H2O/h1H2\tXLYOFNOQVPJJNP-UHFFFAOYSA-N\tO",
            "WATER\twater\tchebi:15377\tH2O\t0\t18\tInChI=1S/H2O/h1H2\tXLYOFNOQVPJJNP-UHFFFAOYSA-N\tO",
            "MNXMunknown\tclass\t\tR\t0\t\t\tNA\t-",
        ]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "metanetx.parquet"
            writer = HubParquetWriter(path)
            emit_metanetx(
                writer, lines=["bigg.metabolite:h2o\tMNXM2\twater"], property_lines=properties
            )
            writer.close()
            rows = pq.read_table(path).to_pylist()
            water = [r for r in rows if r["hub_id"] == "MNXM2"]
            self.assertEqual(sum(r["source_type"] == "metanetx" for r in water), 1)
            self.assertEqual(
                {r["source_type"] for r in water},
                {"metanetx", "name", "inchi", "inchikey", "smiles", "bigg"},
            )
            self.assertTrue(
                any(r["hub_id"] == "WATER" and r["source_type"] == "inchikey" for r in rows)
            )
            self.assertFalse(
                any(
                    r["source_type"] in {"inchi", "inchikey", "smiles"}
                    for r in rows
                    if r["hub_id"] == "MNXMunknown"
                )
            )

    def test_maps_every_xref_onto_metanetx_id(self):
        sample = (
            "# comment\n"
            "chebi:15377\tMNXM2\twater\n"
            "kegg.compound:C00001\tMNXM2\twater\n"
            "bigg.metabolite:h2o\tMNXM2\twater\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "metanetx.parquet"
            writer = HubParquetWriter(path, max_records=20)
            emit_metanetx(writer, lines=sample.splitlines())
            writer.close()
            rows = pq.read_table(path).to_pylist()
            self.assertEqual({row["hub_id"] for row in rows}, {"MNXM2"})
            self.assertTrue(
                {"metanetx", "chebi", "kegg", "bigg"} <= {row["source_type"] for row in rows}
            )
            self.assertTrue(all(row["hub_id"] == "MNXM2" for row in rows))


DIRECTIONS = "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n10000\t10001\t10002\t10003\n"


class TestReactionHubs(unittest.TestCase):
    def rows(self, emit, **sources):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hub.parquet"
            writer = HubParquetWriter(path)
            emit(writer, **sources)
            writer.close()
            return {
                (r["source_type"], r["source_id"], r["hub_id"])
                for r in pq.read_table(path).to_pylist()
            }

    def test_rhea_master_records_directional_ids_xrefs_and_equation(self):
        xrefs = (
            "RHEA_ID\tDIRECTION\tMASTER_ID\tID\tDB\n"
            "10000\tUN\t10000\tR00001\tKEGG_REACTION\n"
            "10001\tLR\t10000\tR-HSA-1.3\tREACTOME\n"
            "10000\tUN\t10000\t3.5.1.50\tEC\n"
            "10000\tUN\t10000\tGO:0050126\tGO\n"
        )
        reactions = (
            "ENTRY       RHEA:10000\nDEFINITION  pentanamide + H2O = pentanoate + NH4(+)\n///\n"
            "ENTRY       RHEA:10001\nDEFINITION  pentanamide + H2O => pentanoate + NH4(+)\n///\n"
        )
        rows = self.rows(
            emit_rhea,
            directions=DIRECTIONS.splitlines(),
            xrefs=xrefs.splitlines(),
            reactions=reactions.splitlines(),
        )
        self.assertEqual(
            rows,
            {
                ("rhea", "10000", "10000"),
                ("rhea", "10001", "10000"),
                ("rhea", "10002", "10000"),
                ("rhea", "10003", "10000"),
                ("kegg_reaction", "R00001", "10000"),
                ("reactome", "R-HSA-1.3", "10000"),
                ("ec", "3.5.1.50", "10000"),
                ("name", "pentanamide + H2O = pentanoate + NH4(+)", "10000"),  # master only, no GO
            },
        )

    def test_metanetx_reactions_map_rhea_ids_to_masters(self):
        lines = (
            "#source\tID\tdescription\n"
            "rhea:10002\tMNXR1\tpentanamide amidase\n"
            "rh:10002\tMNXR1\tduplicate short prefix\n"
            "bigg.reaction:AMIDASE\tMNXR1\t\n"
            "biggR:AMIDASE\tMNXR1\tduplicate short prefix\n"
            "vmhreaction:MODEL\tMNXR2\t\n"
            "mnx:EMPTY\tEMPTY\t\n"
        )
        rows = self.rows(
            emit_metanetx_reaction, lines=lines.splitlines(), directions=DIRECTIONS.splitlines()
        )
        self.assertEqual(
            rows,
            {
                ("metanetx_reaction", "MNXR1", "MNXR1"),
                ("rhea", "10000", "MNXR1"),  # the directional id, written as its master
                ("bigg_reaction", "AMIDASE", "MNXR1"),
                ("metanetx_reaction", "MNXR2", "MNXR2"),
                ("vmh_reaction", "MODEL", "MNXR2"),
            },
        )


class TestBiggStream(unittest.TestCase):
    def test_maps_every_xref_onto_bigg_id(self):
        sample = (
            "bigg_id\tuniversal_bigg_id\tname\tmodel_list\tdatabase_links\n"
            "h2o\th2o\tWater\tiJO1366\t"
            "CHEBI: http://identifiers.org/chebi/CHEBI:15377; "
            "KEGG Compound: http://identifiers.org/kegg.compound/C00001; "
            "MetaNetX (MNX) Chemical: http://identifiers.org/metanetx.chemical/MNXM2\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "bigg.parquet"
            writer = HubParquetWriter(path, max_records=20)
            emit_bigg(writer, lines=sample.splitlines(keepends=True))
            writer.close()
            rows = pq.read_table(path).to_pylist()
            self.assertEqual({row["hub_id"] for row in rows}, {"h2o"})
            self.assertTrue(
                {"bigg", "chebi", "kegg", "metanetx", "name"}
                <= {row["source_type"] for row in rows}
            )


class TestUniprotIdmappingStream(unittest.TestCase):
    SAMPLE = (
        "P04637\tUniProtKB-ID\tP53_HUMAN\n"
        "P04637\tGene_Name\tTP53\n"
        "P04637\tGeneID\t7157\n"
        "P04637\tPDB\t1TUP\n"
        "P04637\tEMBL\tX02469\n"
        "P04637\tEMBL-CDS\tCAA26306.1\n"
        "P04637\tDrugBank\tDB00921\n"
        "P04637\tEnsembl\tENSG00000141510.18\n"
        "P04637\tEnsembl_PRO\tENSP00000269305.4\n"
        "P04637\tHGNC\t11998\n"
        "P04637\tNCBI_TaxID\t9606\n"
        "P38398\tGene_Name\tBRCA1\n"
        "P38398\tPDB\t1JM7\n"
        "P38398\tNCBI_TaxID\t9606\n"
    )

    def test_keeps_protein_reference_types_and_drops_structure_xrefs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "uniprot.parquet"
            writer = HubParquetWriter(path, max_records=500)
            emit_uniprot(writer, lines=self.SAMPLE.splitlines(keepends=True))
            writer.close()
            rows = pq.read_table(path).to_pylist()
            types = {row["source_type"] for row in rows}
            self.assertTrue(
                {"uniprot", "uniprot_entry", "genesymbol", "entrez", "ensg", "ensp", "hgnc"}
                <= types
            )
            self.assertFalse(types & {"pdb", "embl", "embl_id", "drugbank"})
            keyed = {(row["source_type"], row["source_id"]): row for row in rows}
            self.assertEqual(keyed[("ensg", "ENSG00000141510")]["taxonomy_id"], "9606")
            self.assertEqual(keyed[("uniprot", "P04637")]["hub_id"], "P04637")
            self.assertEqual(keyed[("genesymbol", "BRCA1")]["hub_id"], "P38398")
            self.assertNotIn(("ensg", "ENSG00000141510.18"), keyed)

    def test_stops_at_max_records_without_secondary_lookup(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "uniprot.parquet"
            writer = HubParquetWriter(path, max_records=3)
            emit_uniprot(writer, lines=self.SAMPLE.splitlines(keepends=True))
            self.assertEqual(writer.close(), 3)


class TestDictionaries(unittest.TestCase):
    def test_writes_dictionary_parquets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            counts = export_hubs(
                tmpdir,
                hubs=[],
                max_records=1,
                include_dictionaries=True,
            )
            self.assertGreater(counts["id_type"], 10)
            self.assertGreater(counts["organism"], 1)
            manifest = json.loads((Path(tmpdir) / "manifest.json").read_text())
            self.assertIn("export_id", manifest)
            self.assertTrue((Path(tmpdir) / "id_type.parquet").exists())


class TestStreamDecode(unittest.TestCase):
    def test_gzip_lines_stop_without_reading_the_rest(self):
        payload = gzip.compress(b"one\ntwo\nthree\nfour\n")
        lines = iter_decoded_lines(io.BytesIO(payload))
        self.assertEqual(next(lines).rstrip("\n"), "one")
        self.assertEqual(next(lines).rstrip("\n"), "two")
        lines.close()

    def test_zip_first_member_lines(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr("inner.txt", "alpha\nbeta\n")
        lines = [line.rstrip("\n") for line in iter_decoded_lines(io.BytesIO(buf.getvalue()))]
        self.assertEqual(lines, ["alpha", "beta"])

    def test_plain_text_lines(self):
        lines = [line.rstrip("\n") for line in iter_decoded_lines(io.BytesIO(b"alpha\nbeta\n"))]
        self.assertEqual(lines, ["alpha", "beta"])


class TestExportProgress(unittest.TestCase):
    def test_dictionary_export_emits_stage_events(self):
        events: list[dict] = []
        with tempfile.TemporaryDirectory() as tmpdir:
            export_hubs(
                tmpdir,
                hubs=[],
                max_records=1,
                include_dictionaries=True,
                on_progress=events.append,
            )
        stages = [(event["stage"], event["status"]) for event in events]
        self.assertIn(("dictionaries", "running"), stages)
        self.assertIn(("dictionaries", "done"), stages)
        self.assertIn(("manifest", "done"), stages)


class TestHubRegistry(unittest.TestCase):
    def test_one_emitter_per_source(self):
        self.assertEqual(
            list(HUB_NAMES),
            [
                "entrez",
                "uniprot",
                "ensg",
                "chebi",
                "pubchem",
                "refmet",
                "ramp",
                "ramp_gene",
                "chembl",
                "hmdb",
                "lipidmaps",
                "swisslipids",
                "bigg",
                "metanetx",
                "mirbase",
                "rhea",
                "metanetx_reaction",
                "taxon_species",
                "ontology",
            ],
        )


class TestEntrezHub(unittest.TestCase):
    def test_gene_info_adds_genes_without_an_ensembl_link(self):
        from omnipath_build.hubs.sources.entrez import emit as emit_entrez

        gene2ensembl = [
            "#tax_id\tGeneID\tEnsembl_gene_identifier\tRNA\tEnsembl_rna\tprotein\tEnsembl_protein",
            "9606\t7157\tENSG00000141510.18\tNM_000546.6\tENST00000269305.9\tNP_000537.3\tENSP00000269305.4",
        ]
        gene_info = [
            "#tax_id\tGeneID\tSymbol\tLocusTag",
            "9606\t7157\tTP53\t-",
            "9606\t100189199\tTRL-CAA6-1\t-",
            "10090\t115487553\tGm25354\t-",
        ]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "entrez.parquet"
            writer = HubParquetWriter(path)
            emit_entrez(writer, lines=gene2ensembl, gene_info=[gene_info])
            writer.close()
            rows = pq.read_table(path).to_pylist()
        records = {(r["hub_id"], r["taxonomy_id"]) for r in rows if r["source_type"] == "entrez"}
        self.assertEqual(records, {("7157", "9606"), ("100189199", "9606"), ("115487553", "10090")})
        # A gene gene2ensembl already lists is not added twice, but gets its symbol too.
        self.assertEqual(
            sum(r["hub_id"] == "7157" and r["source_type"] == "entrez" for r in rows), 1
        )
        symbols = {
            (r["hub_id"], r["source_id"], r["taxonomy_id"])
            for r in rows
            if r["source_type"] == "genesymbol"
        }
        self.assertEqual(
            symbols,
            {
                ("7157", "TP53", "9606"),
                ("100189199", "TRL-CAA6-1", "9606"),
                ("115487553", "Gm25354", "10090"),
            },
        )
