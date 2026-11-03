"""Review regressions checked across parser, resolution and serving boundaries."""

from types import SimpleNamespace

import duckdb
import pyarrow.parquet as pq
import pytest

from omnipath_resolver.canonical.identifiers import normalize_identifier, normalize_id_sql
from writer_fixture import write_observations
from omnipath_resolver import EntityResolver
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from omnipath_core.measurements import quantity_dict
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.inputs_v2 import (
    cellchat,
    cellphonedb,
    chembl,
    connectomedb,
    icellnet,
    rampdb,
    reactome,
    signor,
    swisslipids,
)
from pypath.inputs_v2.base import ontology_term_to_entity
from pypath.inputs_v2.parsers import chembl as chembl_parser
from pypath.internals.ontology_schema import OntologyTerm
from library_fixture import build_fixture_library


def observations(records, source="review"):
    out = SilverExtractor(source, "test")
    for index, record in enumerate(records):
        out.process_record(record, {"index": index}, str(index), index)
    return out


def final(records, path, source="review"):
    out = observations(records, source)
    resolver = EntityResolver(library_dir=path / "absent")
    try:
        return write_observations(resolver, out.entities, out.relations, out.payloads)
    finally:
        resolver.close()


def test_reactome_opposite_controls_survive_separate_parquet_chunks(tmp_path):
    writer = ParquetWriter(tmp_path / "out")
    base = {
        "control_class": "Control",
        "controller_entity_type": "protein",
        "controller_uniprot": "P00533",
        "controlled_entity_type": "reaction",
        "controlled_reactome_stable_id": "R-HSA-123",
    }
    resolver = EntityResolver(library_dir=tmp_path / "absent")
    try:
        for effect in ("ACTIVATION", "INHIBITION"):
            extractor = observations([reactome.controls_schema({**base, "control_type": effect})])
            writer.append_observations(extractor, resolver)
    finally:
        resolver.close()
    writer.close()
    rows = pq.read_table(tmp_path / "out" / "relation.parquet").to_pylist()
    assert len(rows) == 2
    assert {row["sign"] for row in rows} == {-1, 1}
    assert len({row["relation_key"] for row in rows}) == 2


def test_signor_actual_nested_complex_members_are_not_protein_aliases():
    record = signor.complexes_schema(
        {
            "SIGNOR ID": "SIGNOR-C87",
            "COMPLEX NAME": "MLL/SET subcomplex",
            "LIST OF ENTITIES": "Q15291, Q9C005, Q9UPS6, SIGNOR-C352",
        }
    )
    assert len(record.membership) == 4
    assert [[(str(i.type), i.value) for i in m.member.identifiers] for m in record.membership] == [
        [("uniprot", "Q15291")],
        [("uniprot", "Q9C005")],
        [("uniprot", "Q9UPS6")],
        [("signor", "SIGNOR-C352")],
    ]
    assert [format_term(m.member.type) for m in record.membership] == ["protein"] * 3 + [
        "macromolecular_complex"
    ]
    assert len(observations([record]).relations) == 4
    family = signor.protein_families_schema(
        {"SIGNOR ID": "SIGNOR-PF1", "LIST OF ENTITIES": "SIGNOR-PF2,P00533"}
    )
    assert [format_term(m.member.type) for m in family.membership] == ["protein_family", "protein"]


def test_complex_composition_preserves_species_but_merges_sources(tmp_path):
    records = [
        cellchat._protein_group_entity(name="same complex", genes=["A", "B"], taxon_id=t)
        for t in ("9606", "10090")
    ]
    entities = final(records, tmp_path, "cellchat")["entity"]
    groups = [e for e in entities if e["entity_type"] == "macromolecular_complex"]
    assert len(groups) == 2
    assert {g["taxon"] for g in groups} == {"9606", "10090"}
    other = final([records[0]], tmp_path, "another_resource")["entity"]
    assert next(
        e["entity_key"] for e in other if e["entity_type"] == "macromolecular_complex"
    ) == next(g["entity_key"] for g in groups if g["taxon"] == "9606")


def _cellchat_row(ligand, genes, annotation="Non-protein Signaling"):
    return {
        "taxon_id": "9606",
        "annotation": annotation,
        "interaction_name": f"{ligand}_HTR1A",
        "ligand_name": ligand,
        "ligand_genes": genes,
        "receptor_name": "HTR1A",
        "receptor_genes": ["HTR1A"],
    }


def test_cellchat_metabolite_ligand_is_the_metabolite_not_its_genes(tmp_path):
    relation = cellchat.map_cellchat_interaction(
        _cellchat_row("5-HT-TPH1_SLC18A1", ["TPH1", "SLC18A1"])
    )
    assert format_term(relation.subject.type) == "chemical_entity"
    assert (Namespace.CHEBI, "CHEBI:28790") in {
        (i.type, i.value) for i in relation.subject.identifiers
    }
    assert not relation.subject.membership
    # An abbreviation that does not tell the metabolite keeps it as a name.
    ambiguous = cellchat.map_cellchat_interaction(_cellchat_row("Cort-CYP11B1", ["CYP11B1"]))
    assert [(i.type, i.value) for i in ambiguous.subject.identifiers] == [(Namespace.NAME, "Cort")]
    entities = final([relation, ambiguous], tmp_path, "cellchat")["entity"]
    assert sum(e["entity_type"] == "chemical_entity" for e in entities) == 2


def test_cellchat_protein_ligands_stay_proteins():
    peptide = cellchat.map_cellchat_interaction(_cellchat_row("b-Endorphin-POMC", ["POMC"]))
    hormone = cellchat.map_cellchat_interaction(_cellchat_row("CGA_TSHB", ["CGA", "TSHB"]))
    secreted = cellchat.map_cellchat_interaction(
        _cellchat_row("PGE2-PTGES", ["PTGES"], annotation="Secreted Signaling")
    )
    assert format_term(peptide.subject.type) == "protein"
    assert format_term(hormone.subject.type) == "macromolecular_complex"
    assert format_term(secreted.subject.type) == "protein"


def test_cellchat_symbol_in_another_case_takes_the_gene_info_spelling():
    import pandas as pd
    from pypath.inputs_v2.parsers import cellchat as cellchat_parser

    db = {
        "interaction": pd.DataFrame(
            [
                {"ligand": "Cd99", "receptor": "Pilrb", "ligand.symbol": "Cd99", "receptor.symbol": "Pilrb"},
                {"ligand": "C1orf112", "receptor": "Abc", "ligand.symbol": "C1orf112", "receptor.symbol": "Abc"},
            ]
        ),
        "complex": pd.DataFrame(),
        # ABC and Abc differ only in case: a symbol matching both keeps its spelling.
        "geneInfo": pd.DataFrame({"Symbol": ["CD99", "PILRB", "C1orf112", "ABC", "Abc"]}),
    }
    db["interaction"].loc[1, "receptor.symbol"] = "aBC"
    rows = list(cellchat_parser._iter_interaction_rows(db, "9606"))
    assert [(r["ligand_genes"], r["receptor_genes"]) for r in rows] == [
        (["CD99"], ["PILRB"]),
        (["C1orf112"], ["aBC"]),
    ]


def test_cellchat_metabolite_abbreviation_is_the_longest_prefix():
    assert cellchat.metabolite("5-HT-DOPA-DDC") == "5-HT-DOPA"
    assert cellchat.metabolite("5-HT-TPH1_SLC6A4") == "5-HT"
    assert cellchat.metabolite("DHEAS-SULT1E1") == "DHEAS"
    assert cellchat.metabolite("DHEA-STS") == "DHEA"
    assert cellchat.metabolite("RARA") is None


def test_cellchat_metabolite_genes_are_no_complex_but_associated_once(monkeypatch):
    group = {"name": "5-HT-TPH1_SLC18A1", "genes": ["TPH1", "SLC18A1"], "taxon_id": "9606"}
    assert cellchat.map_cellchat_group({**group, "non_protein_ligand": True}) is None
    assert cellchat.map_cellchat_group({**group, "non_protein_ligand": False}) is not None

    rows = [
        _cellchat_row("5-HT-TPH1_SLC18A1", ["TPH1", "SLC18A1"]),
        _cellchat_row("5-HT-TPH1_SLC18A2", ["TPH1", "SLC18A2"]),
        _cellchat_row("CGA_TSHB", ["CGA", "TSHB"]),
    ]
    monkeypatch.setattr(cellchat, "iter_cellchat_interactions", lambda *a, **k: iter(rows))
    genes = list(cellchat.iter_metabolite_genes(None, taxon_id="9606"))
    assert [(g["gene"], g["metabolite"]) for g in genes] == [
        ("TPH1", "5-HT"), ("SLC18A1", "5-HT"), ("SLC18A2", "5-HT"),
    ]
    relation = cellchat.map_metabolite_gene(genes[0])
    assert format_term(relation.subject.type) == "protein"
    assert format_term(relation.object.type) == "chemical_entity"


def test_cellphonedb_complex_reference_and_definition_join(tmp_path):
    interaction = cellphonedb.interactions_schema({"partner_a": "complex A", "partner_b": "P00533"})
    definition = cellphonedb.complexes_schema({"complex_name": "complex A", "uniprot_1": "P04637"})
    entities = final([interaction, definition], tmp_path, "cellphonedb")["entity"]
    assert len([e for e in entities if e["entity_type"] == "macromolecular_complex"]) == 1


def test_cellphonedb_protein_written_by_its_symbol_is_a_protein(monkeypatch):
    tables = {
        "complexes": [{"complex_name": "IFNAR_complex"}],
        "proteins": [{"uniprot": "IFNA1"}, {"uniprot": "HLAA"}, {"uniprot": "P01562"}],
        "interactions": [{"partner_a": "IFNA1", "partner_b": "IFNAR_complex"}],
    }
    for name in ("complexes", "proteins"):
        monkeypatch.setattr(
            cellphonedb, f"download_{name}", SimpleNamespace(open=lambda name=name: name)
        )
    monkeypatch.setattr(cellphonedb, "iter_csv", lambda opener, **kw: iter(tables[opener]))
    row = next(cellphonedb.iter_interactions("interactions"))
    relation = cellphonedb.interactions_schema(row)
    assert format_term(relation.subject.type) == "protein"
    assert [(i.type, i.value) for i in relation.subject.identifiers] == [
        (Namespace.GENESYMBOL, "IFNA1")
    ]
    assert format_term(relation.object.type) == "macromolecular_complex"
    assert [(i.type, i.value) for i in relation.object.identifiers] == [
        (Namespace.NAME, "IFNAR_complex")
    ]
    hla = cellphonedb.interactions_schema({"partner_a": "HLAA", "partner_b": "P01562"})
    assert [(i.type, i.value) for i in hla.subject.identifiers] == [
        (Namespace.GENESYMBOL, "HLA-A")
    ]


def test_cellphonedb_metabolite_systems_carry_their_chebi_term():
    assert cellphonedb.metabolite_chebi("Cholesterol_byLIPA", "CHEBI:16113") == "CHEBI:16113"
    # Without a term in the table, the metabolite's name decides.
    assert cellphonedb.metabolite_chebi("Glutamate_byGLS_and_SLC1A1", "") == "CHEBI:16015"
    # Several terms, or a system of two metabolites, are no single metabolite.
    assert cellphonedb.metabolite_chebi(
        "Progesterone_byHSD3B1", "CHEBI:17026;CHEBI:16422;CHEBI:17252"
    ) is None
    assert cellphonedb.metabolite_chebi("SerotoninDopamin_byDDC", "CHEBI:28790") is None
    assert cellphonedb.metabolite_chebi("P00533", "") is None

    complex_ = cellphonedb.complexes_schema(
        {"complex_name": "Cholesterol_byLIPA", "secreted_desc": "CHEBI:16113"}
    )
    relation = cellphonedb.interactions_schema(
        {"partner_a": "Cholesterol_byLIPA", "partner_a_chebi": "CHEBI:16113", "partner_b": "P00533"}
    )
    for entity in (complex_, relation.subject):
        assert format_term(entity.type) == "chemical_entity"
        assert (Namespace.CHEBI, "CHEBI:16113") in {(i.type, i.value) for i in entity.identifiers}


def test_cellphonedb_roles_are_classes_below_their_role():
    protein = {
        "uniprot": "P01579", "secreted": "TRUE", "secreted_desc": "Growthfactor | Cytokine",
        "receptor": "FALSE", "transmembrane": "FALSE",
    }
    assert list(cellphonedb._roles(protein)) == [("secreted", "growthfactor"), ("secreted", "cytokine")]
    receptor = {"uniprot": "HLAA", "receptor": "TRUE", "receptor_desc": "", "transmembrane": "TRUE"}
    assert list(cellphonedb._roles(receptor)) == [("receptor", None), ("transmembrane", None)]
    relation = cellphonedb.map_protein_role({**receptor, "role": "receptor", "role_class": None})
    assert (relation.subject.identifiers[0].type, relation.subject.identifiers[0].value) == (
        Namespace.GENESYMBOL, "HLA-A",
    )
    assert relation.object.identifier == "receptor"
    # A class name can be a receptor class and a secreted class.
    receptor_class = cellphonedb._role_ref("receptor", "cytokine_receptor_il6_family")
    secreted_class = cellphonedb._role_ref("secreted", "cytokine_receptor_il6_family")
    assert receptor_class.identifier_type != secreted_class.identifier_type
    # The description of a metabolite system is its ChEBI term, not a class.
    system = {"complex_name": "Cholesterol_byLIPA", "secreted": "TRUE", "secreted_desc": "CHEBI:16113"}
    assert list(cellphonedb._roles(system)) == [("secreted", None)]
    assert cellphonedb._role_classes("Cytokine_receptor_add IL6 family") == ["cytokine_receptor_il6_family"]
    assert cellphonedb._role_classes("nuclear_receptor;needs transporter") == ["nuclear_receptor"]


def _macrophage_row(cell_1, ids_1, type_1, cell_2, ids_2, type_2, interaction="Binding"):
    return {
        "Partner 1 Official Gene Symbol": cell_1, "Partner 1 Gene ID": ids_1,
        "Partner 1 Interactant Type": type_1,
        "Partner 2 Official Gene Symbol": cell_2, "Partner 2 Gene ID": ids_2,
        "Partner 2 Interactant Type": type_2,
        "Interaction Type": interaction, "Species": None, "NCBI-PubMed ID": "16846590",
    }


def test_macrophage_partners_are_read_as_legacy_pypath_read_them(monkeypatch):
    rows = [
        # A complex of four members with five Gene IDs: the members by symbol.
        _macrophage_row("ATM:CHUK:IKBKB:IKBKG", "472:1147:3551:4214:23085", "Complex",
                        "NFKBIA", "4792", "Protein", "Phosphorylation"),
        # A protein cell with a form in brackets and no Gene ID.
        _macrophage_row("NCOR2(SMRT)", None, "Protein", "NFKB1(p50):RELA(p65)", "N/A", "Complex"),
        # A homodimer is no pair of a gene with itself.
        _macrophage_row("ATM:ATM", "472:472", "Complex", "ATM", "472", "Protein"),
    ]
    from pypath.inputs_v2 import macrophage

    monkeypatch.setattr(macrophage, "iter_interactions", lambda *a, **k: iter(rows))
    expanded = list(macrophage.iter_member_interactions(None))
    pairs = [(r["Partner 1 Official Gene Symbol"], r["Partner 2 Official Gene Symbol"]) for r in expanded]
    assert pairs == [
        ("ATM", "NFKBIA"), ("CHUK", "NFKBIA"), ("IKBKB", "NFKBIA"), ("IKBKG", "NFKBIA"),
        ("NCOR2", "NFKB1"), ("NCOR2", "RELA"),
    ]
    assert expanded[0]["Partner 1 Gene ID"] is None  # five IDs for four symbols
    relation = macrophage.interactions_mapper(expanded[4])
    assert [(i.type, i.value) for i in relation.subject.identifiers] == [(Namespace.GENESYMBOL, "NCOR2")]
    assert ("original_object", "NFKB1(p50):RELA(p65)") in {
        (format_term(a.term).split(":")[-1], a.value) for a in relation.annotations
    }


def test_oreganno_cells_of_several_names_give_the_first_as_legacy_did():
    from pypath.inputs_v2 import oreganno

    assert oreganno._symbol("STAT1/STAT3") == "STAT1"
    assert oreganno._symbol("E4TF1, Fli-1 (RBF Motif)") == "E4TF1"
    assert oreganno._symbol("NF-kappa B") == "NFkappa"
    # A cell of one name keeps its dash, which legacy removed.
    assert oreganno._symbol("HLA-A") == "HLA-A"
    relation = oreganno.interactions_schema(
        {"tf_symbol": "STAT1/STAT3", "gene_symbol": "OPRM1", "species": "Homo sapiens"}
    )
    assert {(i.type, i.value) for i in relation.subject.identifiers} >= {
        (Namespace.GENESYMBOL, "STAT1"), (Namespace.NAME, "STAT1/STAT3"),
    }


def test_mirdeathdb_symbol_with_several_products_gives_each_as_legacy_did():
    from pypath.inputs_v2.parsers.mirdeathdb import each_product, mirbase_products

    products = mirbase_products("hsa-let-7a-5p,hsa-let-7b-5p", "MIMAT0000062,MIMAT0000063")
    rows = each_product(products)
    assert [(r["mirbase_mature"], r["mirbase_mature_name"]) for r in rows] == [
        ("MIMAT0000062", "hsa-let-7a-5p"), ("MIMAT0000063", "hsa-let-7b-5p"),
    ]
    assert rows[0]["mature_organism_codes"] == ["hsa"]
    single = mirbase_products("hsa-miR-21-5p", "MIMAT0000076")
    assert each_product(single) == [single]
    # A malformed accession beside a valid one: the valid one, with its name.
    truncated = each_product(mirbase_products("hsa-miR-181b-5p,hsa-miR-181b-3p", "MIMAT0000257,MIMAT00226921"))
    assert [(r["mirbase_mature"], r["mirbase_mature_name"]) for r in truncated] == [("MIMAT0000257", "hsa-miR-181b-5p")]


def test_cellphonedb_hla_partner_is_the_protein_of_its_gene():
    relation = cellphonedb.interactions_schema({"partner_a": "HLAA", "partner_b": "P43629"})
    assert format_term(relation.subject.type) == "protein"
    assert [(i.type, i.value) for i in relation.subject.identifiers] == [(Namespace.GENESYMBOL, "HLA-A")]


def test_hmdb_keeps_the_kegg_compound():
    from pypath.inputs_v2 import hmdb

    kegg = hmdb.f("kegg_id", extract="kegg_compound", transform="kegg_compound")
    assert kegg.extract({"kegg_id": "C01152"}) == ["C01152"]
    assert kegg.extract({"kegg_id": "c1152"}) == ["C01152"]


def test_guidetopharma_peptide_ligand_has_the_gene_symbol_of_the_interaction():
    from pypath.inputs_v2 import guidetopharma

    ligand = {"Ligand ID": "5081", "Name": "4-1BB ligand", "UniProt ID": "P41273", "Type": "Peptide"}
    with_symbol = guidetopharma._with_gene_symbol(ligand, "TNFSF9")
    entity = guidetopharma.ligands_schema.build(with_symbol)
    assert (Namespace.GENESYMBOL, "TNFSF9") in {(i.type, i.value) for i in entity.identifiers}
    assert guidetopharma._with_gene_symbol(ligand, "Calca|Calca")["Gene Symbol"] == "Calca"
    assert "Gene Symbol" not in guidetopharma._with_gene_symbol(ligand, "INHBA|INHBB")


def test_chembl_units_are_selected_from_source_and_reach_quantity():
    con = duckdb.connect()
    try:
        con.execute("CREATE SCHEMA s")
        con.execute(
            "CREATE TABLE s.activities(activity_id BIGINT, assay_id BIGINT, molregno BIGINT, doc_id BIGINT, standard_type VARCHAR, standard_relation VARCHAR, standard_value DOUBLE, standard_units VARCHAR, pchembl_value DOUBLE, data_validity_comment VARCHAR, action_type VARCHAR)"
        )
        con.execute(
            "INSERT INTO s.activities VALUES (1, 1, 1, 1, 'IC50', '=', 25, 'uM', 5.6, NULL, NULL)"
        )
        query = con.execute(chembl_parser.DUCKDB_TABLES["activities"])
        row = dict(zip([d[0] for d in query.description], query.fetchone()))
        assert row["standard_units"] == "uM"
        # Both output backends must select this staging field as well.
        assert "act.standard_units" in chembl_parser.PARQUET_QUERIES["activities"]
        row.update(
            molecule_chembl_id="CHEMBL1", target_chembl_id="CHEMBL2", target_type="SINGLE PROTEIN"
        )
        record = chembl.activities_schema(row)
        quantity = next(
            quantity_dict(a.value)
            for a in record.annotations
            if format_term(a.term) == "BAO:0000190"
        )
        assert quantity["has_unit"] == "uM"
        assert quantity["has_numeric_value"] == 25
        assert observations([record]).relations[-1].upstream_id == "1"
    finally:
        con.close()


def test_swisslipids_multiple_crossrefs_and_formula():
    record = swisslipids.lipids_schema(
        {
            "Lipid ID": "SLM:000000784",
            "CHEBI": "74546 | 82922",
            "LIPID MAPS": "LM1 | LM2",
            "Formula (pH7.3)": "C2H4",
        }
    )
    pairs = {(str(i.type), i.value) for i in record.identifiers}
    assert {
        ("chebi", "74546"),
        ("chebi", "82922"),
        ("lipidmaps", "LM1"),
        ("lipidmaps", "LM2"),
    } <= pairs
    assert ("has_chemical_formula", "C2H4") in {
        (format_term(a.term), a.value) for a in record.annotations
    }


def test_ramp_pathway_name_and_source_domain():
    record = rampdb.pathway_schema(
        {
            "pathwayRampId": "RAMP_P_000000001",
            "pathwayName": "1-Methylhistidine Metabolism",
            "sourceId": "SMP0124716",
            "type": "hmdb",
        }
    )
    assert ("name", "1-Methylhistidine Metabolism") in {
        (str(i.type), i.value) for i in record.identifiers
    }
    assert ("xref", "SMPDB:SMP0124716") in {
        (format_term(a.term), a.value) for a in record.annotations
    }


def test_ontology_imports_keep_their_own_namespace():
    term = OntologyTerm(id="MONDO:0000001", name="example", is_a=["BFO:0000001"])
    record = ontology_term_to_entity(term, ontology_id="mondo", identifier_type=Namespace.MONDO)
    assert str(record.ontology_relations[0].object.identifier_type) == "bfo"
    imported = ontology_term_to_entity(
        term._replace(id="BFO:0000001", is_a=[]),
        ontology_id="mondo",
        identifier_type=Namespace.MONDO,
    )
    assert str(imported.identifiers[0].type) == "bfo"


def test_ensembl_protein_and_transcript_lookup_from_built_library(tmp_path):
    # ENSP/ENST/RefSeq-protein aliases of P04637 are rows of the shared fixture
    # hubs, so this reads them through the production-built shared reference.
    library = build_fixture_library(tmp_path / "reference")
    resolver = EntityResolver(library_dir=library)
    try:
        from omnipath_core.silver_schema import Entity, Identifier
        from biolink_model.datamodel.model import Protein

        records = [
            Entity(type=Protein, identifiers=[Identifier(type=ns, value=ident)])
            for ns, ident in [
                ("ensp", "ENSP00000269305"),
                ("ensembl", "ENSP00000269305.4"),
                ("enst", "ENST00000269305"),
                ("refseq_protein", "NP_000537.3"),
            ]
        ]
        out = observations(records)
        resolved = resolver.resolve_entities(out.entities, progress=False)
        assert len(resolved) == 4
        assert all(
            r.matched and (r.canonical_namespace, r.canonical_identifier) == ("entrez", "7157")
            for r in resolved.values()
        )
        assert all(r.entity_type == "protein" for r in out.entities.values())
        by_namespace = {out.entities[key].namespace: value for key, value in resolved.items()}
        for namespace in ("ensp", "ensembl", "refseq_protein"):
            product = by_namespace[namespace]
            assert (product.protein_namespace, product.protein_identifier) == ("uniprot", "P04637")
            assert product.protein_gene_candidates == ("entrez:7157",)
        transcript = by_namespace["enst"]
        assert transcript.protein_identifier is None
        assert (transcript.transcript_namespace, transcript.transcript_identifier) == (
            "enst",
            "ENST00000269305",
        )
    finally:
        resolver.close()


@pytest.mark.parametrize(
    "ns,value", [("ensg", "ENSG000001.3"), ("ensp", "ENSP000001.2"), ("enst", "ENSMUST000001.4")]
)
def test_ensembl_normalization_matches_library_sql(ns, value):
    con = duckdb.connect()
    try:
        normalized = con.execute(
            "SELECT " + normalize_id_sql("ns", "id") + " FROM (SELECT ? AS ns, ? AS id)",
            [ns, value],
        ).fetchone()[0]
        assert normalize_identifier(ns, value) == [(ns, normalized)]
    finally:
        con.close()


def _participant_values(tables, term):
    """{(endpoint label, value)} of a participant-scoped relation annotation."""
    relations = {r["relation_id"]: r for r in tables["relation"]}
    return {
        (relations[a["relation_id"]][f"{a['scope']}_label"], a["value"])
        for a in tables["relation_annotation"]
        if a["term"] == term and a["scope"] in ("subject", "object")
    }


def _relation_values(tables, term):
    return {
        a["value"]
        for a in tables["relation_annotation"]
        if a["term"] == term and a["scope"] == "relation"
    }


def test_connectomedb_locations_stay_with_their_participant(tmp_path):
    row = {
        "Interaction ID": "CDB15:0000001",
        "LR Pair": "A2M LRP1",
        "Ligand Symbols": "A2M (FWP007)",
        "Receptor Symbols": "LRP1 (APR, CD91)",
        "Ligand Species ID": "HGNC:7",
        "Receptor Species ID": "HGNC:6692",
        "Ligand Location": "secreted, cell membrane",
        "Receptor Location": "cell membrane",
        "Species": "Human",
    }
    tables = final([connectomedb.interactions_schema(row)], tmp_path, "connectomedb")
    assert _participant_values(tables, connectomedb.LOCATION) == {
        ("A2M", "secreted, cell membrane"),
        ("LRP1", "cell membrane"),
    }
    assert _participant_values(tables, "connectomedb:participant_role") == {
        ("A2M", "ligand"),
        ("LRP1", "receptor"),
    }


def test_cellchat_interaction_class_and_partner_attributes(tmp_path):
    protein = {
        **_cellchat_row("TGFB1", ["TGFB1"], annotation="Secreted Signaling"),
        "interaction_name": "TGFB1_TGFBR1_TGFBR2",
        "receptor_name": "TGFbR1_R2",
        "receptor_genes": ["TGFBR1", "TGFBR2"],
        "ligand_location": "Extracellular",
        "ligand_transmembrane": "FALSE",
        "ligand_secreted_type": "growth factor",
        "receptor_location": "Cell membrane",
        "receptor_transmembrane": "TRUE",
        "receptor_secreted_type": "",
    }
    metabolite = {
        **_cellchat_row("5-HT-TPH1_SLC18A1", ["TPH1", "SLC18A1"]),
        "receptor_location": "Cell membrane",
    }
    relations = [cellchat.map_cellchat_interaction(r) for r in (protein, metabolite)]
    tables = final(relations, tmp_path, "cellchat")
    assert _relation_values(tables, cellchat.INTERACTION_CLASS) == {
        "Secreted Signaling",
        "Non-protein Signaling",
    }
    roles = _participant_values(tables, "connectomedb:participant_role")
    assert ("TGFB1", "ligand") in roles and ("HTR1A", "receptor") in roles
    assert {value for _, value in roles} == {"ligand", "receptor"} and len(roles) == 4
    assert ("TGFB1", "growth factor") in _participant_values(tables, cellchat.SECRETED_TYPE)
    # An empty cell gives no annotation.
    assert {v for _, v in _participant_values(tables, cellchat.SECRETED_TYPE)} == {"growth factor"}
    assert {v for _, v in _participant_values(tables, cellchat.TRANSMEMBRANE)} == {"FALSE", "TRUE"}
    assert {v for _, v in _participant_values(tables, cellchat.LOCATION)} == {
        "Extracellular",
        "Cell membrane",
    }


def test_cellphonedb_interaction_attributes(tmp_path):
    row = {
        "partner_a": "P01137",
        "partner_b": "P36897",
        "directionality": "Ligand-Receptor",
        "classification": "Signaling by Transforming growth factor",
        "is_ppi": "True",
    }
    tables = final([cellphonedb.interactions_schema(row)], tmp_path, "cellphonedb")
    assert _relation_values(tables, cellphonedb.DIRECTIONALITY) == {"Ligand-Receptor"}
    assert _relation_values(tables, cellphonedb.CLASSIFICATION) == {
        "Signaling by Transforming growth factor"
    }
    assert _relation_values(tables, cellphonedb.IS_PPI) == {"True"}


def _icellnet_row(ligands, receptors):
    row = {f"Ligand {i}": "" for i in range(1, 5)} | {f"Receptor {i}": "" for i in range(1, 6)}
    row |= {f"Ligand {i}": gene for i, gene in enumerate(ligands, 1)}
    row |= {f"Receptor {i}": gene for i, gene in enumerate(receptors, 1)}
    return row | {"Reference": "1335692"}


def test_icellnet_spreadsheet_and_misspelt_cells_are_gene_symbols():
    assert icellnet._clean_gene("1,00\xa0DKK") == "DKK1"
    assert icellnet._clean_gene("3,00\xa0NPR") == "NPR3"
    assert icellnet._clean_gene("SCD2") == "SDC2"
    assert icellnet._clean_gene("NRP") == "NDP"


def test_icellnet_roles_follow_the_ligand_and_receptor_columns(tmp_path):
    pairs = [("OSTN", "3,00\xa0NPR"), ("HGF", "SDC1"), ("APOE", "SCARB1"), ("SDC2", "FGFR4")]
    relations = [icellnet.map_icellnet_interaction(_icellnet_row([l], [r])) for l, r in pairs]
    tables = final(relations, tmp_path, "icellnet")
    expected = {(icellnet._clean_gene(l), "ligand") for l, _ in pairs}
    expected |= {(icellnet._clean_gene(r), "receptor") for _, r in pairs}
    assert _participant_values(tables, "connectomedb:participant_role") == expected
    # The predicate is symmetric: some pairs are stored receptor first, and the
    # roles must still name the right endpoint.
    ligands = {l for l, _ in pairs}
    assert any(r["subject_label"] not in ligands for r in tables["relation"])
