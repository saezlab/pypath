"""Resource forms survive real resolution, Parquet reduction and direct serving."""

import csv
import io
import json
from types import SimpleNamespace


from library_fixture import build_fixture_library
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from tables_fixture import read_tables
from omnipath_resolver.resolver import EntityResolver
from pypath.inputs_v2 import bindingdb, brenda, chembl, mirbase, uniprot
from pypath.inputs_v2.parsers import brenda as brenda_parser


def opener(text):
    return SimpleNamespace(result={"data": io.StringIO(text)})


def catalogue_rows():
    data = io.StringIO()
    row = {
        "Entry": "P04637",
        "Organism (ID)": "9606",
        "Sequence": "MASG",
        "Sequence version": "4",
        "Mutagenesis": 'MUTAGEN 2; /note="A->G,V: alternatives"',
        "Modified residue": 'MOD_RES 3; /note="Phosphoserine"',
    }
    writer = csv.DictWriter(data, fieldnames=list(row), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
    return [
        uniprot.catalogue_features_schema(r)
        for r in uniprot._catalogue_feature_rows(
            opener(data.getvalue()),
            ptmlist=opener("ID   Phosphoserine\nFT   MOD_RES\nDR   PSI-MOD; MOD:00046.\n//\n"),
        )
    ]


def test_molecular_inputs_roundtrip_without_combinatorial_entities(tmp_path):
    folder = tmp_path / "resources" / "fixture" / "1"
    rows = catalogue_rows()
    for row in brenda_parser.iter_molecular_forms(
        opener(
            "ID\t1.1.1.1\nPR\t#1# Homo sapiens {P04637; source: UniProt} <1>\n"
            "EN\t#1# A2G/S3D <1>\nEN\t#1# A2V <1>\n///\n"
        )
    ):
        rows.append(brenda.molecular_forms_schema(row))
    rows.append(
        chembl.activities_schema(
            {
                "molecule_chembl_id": "CHEMBL1",
                "target_chembl_id": "CHEMBLT1",
                "target_type": "SINGLE PROTEIN",
                "target_component_uniprot_accessions": "P04637",
                "target_tax_id": "9606",
                "variant_id": 12,
                "variant_accession": "P04637",
                "variant_isoform": 2,
                "variant_version": 4,
                "variant_mutation": "A2G/S3D",
                "variant_sequence": "MGDG",
            }
        )
    )
    rows.append(
        bindingdb._target(
            {
                "Target Name": "two source chains",
                "Number of Protein Chains in Target (>1 implies a multichain complex)": "2",
                "UniProt (SwissProt) Primary ID of Target Chain 1": "P04637-2",
                "UniProt (SwissProt) Primary ID of Target Chain 2": "P04637-3",
                "BindingDB Target Chain 1 Sequence": "MAAK",
                "BindingDB Target Chain 2 Sequence": "MAAT",
            }
        )
    )
    rows.append(
        mirbase.matures_schema(
            {
                "mirbase_mat": "MIMAT1",
                "name": "a mature RNA",
                "sequence": "ACGU",
                "precursors": ["MI1"],
                "precursor_regions": [
                    {
                        "location": "2..5",
                        "position": 2,
                        "end_position": 5,
                        "coordinate_reference": {
                            "identifier": {"ns": "mirbase_precursor_release", "id": "MI1@22"},
                            "coordinate_system": "transcript",
                            "position_base": 1,
                        },
                    }
                ],
            }
        )
    )
    resolver = EntityResolver(build_fixture_library(tmp_path))
    writer = ParquetWriter(folder, library_dir=resolver.library_dir)
    try:
        for number, record in enumerate(rows):
            extractor = SilverExtractor("fixture", "molecular")
            extractor.process_record(record, {"source_row": number}, str(number), number)
            writer.append_observations(extractor, resolver)
        tables = read_tables(writer.close()["files"])
    finally:
        resolver.close()
    entities = tables["entity"]
    (anchor,) = [
        entity
        for entity in entities
        if entity["namespace"] == "entrez" and entity["entity_type"] == "protein"
    ]
    assert anchor["reference_entity_key"] == "entrez:7157"
    assert not any(
        identifier["ns"].endswith("_sequence_sha256")
        for identifier in tables.children("entity_identifier", anchor)
    )
    (product,) = [
        entity
        for entity in entities
        if entity["namespace"] == "uniprot" and entity["identifier"] == "P04637"
    ]
    assert product["reference_entity_key"] == anchor["reference_entity_key"]
    assert not any(
        "-" in identifier["id"]
        for identifier in tables.children("entity_identifier", anchor)
        if identifier["ns"] == "uniprot"
    )
    alternatives = [
        evidence["molecular_form"]
        for evidence in tables.children("entity_evidence", anchor)
        if evidence["molecular_form"]
    ]
    assert len(alternatives) == 7
    assert sorted(len(form["variants"] or []) for form in alternatives) == [0, 0, 0, 1, 1, 1, 2]
    assert all(form["protein_entity_key"] == product["entity_key"] for form in alternatives)
    assert any(
        modification["term"] == "MOD:00046"
        for form in alternatives
        for modification in form["modifications"] or []
    )
    (mature,) = [entity for entity in entities if entity["namespace"] == "mirbase_mature"]
    (region,) = tables.children("entity_evidence", mature)[0]["molecular_form"]["regions"]
    assert (region["type"], region["position"], region["end_position"]) == ("Mature miRNA", 2, 5)
    mutant = next(r for r in tables["relation"] if r["predicate"] == "interacts_with")
    form = tables.children("relation_evidence", mutant)[0]["object_molecular_form"]
    assert form["protein_entity_key"] == product["entity_key"]
    assert len(form["variants"]) == 2
    assert {identifier["ns"] for identifier in form["sequence_identifiers"]} >= {
        "uniprot",
        "chembl_representative_sequence_sha256",
    }
    engine = ParquetQueries(tmp_path)
    try:
        served = engine.get_entity_core(anchor["entity_key"])["entity"]
        assert served["molecularEvidenceTotal"] == 7
        evidence = engine.get_relation_evidence(mutant["relation_key"])["evidence"]
        assert len(evidence) == 1
        assert len(evidence[0]["objectMolecularForm"]["variants"]) == 2
        context = engine.get_molecular_context(
            product["entity_key"], view="product", isoform_identifier="uniprot:P04637-2"
        )
        assert context["relationsTotal"] == 2
        assert all(
            "uniprot:P04637-3" not in json.dumps(row["evidence"]) for row in context["relations"]
        )
    finally:
        engine.close()
