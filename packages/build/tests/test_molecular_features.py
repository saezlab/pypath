"""Source sequence features become structured molecular forms, not raw JSON."""

import csv
import io
from types import SimpleNamespace

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF

from pypath.inputs_v2 import mirbase, uniprot
from pypath.inputs_v2.parsers.reactome import BP, _participant_molecular_form

PTMLIST = """ID   Phosphoserine
AC   PTM-0253
FT   MOD_RES
DR   PSI-MOD; MOD:00046.
//
ID   N-linked (GlcNAc...) asparagine
AC   PTM-0526
FT   CARBOHYD
//
ID   Glycyl lysine isopeptide (Lys-Gly) (interchain with G-...)
AC   PTM-0640
FT   CROSSLNK
DR   PSI-MOD; MOD:00134.
//
"""
SEQUENCE = "MNCSKCSTAAAANAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"


def opener(text):
    return SimpleNamespace(result={"data": io.StringIO(text)})


def catalogue(**features):
    data = io.StringIO()
    row = {
        "Entry": "P04637",
        "Organism (ID)": "9606",
        "Sequence": SEQUENCE,
        "Sequence version": "4",
        **features,
    }
    writer = csv.DictWriter(data, fieldnames=list(row), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
    rows = uniprot._catalogue_feature_rows(opener(data.getvalue()), ptmlist=opener(PTMLIST))
    return [uniprot.catalogue_features_schema(row) for row in rows]


def test_uniprot_modifications_use_psi_mod_and_keep_source_notes():
    records = catalogue(
        **{
            "Modified residue": 'MOD_RES 4; /note="Phosphoserine; by PKA"; '
            '/evidence="ECO:0000269|PubMed:123"; MOD_RES P04637-2:7; /note="Phosphoserine"',
            "Glycosylation": 'CARBOHYD 13; /note="N-linked (GlcNAc...) asparagine"',
            "Cross-link": 'CROSSLNK 5; /note="Glycyl lysine isopeptide (Lys-Gly) '
            '(interchain with G-Cter in SUMO2); alternate"',
            "Disulfide bond": 'DISULFID 3..6; /evidence="ECO:0000250"; DISULFID ?..6; '
            '/note="Interchain"',
        }
    )
    assert all(record.annotations[0].term == "in_taxon" for record in records)
    assert not any(
        "catalogue_feature" in annotation.term
        for record in records
        for annotation in record.annotations
    )
    assert [a.value for a in records[0].annotations if a.term == "publications"] == ["PMID:123"]
    forms = [record.molecular_form for record in records]
    (canonical,) = forms[0]["modifications"]
    assert (canonical["term"], canonical["residue"], canonical["position"]) == ("MOD:00046", "S", 4)
    assert canonical["description"] == "Phosphoserine; by PKA"
    assert canonical["coordinate_reference"]["identifier"]["ns"] == "protein_sequence_sha256"
    # Isoform coordinates keep the isoform reference; its residue is not looked up.
    assert forms[1]["isoform_identifier"] == {"ns": "uniprot", "id": "P04637-2"}
    (isoform,) = forms[1]["modifications"]
    assert isoform["coordinate_reference"]["identifier"] == {"ns": "uniprot", "id": "P04637-2"}
    assert (isoform["position"], isoform["residue"]) == (7, None)
    # A vocabulary entry without a PSI-MOD cross-reference keeps its UniProt name.
    assert forms[2]["modifications"][0]["term"] == "N-linked (GlcNAc...) asparagine"
    # A disulfide 'a..b' links two cysteines; it is not a modified range.
    bond = forms[3]["modifications"]
    assert [(m["term"], m["residue"], m["position"], m["end_position"]) for m in bond] == [
        ("Disulfide bond", "C", 3, 3),
        ("Disulfide bond", "C", 6, 6),
    ]
    uncertain = forms[4]["modifications"]
    assert [m["position"] for m in uncertain] == [None, 6]
    assert uncertain[0]["description"] == "Interchain; source location ?..6"
    # Interchain crosslink names match UniProt's partner-residue template.
    assert forms[5]["modifications"][0]["term"] == "MOD:00134"


def test_uniprot_variants_keep_ids_alleles_and_deletions():
    records = catalogue(
        **{
            "Natural variant": 'VARIANT 3; /note="C -> R (in dbSNP:rs1)"; /id="VAR_000001"',
            "Alternative sequence": 'VAR_SEQ 1..4; /note="Missing (in isoform 2)"; '
            '/id="VSP_060791"; VAR_SEQ 9..10; /note="AA -> MPR (in isoform 3)"; /id="VSP_2"',
        }
    )
    variants = [record.molecular_form["variants"][0] for record in records]
    assert [
        (v["identifier"]["id"], v["reference"], v["alternate"], v["position"], v["end_position"])
        for v in variants
    ] == [
        ("VAR_000001", "C", "R", 3, 3),
        ("VSP_060791", None, "", 1, 4),
        ("VSP_2", "AA", "MPR", 9, 10),
    ]
    assert variants[1]["description"] == "Missing (in isoform 2)"
    assert all(record.molecular_form["isoform_identifier"] is None for record in records)


def test_uniprot_ranged_features_become_regions_without_invented_ends():
    records = catalogue(
        **{
            "Chain": 'CHAIN 2..40; /note="Cellular tumor antigen p53"; /id="PRO_0000185703"; '
            'CHAIN ?..30; /note="Processed form"; /id="PRO_0000000001"',
            "Signal peptide": 'SIGNAL 1..>20; /evidence="ECO:0000255"',
            "Transmembrane": 'TRANSMEM 21..41; /note="Helical"',
        }
    )
    forms = [record.molecular_form for record in records]
    assert forms[0]["sequence_identifiers"] == [{"ns": "uniprot", "id": "P04637-PRO_0000185703"}]
    (chain,) = forms[0]["regions"]
    assert chain["type"] == "Chain"
    assert chain["identifier"] == {"ns": "uniprot_feature", "id": "PRO_0000185703"}
    assert (chain["position"], chain["end_position"]) == (2, 40)
    assert chain["description"] == "Cellular tumor antigen p53"
    assert chain["coordinate_reference"]["position_base"] == 1
    fuzzy, signal, helix = forms[1]["regions"][0], forms[2]["regions"][0], forms[3]["regions"][0]
    assert (fuzzy["position"], fuzzy["end_position"]) == (None, 30)
    assert (signal["type"], signal["position"], signal["end_position"]) == (
        "Signal peptide",
        1,
        None,
    )
    assert signal["description"] == "source location 1..>20"
    assert (helix["type"], helix["description"]) == ("Transmembrane", "Helical")
    assert all(not form["modifications"] and not form["variants"] for form in forms)


def test_reactome_fragments_become_regions_and_descriptions_are_readable():
    graph = Graph()
    protein = URIRef("urn:protein")
    sites = {}
    for name, position, status in [
        ("begin", 29, "EQUAL"),
        ("end", 175, "EQUAL"),
        ("site", 15, "EQUAL"),
        ("fuzzy", 20, "LESS-THAN"),
    ]:
        sites[name] = URIRef(f"urn:{name}")
        graph.add((sites[name], BP.sequencePosition, Literal(position)))
        graph.add((sites[name], BP.positionStatus, Literal(status)))
    fragment, interval = URIRef("urn:fragment"), URIRef("urn:interval")
    graph.add((protein, BP.feature, fragment))
    graph.add((fragment, RDF.type, BP.FragmentFeature))
    graph.add((fragment, BP.featureLocation, interval))
    graph.add((interval, BP.sequenceIntervalBegin, sites["begin"]))
    graph.add((interval, BP.sequenceIntervalEnd, sites["end"]))
    partial, partial_interval = URIRef("urn:partial"), URIRef("urn:partial-interval")
    graph.add((protein, BP.feature, partial))
    graph.add((partial, RDF.type, BP.FragmentFeature))
    graph.add((partial, BP.featureLocation, partial_interval))
    graph.add((partial_interval, BP.sequenceIntervalBegin, sites["begin"]))
    graph.add((partial_interval, BP.sequenceIntervalEnd, sites["fuzzy"]))
    modification, vocabulary, xref = URIRef("urn:mod"), URIRef("urn:vocabulary"), URIRef("urn:x")
    graph.add((protein, BP.feature, modification))
    graph.add((modification, BP.modificationType, vocabulary))
    graph.add((modification, BP.featureLocation, sites["site"]))
    graph.add((vocabulary, BP["term"], Literal("O-phospho-L-serine")))
    graph.add((vocabulary, BP.xref, xref))
    graph.add((xref, BP.db, Literal("MOD")))
    graph.add((xref, BP.id, Literal("MOD:00046")))
    form = _participant_molecular_form(
        graph, protein, {"entity_type": "protein", "uniprot": "P42574-1"}
    )
    (mod,) = form["modifications"]
    assert (mod["term"], mod["position"], mod["description"]) == (
        "MOD:00046",
        15,
        "O-phospho-L-serine",
    )
    regions = sorted(form["regions"], key=lambda region: region["description"] or "")
    assert [(r["type"], r["position"], r["end_position"]) for r in regions] == [
        ("Fragment", 29, 175),
        ("Fragment", 29, None),
    ]
    assert regions[0]["description"] is None
    assert regions[1]["description"] == "source position 29 (EQUAL)..20 (LESS-THAN)"
    assert regions[0]["coordinate_reference"]["identifier"] == {"ns": "uniprot", "id": "P42574-1"}


def test_mirbase_precursor_regions_are_mature_form_regions():
    entry = """ID   hsa-mir-1 standard; RNA;
AC   MI0000001;
FT   miRNA           3..6
FT                   /accession="MIMAT0000001"
FT                   /product="hsa-miR-1"
SQ   Sequence 8 BP;
     acguacgu 8
//
"""
    (row,) = mirbase._matures_raw(opener(entry))
    record = mirbase.matures_schema(row)
    assert not record.annotations
    form = record.molecular_form
    assert form["sequence_identifiers"][0]["ns"] == "transcript_sequence_sha256"
    (region,) = form["regions"]
    assert (region["type"], region["position"], region["end_position"]) == ("Mature miRNA", 3, 6)
    assert region["coordinate_reference"] == {
        "identifier": {"ns": "mirbase_precursor_release", "id": "MI0000001@22"},
        "coordinate_system": "transcript",
        "position_base": 1,
    }
