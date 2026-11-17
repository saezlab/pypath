"""What the MetSigDB Parquet product publishes, set by set, over small purpose-made releases."""

import json

import duckdb
import pytest

from omnipath_subsets.parquet_products import build_release_subsets
from parquet_release import entity_row, relation_row, write_resources

CHEMICAL = "small_molecule"


def chemical(key, identifier="CHEBI:1", identifiers=(), label=None, group=None):
    return entity_row(
        key,
        CHEMICAL,
        "chebi",
        identifier,
        label or key,
        [("chebi", identifier, True), *identifiers],
        group,
    )


def term(
    key, identifier, kind="ontology_class", taxon=None, label=None, names=(), namespace="omnipath"
):
    row = entity_row(
        key,
        kind,
        namespace,
        identifier,
        label,
        [(namespace, identifier, True)] + [("name", name, False) for name in names],
    )
    row["taxon"] = taxon
    return row


def membership(source, number, subject, predicate, obj, subject_type, object_type):
    return relation_row(
        f"{source}:{number}",
        subject,
        predicate,
        obj,
        kind="relation",
        source=source,
        subject_type=subject_type,
        object_type=object_type,
        row=number,
    )


def ontology_edge(source, number, child, parent):
    return relation_row(
        f"{source}:o{number}",
        child,
        "subclass_of",
        parent,
        kind="ontology",
        source=source,
        subject_type="ontology_class",
    )


def rows_of(root, version, where="true", columns="*"):
    path = root / "subsets" / version / "metsigdb" / "metsigdb_membership.parquet"
    con = duckdb.connect()
    cursor = con.execute(f"SELECT {columns} FROM '{path}' WHERE {where} ORDER BY 1, 2, 3")
    names = [column[0] for column in cursor.description]
    return [dict(zip(names, row)) for row in cursor.fetchall()]


def build(root, release):
    build_release_subsets(root, release)
    return release


# ---------------------------------------------------------------------------
# ClassyFire sets are ChemOnt terms only


@pytest.fixture
def classyfire_release(tmp_path):
    cpd = chemical("cpd")
    other = chemical("other", "CHEBI:2")
    drug_class = term("class", "CHEMONTID:0000001", label="Organic compounds")
    disease = term("disease", "OMIM:123456", label="Obesity")
    version = write_resources(
        tmp_path,
        {
            "hmdb": (
                [cpd, other, drug_class, disease],
                [
                    membership(
                        "hmdb", 1, "cpd", "associated_with", "class", CHEMICAL, "ontology_class"
                    ),
                    membership(
                        "hmdb", 2, "other", "associated_with", "disease", CHEMICAL, "ontology_class"
                    ),
                ],
            ),
            "chemont": ([drug_class], []),
            "chebi": ([cpd, other], []),
        },
    )
    return tmp_path, build(tmp_path, version)


def test_classyfire_sets_are_chemont_terms_only(classyfire_release):
    root, version = classyfire_release
    sets = {r["set_source_id"] for r in rows_of(root, version, "resource = 'ClassyFire'")}
    assert sets == {"CHEMONTID:0000001"}


def test_a_class_that_chemont_does_not_publish_is_not_a_classyfire_set(tmp_path):
    cpd = chemical("cpd")
    not_chemont = term("tagged", "CHEMONTID:0000099", label="Not in the ontology")
    chemont_class = term("real", "CHEMONT-REAL", label="Real class")
    version = write_resources(
        tmp_path,
        {
            "hmdb": (
                [cpd, not_chemont, chemont_class],
                [
                    membership(
                        "hmdb", 1, "cpd", "associated_with", "tagged", CHEMICAL, "ontology_class"
                    ),
                    membership(
                        "hmdb", 2, "cpd", "associated_with", "real", CHEMICAL, "ontology_class"
                    ),
                ],
            ),
            "chemont": ([chemont_class], []),
        },
    )
    build(tmp_path, version)
    sets = {r["set_source_id"] for r in rows_of(tmp_path, version, "resource = 'ClassyFire'")}
    # ChemOnt's own list of terms decides, whatever the accession looks like.
    assert sets == {"CHEMONT-REAL"}


# ---------------------------------------------------------------------------
# KEGG reference maps have no organism; organism maps are overview maps by number


def kegg_release(root, maps):
    """A KEGG release: each map is (identifier, taxon); every map has one reaction and compound."""
    cpd = chemical("cpd")
    entities, relations = [cpd], []
    for number, (identifier, taxon) in enumerate(maps, start=1):
        entities.append(term(f"map{number}", identifier, kind="pathway", taxon=taxon))
        entities.append(term(f"rx{number}", f"R{number:05d}", kind="molecular_activity"))
        relations.append(
            membership(
                "kegg",
                number,
                f"map{number}",
                "associated_with",
                f"rx{number}",
                "pathway",
                "molecular_activity",
            )
        )
        relations.append(
            membership(
                "kegg",
                100 + number,
                f"rx{number}",
                "has_input",
                "cpd",
                "molecular_activity",
                CHEMICAL,
            )
        )
    version = write_resources(root, {"kegg": (entities, relations), "chebi": ([cpd], [])})
    return build(root, version)


def test_a_kegg_reference_map_has_no_organism(tmp_path):
    version = kegg_release(
        tmp_path, [("rn00010", "10090"), ("hsa00010", "9606"), ("mmu00010", "10090")]
    )
    organisms = {r["set_source_id"]: r["organism"] for r in rows_of(tmp_path, version)}
    assert organisms == {"rn00010": None, "hsa00010": 9606, "mmu00010": 10090}


def test_an_organism_specific_overview_map_is_an_overview_map(tmp_path):
    version = kegg_release(
        tmp_path,
        [
            ("rn01100", None),
            ("hsa01100", "9606"),
            ("mmu01240", "10090"),
            ("hsa00010", "9606"),
            ("rn00010", None),
            ("hsa11000", "9606"),
        ],
    )
    kinds = {r["set_source_id"]: r["set_sub_type"] for r in rows_of(tmp_path, version)}
    assert kinds == {
        "rn01100": "overview_map",
        "hsa01100": "overview_map",
        "mmu01240": "overview_map",
        "hsa00010": "metabolic_map",
        "rn00010": "metabolic_map",
        # Five digits that only end alike are not a map number.
        "hsa11000": "metabolic_map",
    }


# ---------------------------------------------------------------------------
# Identifier columns


@pytest.mark.parametrize(
    "held, shown",
    [
        ("HMDB0000008", "HMDB0000008"),
        ("HMDB00008", "HMDB0000008"),
        ("HMDB000008", "HMDB0000008"),
        ("hmdb0000008", "HMDB0000008"),
        ("HMDB:0000008", "HMDB0000008"),
        (" HMDB0000008 ", "HMDB0000008"),
        ("HMDB00000008", "HMDB0000008"),
    ],
)
def test_every_spelling_of_an_hmdb_accession_is_shown_with_seven_digits(tmp_path, held, shown):
    cpd = chemical("cpd", identifiers=[("hmdb", held, False)])
    pathway = term("wp", "WP1", kind="pathway", label="A pathway")
    version = write_resources(
        tmp_path,
        {
            "wikipathways": (
                [cpd, pathway],
                [membership("wikipathways", 1, "wp", "has_part", "cpd", "pathway", CHEMICAL)],
            ),
            "chebi": ([cpd], []),
        },
    )
    build(tmp_path, version)
    assert [r["hmdb"] for r in rows_of(tmp_path, version, "resource = 'WikiPathways'")] == [shown]


def test_an_identifier_stated_by_another_resource_fills_the_column(tmp_path):
    cpd = chemical("cpd")
    from_hmdb = chemical("cpd", identifiers=[("hmdb", "HMDB0000008", False)])
    from_kegg = chemical("cpd", identifiers=[("kegg", "C00008", False), ("pubchem", "8", False)])
    pathway = term("wp", "WP1", kind="pathway", label="A pathway")
    version = write_resources(
        tmp_path,
        {
            "wikipathways": (
                [cpd, pathway],
                [membership("wikipathways", 1, "wp", "has_part", "cpd", "pathway", CHEMICAL)],
            ),
            "hmdb": ([from_hmdb], []),
            "kegg": ([from_kegg], []),
        },
    )
    build(tmp_path, version)
    (row,) = rows_of(tmp_path, version, "resource = 'WikiPathways'")
    assert (row["hmdb"], row["kegg"], row["pubchem"], row["chebi"]) == (
        "HMDB0000008",
        "C00008",
        "8",
        "1",
    )


# ---------------------------------------------------------------------------
# Exclusion classes, structure presence, MACdb ids


def release_with_classes(root, seeds):
    cofactor = chemical(
        "cofactor", "CHEBI:57540", identifiers=[("inchikey", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", False)]
    )
    plain = chemical("plain", "CHEBI:7")
    ion_root = term("ions", "CHEBI:24867", kind=CHEMICAL, namespace="chebi", label="monoatomic ion")
    pathway = term("wp", "WP1", kind="pathway", label="A pathway")
    version = write_resources(
        root,
        {
            "wikipathways": (
                [cofactor, plain, pathway],
                [
                    membership(
                        "wikipathways", 1, "wp", "has_part", "cofactor", "pathway", CHEMICAL
                    ),
                    membership("wikipathways", 2, "wp", "has_part", "plain", "pathway", CHEMICAL),
                ],
            ),
            "chebi": (
                [cofactor, plain, ion_root],
                [
                    relation_row(
                        "chebi:o1",
                        "cofactor",
                        "subclass_of",
                        "ions",
                        kind="ontology",
                        source="chebi",
                    ),
                ],
            ),
        },
    )
    build_release_subsets(root, version, seed_directory=seeds)
    return version


@pytest.fixture
def classed(tmp_path):
    from test_exclusion_metabolites import seed_directory

    seeds = tmp_path / "seeds"
    seeds.mkdir()
    seed_directory(
        seeds,
        metaproviz="seeds: []\n",
        overrides="add: []\nveto: []\n",
    )
    return release_with_classes(tmp_path, seeds), tmp_path


def test_each_row_carries_the_exclusion_classes_of_its_metabolite(classed):
    version, root = classed
    classes = {r["metabolite_entity_key"]: r["exclusion_classes"] for r in rows_of(root, version)}
    assert classes == {"cofactor": ["cofactors", "ions"], "plain": []}


def test_a_row_says_whether_its_metabolite_has_a_structure(classed):
    version, root = classed
    structure = {r["metabolite_entity_key"]: r["has_structure"] for r in rows_of(root, version)}
    assert structure == {"cofactor": True, "plain": False}


def test_a_structure_key_alone_is_a_structure(tmp_path):
    cpd = chemical("cpd", group="CCCCCCCCCCCCCC")
    pathway = term("wp", "WP1", kind="pathway", label="A pathway")
    version = write_resources(
        tmp_path,
        {
            "wikipathways": (
                [cpd, pathway],
                [membership("wikipathways", 1, "wp", "has_part", "cpd", "pathway", CHEMICAL)],
            )
        },
    )
    build(tmp_path, version)
    assert [r["has_structure"] for r in rows_of(tmp_path, version)] == [True]


def test_macdb_set_ids_are_prefixed_and_do_not_collide_with_chebi(tmp_path):
    cpd = chemical("cpd", "CHEBI:15378")
    trait = term("trait", "15378", label="A trait")
    version = write_resources(
        tmp_path,
        {
            "macdb": (
                [cpd, trait],
                [
                    membership(
                        "macdb", 1, "cpd", "associated_with", "trait", CHEMICAL, "ontology_class"
                    )
                ],
            )
        },
    )
    build(tmp_path, version)
    assert [r["set_source_id"] for r in rows_of(tmp_path, version, "resource = 'MACdb'")] == [
        "MACDB:15378"
    ]


def test_the_manifest_names_the_row_contract_of_the_product(tmp_path):
    version = release_with_classes_manifest(tmp_path)
    manifest = json.loads(
        (tmp_path / "subsets" / version / "metsigdb" / "subset_manifest.json").read_text()
    )
    assert manifest["row_contract"] == 2
    assert manifest["stats"]["rows"] == 1


def release_with_classes_manifest(root):
    cpd = chemical("cpd")
    pathway = term("wp", "WP1", kind="pathway", label="A pathway")
    version = write_resources(
        root,
        {
            "wikipathways": (
                [cpd, pathway],
                [membership("wikipathways", 1, "wp", "has_part", "cpd", "pathway", CHEMICAL)],
            )
        },
    )
    return build(root, version)


# ---------------------------------------------------------------------------
# Coverage, measured at build time


def test_the_build_measures_what_share_of_each_resource_carries_each_identifier(tmp_path):
    full = chemical(
        "full",
        identifiers=[
            ("hmdb", "HMDB0000008", False),
            ("kegg", "C00008", False),
            ("pubchem", "8", False),
            ("inchikey", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", False),
        ],
    )
    bare = chemical("bare", "CHEBI:2")
    pathway = term("wp", "WP1", kind="pathway", label="A pathway")
    version = write_resources(
        tmp_path,
        {
            "wikipathways": (
                [full, bare, pathway],
                [
                    membership("wikipathways", 1, "wp", "has_part", "full", "pathway", CHEMICAL),
                    membership("wikipathways", 2, "wp", "has_part", "bare", "pathway", CHEMICAL),
                ],
            )
        },
    )
    build(tmp_path, version)
    manifest = json.loads(
        (tmp_path / "subsets" / version / "metsigdb" / "subset_manifest.json").read_text()
    )
    coverage = manifest["stats"]["coverage"]["WikiPathways"]
    assert coverage["metabolites"] == 2 and coverage["sets"] == 1 and coverage["rows"] == 2
    assert coverage["identifiers"] == {
        "hmdb": 0.5,
        "chebi": 1.0,
        "kegg": 0.5,
        "pubchem": 0.5,
        "inchikey": 0.5,
    }
