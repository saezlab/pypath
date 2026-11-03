"""Independent gene and product identity against explicit reference links."""

from contextlib import closing

import pytest

from library_fixture import build_fixture_library
from test_canonical import obs
from omnipath_resolver.canonical.match import LibraryMatcher
from omnipath_resolver._omnipath_resolver import resolve_molecular_batch


def decide(votes, mappings):
    facts = [
        (1, "uniprot:P04637", 2, "uniprot:P04637", False, True, ["7157"]),
        (2, "uniprot:A0A0U1RQF1", 2, "uniprot:A0A0U1RQF1", False, False, ["7157"]),
        (3, "entrez:7157", 3, None, False, False, ["7157"]),
        (4, "uniprot:P0DP23", 2, "uniprot:P0DP23", False, True, ["801", "805", "808"]),
        (5, "entrez:55", 3, None, False, False, ["55"]),
        (6, "uniprot:Q99999", 2, "uniprot:Q99999", False, False, []),
        (7, "entrez:805", 3, None, False, False, ["805"]),
    ]
    return resolve_molecular_batch([("input", 2, votes)], mappings, facts)[0]


def v(key, gene=False, route=1):
    return key, "", route, 0, gene, False


def test_gene_only_never_asserts_catalogue_products():
    for ids in ([1], [1, 2], [3]):
        result = decide([v(b"gene", True)], [(b"gene", ids)])
        assert result[2] == ["entrez:7157"]
        assert result[4] is None
        assert result[6] == "resolved"


def test_distinct_products_share_gene_without_merging_product_identity():
    for candidate, accession in ((1, "P04637"), (2, "A0A0U1RQF1")):
        result = decide([v(b"protein")], [(b"protein", [candidate])])
        assert result[2] == ["entrez:7157"]
        assert result[4] == "uniprot:" + accession


def test_product_with_ambiguous_gene_retains_primary_protein():
    result = decide([v(b"protein")], [(b"protein", [4])])
    assert result[2] == ["uniprot:P0DP23"]
    assert result[4] == "uniprot:P0DP23"
    assert result[5] == ["entrez:801", "entrez:805", "entrez:808"]
    assert result[6] == "ambiguous"


def test_gene_identifier_disambiguates_product_gene_links():
    result = decide([v(b"protein"), v(b"gene", True)], [(b"protein", [4]), (b"gene", [7])])
    assert result[2] == ["entrez:805"]
    assert result[4] == "uniprot:P0DP23"
    assert result[6] == "resolved"


def test_gene_product_mismatch_is_diagnosable_supported_protein_fallback():
    result = decide([v(b"protein"), v(b"gene", True)], [(b"protein", [1]), (b"gene", [5])])
    assert result[2] == ["uniprot:P04637"]
    assert result[4] == "uniprot:P04637"
    assert result[5] == ["entrez:55", "entrez:7157"]
    assert result[6] == "conflict"


def test_missing_gene_does_not_discard_supported_protein():
    result = decide([v(b"protein")], [(b"protein", [6])])
    assert result[2] == ["uniprot:Q99999"]
    assert result[4] == "uniprot:Q99999"
    assert result[6] == "missing"


@pytest.mark.parametrize("defer", [False, True])
def test_coordinated_result_preserves_source_type_and_product_enrichment(tmp_path, defer):
    with closing(LibraryMatcher(build_fixture_library(tmp_path), defer_aliases=defer)) as matcher:
        result = matcher.match(
            {
                "p": obs("protein", "uniprot", "P04637"),
                "u": obs("protein", "uniprot", "A0A0U1RQF1"),
                "g": obs("gene", "entrez", "7157"),
                "pg": obs("protein", "entrez", "7157"),
                "r": obs("rna_product", "uniprot", "P04637"),
            }
        )
        assert {m.canonical_identifier for m in result.values()} == {"7157"}
        assert [result[k].entity_type for k in ("p", "g", "pg", "r")] == [
            "protein",
            "gene",
            "protein",
            "rna_product",
        ]
        assert result["p"].protein_identifier == "P04637"
        assert result["u"].protein_identifier == "A0A0U1RQF1"
        assert result["p"].protein_gene_candidates == ("entrez:7157",)
        assert result["p"].protein_aliases["uniprot"] == ["P04637"]
        assert all(result[k].protein_identifier is None for k in ("g", "pg", "r"))


def test_missing_supplied_gene_still_detects_conflicting_known_product():
    result = decide(
        [v(b"protein"), (b"gene", "entrez:999999", 2, 1, True, False)],
        [(b"protein", [1]), (b"gene", [])],
    )
    assert result[2] == ["uniprot:P04637"]
    assert result[6] == "conflict"


def test_entrez_hub_retains_noncoding_transcript_links(monkeypatch):
    from omnipath_build.hubs.sources import entrez

    monkeypatch.setattr(
        entrez,
        "iter_url_lines",
        lambda url: iter(
            [
                "9606\t1\tENSG000001\tNR_001.2\tENST000001.3\t-\t-\n",
            ]
        ),
    )

    class Writer:
        full = False

        def __init__(self):
            self.rows = []

        def add_identity(self, *args):
            return True

        def add(self, namespace, identifier, *rest):
            self.rows.append((namespace, identifier))
            return True

    writer = Writer()
    entrez.emit(writer)
    assert ("enst", "ENST000001") in writer.rows
    assert ("refseq", "NR_001.2") in writer.rows


@pytest.mark.parametrize(
    "entity_type",
    [
        "rna_product",
        "rna_product_isoform",
        "noncoding_rna_product",
        "microrna",
        "sirna",
        "transcript",
    ],
)
def test_rna_subtypes_have_gene_policy_without_asserting_protein(tmp_path, entity_type):
    from omnipath_resolver.canonical.policy import get_policy, RNA_ENTITY_TYPES

    assert entity_type in RNA_ENTITY_TYPES
    assert get_policy(entity_type).library == "gene_protein"
    with closing(LibraryMatcher(build_fixture_library(tmp_path))) as matcher:
        result = matcher.match({"r": obs(entity_type, "uniprot", "P04637")})["r"]
    assert result.canonical_namespace == "entrez"
    assert result.canonical_identifier == "7157"
    assert result.entity_type == entity_type
    assert result.protein_identifier is None


@pytest.mark.parametrize(
    "namespace,identifier,expected_ns",
    [
        ("ensembl", "ENST00000141510.7", "enst"),
        ("enst", "ENST00000141510.7", "enst"),
        ("refseq", "NR_001.4", "refseq"),
        ("refseq", "NM_001.2", "refseq"),
    ],
)
def test_unmatched_transcript_retains_asserted_sequence_version(namespace, identifier, expected_ns):
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"r": obs("transcript", namespace, identifier)})["r"]
    assert result.entity_type == "transcript"
    assert result.transcript_namespace == expected_ns
    assert result.transcript_identifier == identifier
    assert result.protein_identifier is None
    assert result.gene_mapping_status == "missing"


def test_unmatched_transcript_form_is_retained_without_crossref_inference():
    molecule = obs("rna_product", "entrez", "missing")
    molecule.molecular_form = {"sequence_identifiers": [{"ns": "ensembl", "id": "ENST000001.5"}]}
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"r": molecule})["r"]
    assert result.transcript_identifier == "ENST000001.5"
    molecule.molecular_form = None
    molecule.identifiers = [{"ns": "enst", "id": "ENST000002.6"}]
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"r": molecule})["r"]
    assert result.transcript_identifier is None


@pytest.mark.parametrize("sequences", [None, []])
@pytest.mark.parametrize(
    "entity_type,namespace,identifier,form,protein,transcript",
    [
        (
            "protein",
            "uniprot",
            "Q9Y6K9",
            {"modifications": [{"residue": "S", "position": 15}]},
            ("uniprot", "Q9Y6K9"),
            None,
        ),
        (
            "protein",
            "uniprot",
            "Q9Y6K9-2",
            {"isoform_identifier": {"ns": "uniprot", "id": "Q9Y6K9-2"}},
            ("uniprot", "Q9Y6K9-2"),
            None,
        ),
        (
            "transcript",
            "ensembl",
            "ENST00000141510.7",
            {"variants": [{"reference": "A", "alternate": "G", "position": 12}]},
            None,
            ("enst", "ENST00000141510.7"),
        ),
        (
            "rna_product",
            "refseq",
            "NR_001.4",
            {"modifications": [{"description": "reported modification"}]},
            None,
            ("refseq", "NR_001.4"),
        ),
    ],
)
def test_serialized_nullable_form_collections_preserve_native_specificity(
    tmp_path,
    sequences,
    entity_type,
    namespace,
    identifier,
    form,
    protein,
    transcript,
):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from omnipath_core.molecular_forms import MOLECULAR_FORM_STRUCT

    # Arrow materializes omitted struct fields as explicit None. Test the
    # round-tripped shape consumed by the matcher, including null feature lists.
    table = pa.Table.from_pylist(
        [{"molecular_form": {**form, "sequence_identifiers": sequences}}],
        schema=pa.schema([("molecular_form", MOLECULAR_FORM_STRUCT)]),
    )
    path = tmp_path / "serialized-form.parquet"
    pq.write_table(table, path)
    serialized = pq.read_table(path).to_pylist()[0]["molecular_form"]
    assert serialized["sequence_identifiers"] == sequences
    assert serialized.get("variants") is None or serialized.get("modifications") is None
    molecule = obs(entity_type, namespace, identifier)
    molecule.molecular_form = serialized
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"participant": molecule})["participant"]
    assert not result.matched
    assert result.entity_type == entity_type
    assert (result.protein_namespace, result.protein_identifier) == (protein or (None, None))
    assert (result.transcript_namespace, result.transcript_identifier) == (
        transcript or (None, None)
    )
    assert result.protein_node_id is None
    assert result.protein_gene_candidates == ()
    assert molecule.molecular_form == serialized


def test_catalogued_protein_accepts_serialized_null_sequence_identifiers(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from omnipath_core.molecular_forms import MOLECULAR_FORM_STRUCT

    path = tmp_path / "protein-form.parquet"
    pq.write_table(
        pa.Table.from_pylist(
            [{"molecular_form": {"isoform_identifier": {"ns": "uniprot", "id": "P04637-2"}}}],
            schema=pa.schema([("molecular_form", MOLECULAR_FORM_STRUCT)]),
        ),
        path,
    )
    molecule = obs("protein", "uniprot", "P04637")
    molecule.molecular_form = pq.read_table(path).to_pylist()[0]["molecular_form"]
    assert molecule.molecular_form["sequence_identifiers"] is None
    with closing(LibraryMatcher(build_fixture_library(tmp_path))) as matcher:
        result = matcher.match({"participant": molecule})["participant"]
    assert result.matched
    assert (result.canonical_namespace, result.canonical_identifier) == ("entrez", "7157")
    assert (result.protein_namespace, result.protein_identifier) == ("uniprot", "P04637")
    assert result.protein_node_id == "uniprot:P04637"
    assert result.transcript_identifier is None
    assert molecule.molecular_form["isoform_identifier"]["id"] == "P04637-2"


def test_isoform_anchor_without_primary_parent_is_not_primary_product():
    facts = [(1, "uniprot:P04637-2", 2, "uniprot:P04637-2", False, True, [])]
    result = resolve_molecular_batch([("i", 2, [v(b"isoform")])], [(b"isoform", [1])], facts)[0]
    assert result[2] == []
    assert result[4] is None


def test_known_gene_and_missing_conflicting_gene_are_not_treated_as_alternatives():
    result = decide(
        [
            v(b"protein"),
            (b"known", "entrez:7157", 2, 1, True, False),
            (b"missing", "entrez:999999", 2, 2, True, False),
        ],
        [(b"protein", [1]), (b"known", [3]), (b"missing", [])],
    )
    assert result[2] == ["uniprot:P04637"]
    assert result[6] == "conflict"


def test_generic_refseq_protein_prefix_uses_product_lookup():
    from omnipath_resolver.canonical.match import votes_for
    from omnipath_resolver.canonical.policy import get_policy

    normalized, observed = votes_for(obs("protein", "refseq", "NP_001.4"), get_policy("protein"))
    assert [(v.ns, v.id) for v in normalized] == [("refseq_protein", "NP_001")]
    assert observed["refseq_protein"] == ["NP_001"]


def test_refseq_transcript_gene_mapping_uses_stable_accession_retaining_sequence_version(tmp_path):
    with closing(LibraryMatcher(build_fixture_library(tmp_path))) as matcher:
        result = matcher.match({"r": obs("noncoding_rna_product", "refseq", "NR_000055.9")})["r"]
    assert (result.canonical_namespace, result.canonical_identifier) == ("entrez", "55")
    assert (result.transcript_namespace, result.transcript_identifier) == ("refseq", "NR_000055.9")
    assert result.protein_identifier is None
    assert result.entity_type == "noncoding_rna_product"


def test_rna_subtype_replay_bundle_cannot_assert_catalogue_protein():
    from omnipath_resolver.canonical.match import votes_for
    from omnipath_resolver.canonical.policy import get_policy
    from omnipath_resolver.observations import observation_bundle

    molecule = obs("microrna", "uniprot", "P04637")
    normalized, observed = votes_for(molecule, get_policy(molecule.entity_type))
    _, rows = observation_bundle("r", molecule, normalized, observed, "gene_protein")
    assert rows and all(row["gene_only"] for row in rows)


def test_gene_alias_with_many_products_resolves_to_the_gene(tmp_path):
    with closing(LibraryMatcher(build_fixture_library(tmp_path))) as matcher:
        for namespace, identifier in [
            ("genesymbol", "MANYPRODUCTS"),
            ("ensg", "ENSG00000000777"),
            ("hgnc", "HGNC:777"),
            ("entrez", "777"),
        ]:
            result = matcher.match({"g": obs("protein", namespace, identifier)})["g"]
            assert result.canonical_namespace == "entrez"
            assert result.canonical_identifier == "777"
            assert result.protein_identifier is None
        result = matcher.match({"p": obs("protein", "uniprot", "P80000")})["p"]
        assert result.canonical_identifier == "777"
        assert result.protein_identifier == "P80000"


@pytest.mark.parametrize(
    "namespace,identifier,expected_namespace",
    [
        ("uniprot", "Q9Y6K9", "uniprot"),
        ("uniprot", "Q9Y6K9-2", "uniprot"),
        ("uniprot", "Q9Y6K9-PRO_000001", "uniprot"),
        ("uniprot", "Q9Y6K9-2-PRO_000001", "uniprot"),
        ("ensembl", "ENSP999999999999.73", "ensp"),
        ("ensp", "ENSP999999999999.73", "ensp"),
        ("refseq", "NP_999999999999.73", "refseq"),
        ("refseq_protein", "NP_999999999999.73", "refseq_protein"),
    ],
)
def test_unmatched_reported_protein_retains_exact_native_identity(
    namespace, identifier, expected_namespace
):
    from omnipath_core.molecular_forms import molecular_form_from_identifiers

    molecule = obs("protein", namespace, identifier)
    molecule.molecular_form = molecular_form_from_identifiers([{"ns": namespace, "id": identifier}])
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"p": molecule})["p"]
    assert not result.matched
    assert result.resolved_by == "unmatched"
    assert (result.protein_namespace, result.protein_identifier) == (expected_namespace, identifier)
    assert result.protein_node_id is None
    assert result.protein_gene_candidates == ()
    assert result.protein_aliases == {expected_namespace: [identifier]}
    assert result.protein_label == identifier


@pytest.mark.parametrize(
    "entity_type,namespace,identifier",
    [
        ("gene", "uniprot", "Q9Y6K9"),
        ("physical_entity", "uniprot", "Q9Y6K9"),
        ("rna_product", "uniprot", "Q9Y6K9"),
        ("protein", "entrez", "7157"),
        ("protein", "genesymbol", "Q9Y6K9"),
        ("protein", "uniprot", "NOT_A_VALID_ACCESSION"),
        ("protein", "enst", "ENST999999999999.73"),
        ("protein", "refseq", "NM_999999999999.73"),
    ],
)
def test_unmatched_gene_symbol_or_generic_type_does_not_assert_protein(
    entity_type, namespace, identifier
):
    molecule = obs(entity_type, namespace, identifier)
    molecule.identifiers = [{"ns": "uniprot", "id": "Q9Y6K9"}]
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"p": molecule})["p"]
    assert result.protein_identifier is None


def test_ambiguous_explicit_native_protein_forms_do_not_choose_representative():
    molecule = obs("protein", "entrez", "7157")
    molecule.molecular_form = {
        "sequence_identifiers": [
            {"ns": "uniprot", "id": "Q9Y6K9"},
            {"ns": "uniprot", "id": "Q9Y6K8"},
        ]
    }
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"p": molecule})["p"]
    assert result.protein_identifier is None


def test_native_product_and_gene_matching_have_independent_catalogue_status(tmp_path):
    molecule = obs("protein", "uniprot", "Q9Y6K9")
    molecule.identifiers = [{"ns": "entrez", "id": "7157"}]
    with closing(LibraryMatcher(build_fixture_library(tmp_path))) as matcher:
        results = matcher.match(
            {
                "unknown": obs("protein", "uniprot", "Q9Y6K9"),
                "unknown_gene": molecule,
                "known_without_gene": obs("protein", "uniprot", "Q99999"),
            }
        )
    assert not results["unknown"].matched
    assert results["unknown"].protein_identifier == "Q9Y6K9"
    assert results["unknown"].protein_node_id is None
    assert results["unknown_gene"].matched
    assert results["unknown_gene"].canonical_identifier == "7157"
    assert results["unknown_gene"].protein_identifier == "Q9Y6K9"
    assert results["unknown_gene"].protein_node_id is None
    assert results["unknown_gene"].protein_gene_candidates == ()
    known = results["known_without_gene"]
    assert known.matched and known.protein_identifier == "Q99999"
    assert known.protein_node_id == "uniprot:Q99999"
    assert known.gene_mapping_status == "missing"


@pytest.mark.parametrize("isoform", ["Q9Y6K9-2", "Q9Y6K8-2"])
def test_uncatalogued_principal_and_explicit_isoform_remain_conservative(isoform):
    molecule = obs("protein", "uniprot", "Q9Y6K9")
    molecule.molecular_form = {"isoform_identifier": {"ns": "uniprot", "id": isoform}}
    with closing(LibraryMatcher(None)) as matcher:
        result = matcher.match({"p": molecule})["p"]
    # A same-base entry/isoform pair is retained without claiming that the
    # entry is the registry's primary product; different products also stay apart.
    assert result.protein_identifier is None
    assert molecule.molecular_form["isoform_identifier"]["id"] == isoform
