import pytest
from pypath.inputs_v2.tcdb import _transporters_schema
from omnipath_build.silver import SilverExtractor


@pytest.mark.parametrize(
    "accession,namespace",
    [
        ("P12345", "uniprot"),
        ("A0A011NL42", "uniprot"),
        ("XP_026680479.1", "refseq_protein"),
        ("KQC09679.1", "genbank"),
    ],
)
def test_tcdb_sequence_accessions_are_not_uniprot(accession, namespace):
    row = dict(uniprot=accession, tcid="1.A.1")
    ex = SilverExtractor("tcdb", "transporters")
    ex.process_record(_transporters_schema(row), row, "transporters:0", 0)
    protein = next(iter(ex.entities.values()))
    assert protein.namespace == namespace and protein.identifier == accession


@pytest.mark.parametrize(
    "accession,namespace", [("WP_090136218", "refseq_protein"), ("ABC12345.2", "genbank")]
)
def test_brenda_uses_same_sequence_namespace_rule(accession, namespace):
    from pypath.inputs_v2.brenda import schema

    row = {"UniProt": accession}
    ex = SilverExtractor("brenda", "data")
    ex.process_record(schema(row), row, "data:0", 0)
    protein = next(iter(ex.entities.values()))
    assert protein.namespace == namespace and protein.identifier == accession
