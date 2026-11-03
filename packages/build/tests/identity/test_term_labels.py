"""Ontology hub -> term_labels.parquet -> the resolver names terms cited by id alone."""

import pyarrow.parquet as pq

from omnipath_build.hubs.sources.ontology import emit
from omnipath_build.hubs.writer import HubParquetWriter
from omnipath_build.identity.common import Context
from omnipath_build.identity.decisions import find_term_index, stage_term_labels
from omnipath_build.identity.hubindex import build_hub_index

GO = """format-version: 1.2

[Term]
id: GO:0001675
name: acrosome assembly
alt_id: GO:0000001

[Term]
id: GO:0000001
name: mitochondrion inheritance

[Term]
id: GO:0000002
name: obsolete term
is_obsolete: true
"""
CHEMONT = """[Term]
id: CHEMONTID:0000002
name: Organoheterocyclic compounds
"""
BTO = """[Term]
id: BTO:0000007
name: HEK-293 cell

[Term]
id: BTO:0001416
name: umbilical vein endothelium
alt_id: BTO:0001414
"""
# a merged edition carries the terms of the ontologies it imports: they are not CL's
CL = """[Term]
id: CL:0002518
name: kidney epithelial cell

[Term]
id: UBERON:0002113
name: kidney (as CL imports it)

[Term]
id: CL:0000002
name: obsolete immortal cell line cell
is_obsolete: true
"""
# APO names its former ids (YPO:) as alt_ids; the file holds no foreign terms but the test pins the own-prefix rule
APO = """[Term]
id: APO:0000063
name: colony size
alt_id: YPO:0000063

[Term]
id: GO:0008150
name: biological_process

[Term]
id: APO:0000064
name: liquid culture appearance
is_obsolete: true
"""
UBERON = """[Term]
id: UBERON:0002113
name: kidney
"""
DOID = """[Term]
id: DOID:3702
name: cervical adenocarcinoma

[Term]
id: DOID:8552
name: chronic myeloid leukemia
alt_id: DOID:8551
"""
# pato-base.obo carries an unnamed stub of a Cell Ontology term; pato.obo does not. The rule is the same
PATO = """[Term]
id: PATO:0000169
name: viability

[Term]
id: CL:0000000
relationship: BFO:0000051 GO:0005634
"""
# eco.obo holds ECO's own terms only; the foreign stanza is the one a merged edition would add
ECO = """[Term]
id: ECO:0000178
name: in vivo assay evidence

[Term]
id: ECO:0000181
name: in vitro assay evidence
alt_id: ECO:0005013

[Term]
id: ECO:0000039
name: obsolete protein assay evidence
is_obsolete: true

[Term]
id: OBI:0000623
name: a term of another ontology
"""
MP = """[Term]
id: MP:0001262
name: decreased body weight

[Term]
id: MP:0000012
name: obsolete loss of subcutaneous adipose tissue
is_obsolete: true
"""
# the prefix is mixed case; a term of another case is not WBPhenotype's
WBPHENOTYPE = """[Term]
id: WBPhenotype:0000050
name: embryonic lethal

[Term]
id: WBPHENOTYPE:0000001
name: not a WBPhenotype term
"""


def test_term_labels_from_the_ontology_hub(tmp_path):
    hubs = tmp_path / "hubs"
    hubs.mkdir()
    writer = HubParquetWriter(hubs / "ontology.parquet")
    # streamed lines keep their line ends
    emit(
        writer,
        {
            "go": GO.splitlines(keepends=True),
            "chemont": CHEMONT.splitlines(keepends=True),
            "hpo": [],
            "mondo": [],
            "psi_mi": [],
            "bto": BTO.splitlines(keepends=True),
            "cl": CL.splitlines(keepends=True),
            "uberon": UBERON.splitlines(keepends=True),
            "apo": APO.splitlines(keepends=True),
            "doid": DOID.splitlines(keepends=True),
            "pato": PATO.splitlines(keepends=True),
            "eco": ECO.splitlines(keepends=True),
            "mp": MP.splitlines(keepends=True),
            "wbphenotype": WBPHENOTYPE.splitlines(keepends=True),
        },
    )
    writer.close()
    build_hub_index("ontology", hubs, tmp_path / "index", memory="1GB", threads=1, min_free_gib=0)
    index = find_term_index(tmp_path / "index")
    out = tmp_path / "identity"
    stage_term_labels(Context(out, "1GB", 1, None, 0), out, index, out)
    labels = dict(zip(*pq.read_table(out / "term_labels.parquet").to_pydict().values()))
    # a term's own id wins over another term's alt_id; obsolete terms are not named
    assert labels == {
        "GO:0001675": "acrosome assembly",
        "GO:0000001": "mitochondrion inheritance",
        "CHEMONTID:0000002": "Organoheterocyclic compounds",
        "BTO:0000007": "HEK-293 cell",
        "BTO:0001416": "umbilical vein endothelium",
        "BTO:0001414": "umbilical vein endothelium",
        "CL:0002518": "kidney epithelial cell",
        # only Uberon names Uberon terms, not the CL file that imports them
        "UBERON:0002113": "kidney",
        "APO:0000063": "colony size",
        "YPO:0000063": "colony size",
        "DOID:3702": "cervical adenocarcinoma",
        "DOID:8552": "chronic myeloid leukemia",
        "DOID:8551": "chronic myeloid leukemia",
        "PATO:0000169": "viability",
        "ECO:0000178": "in vivo assay evidence",
        "ECO:0000181": "in vitro assay evidence",
        "ECO:0005013": "in vitro assay evidence",
        "MP:0001262": "decreased body weight",
        "WBPhenotype:0000050": "embryonic lethal",
    }


def test_resolver_names_terms_cited_by_id_alone():
    from omnipath_resolver.canonical.match import LibraryMatcher
    from omnipath_resolver.contracts import RawEntityObservation

    class Runtime:
        libraries = ()

        @staticmethod
        def term_label(term):
            return {"GO:0001675": "acrosome assembly"}.get(term)

        def close(self):
            pass

    matcher = LibraryMatcher(None)
    matcher.runtime = Runtime()
    observations = {
        "bare": RawEntityObservation(
            entity_key="bare", entity_type="ontology_class", namespace="go", identifier="GO:0001675"
        ),
        "unknown": RawEntityObservation(
            entity_key="unknown", entity_type="ontology_class", namespace="go", identifier="GO:9999999"
        ),
    }
    results = matcher.targets(observations)
    assert results["bare"][0].label == "acrosome assembly"
    assert results["unknown"][0].label == "GO:9999999"
