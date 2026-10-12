"""
Parse DEPOD phosphatase-substrate data and emit Relation records.

DEPOD-2019 publishes one worksheet of human phosphatases with their protein
substrates and one with their non-protein substrates. A protein substrate row
lists the dephosphorylated residues of one phosphatase-substrate pair
("Ser-339, Ser-340"), so a row becomes one record per residue. Rows without a
site ("N/A") and rows whose site cell holds free text ("Ser-660 in Isoform
Beta-II") stay single records: the first has no site, the second keeps the
source text without a position.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.depod import (
    ASSAY,
    ENZYME,
    REFERENCES,
    SCORE,
    SUBSTRATE,
    iter_non_protein_substrates,
    iter_protein_substrates,
    reference_labels,
)
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

_BASE_URL = 'https://depod.zellbiologie.uni-bonn.de/download/'

_UNIPROT_REGEX = (
    r'^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?)$'
)
# The literature column holds PubMed IDs or, for a UniProt annotation, this label.
_REFERENCE_DATABASES = {'UniProt'}
# DEPOD is the human dephosphorylation database: every phosphatase and protein
# substrate accession of the 2019 worksheets is a human UniProt entry. The
# worksheets have no organism column.
DEPOD_TAX_ID = '9606'


# Evidence & Conclusion Ontology terms of the assay type column ("Bioassay
# type" on the DEPOD site: "the type of experiments (in vitro or in vivo or
# both) showing dephosphorylation in original publications"). The column
# also holds "In vitro and ex vivo" in the non-protein worksheet. The source
# writes the same values with and without a capital letter. "N/A" gives none.
_IN_VIVO = 'ECO:0000178'  # in vivo assay evidence
_IN_VITRO = 'ECO:0000181'  # in vitro assay evidence
_EX_VIVO = 'ECO:0005502'  # ex vivo assay evidence
_ASSAY_TERMS = {
    'in vitro': (_IN_VITRO,),
    'in vivo': (_IN_VIVO,),
    'in vitro and in vivo': (_IN_VITRO, _IN_VIVO),
    'in vitro and ex vivo': (_IN_VITRO, _EX_VIVO),
    'n/a': (),
}


def _evidence_types(row: dict[str, Any]) -> list[str]:
    assay = str(row.get(ASSAY) or '').strip().lower()
    if assay not in _ASSAY_TERMS:
        raise ValueError(f'Unknown DEPOD assay type: {row.get(ASSAY)!r}')
    return list(_ASSAY_TERMS[assay])


def _data_sources(row: dict[str, Any]) -> list[str]:
    """Databases the literature column names instead of a publication."""
    labels = reference_labels(row.get(REFERENCES))
    unknown = sorted(set(labels) - _REFERENCE_DATABASES)
    if unknown:
        raise ValueError(f'Unknown DEPOD literature source: {unknown}')
    return labels


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The dephosphorylated residue, a feature of the substrate in this record.

    The position is the source's, on the UniProt entry the row names as the
    substrate. DEPOD states: "The dephosphorylation site information in case
    of protein substrates, has been corrected to the corresponding amino acid
    positions in the most recent reviewed entry for that protein in the
    latest release of UniProt (presently November 2018)" (About DEPOD,
    section 2.2), and DEPOD "only considers reviewed entries of proteins".
    It does not say that the position is on the canonical sequence, and the
    position is not checked against any sequence here. A site cell with free
    text has no position; its text is the description. The description is
    the site and, where the sheet lists one window per site, the residues
    around it as written (the dephosphorylated residue in lower case).
    """
    if not row.get('site') or not row.get(SUBSTRATE):
        return None
    position = row.get('position')
    return {
        'modifications': [
            {
                'term': 'de-phosphorylated residue',
                'residue': row.get('residue'),
                'position': position,
                'end_position': position,
                'coordinate_reference': {
                    'identifier': {'ns': 'uniprot', 'id': row[SUBSTRATE]},
                    'coordinate_system': 'protein',
                    'position_base': 1,
                }
                if position is not None
                else None,
                'description': f'{row["site"]} ({row["window"]})'
                if row.get('window')
                else row['site'],
            }
        ],
    }


config = ResourceConfig(
    id=ResourceCv.DEPOD,
    name='DEPOD',
    url='https://depod.zellbiologie.uni-bonn.de/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='31836896',
    topics=(Topic.SIGNALING, Topic.METABOLISM),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human phosphatases with their protein and non-protein substrates, '
        'dephosphorylation sites and in vitro or in vivo assay types.'
    ),
    description=(
        'DEPOD, the human DEPhOsphorylation Database, is a manually curated '
        'collection of human active phosphatases and their experimentally '
        'verified protein and non-protein substrates. A protein is a '
        'substrate only if in vitro or in vivo experiments show its direct '
        'dephosphorylation. DEPOD gives the dephosphorylation sites, the '
        'assay type, the publications and a reliability score for each '
        'phosphatase-substrate pair.'
    ),
)

f = FieldConfig(extract={'uniprot': _UNIPROT_REGEX, 'chebi': r'^CHEBI:(\d+)$'})


def _protein(accession_column: str, name_column: str, molecular_form=None) -> EntityBuilder:
    return EntityBuilder(
        entity_type=model.Protein,
        molecular_form=molecular_form,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f(accession_column, extract='uniprot')),
            CV(term=Namespace.NAME, value=f(name_column)),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{DEPOD_TAX_ID}'),
        ),
    )


def _dephosphorylation_annotations() -> AnnotationsBuilder:
    return AnnotationsBuilder(
        # The phosphatase decreases the substrate's phosphorylation; DEPOD
        # says nothing here about the substrate's activity.
        CV(
            term=slots.object_direction_qualifier,
            value=model.DirectionQualifierEnum.decreased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=model.GeneOrGeneProductOrChemicalEntityAspectEnum.phosphorylation,
        ),
        CV(
            term=slots.publications,
            value=lambda row: [f'PMID:{pmid}' for pmid in row.get('pmids') or ()],
        ),
        CV(term=slots.supporting_data_source, value=_data_sources),
        CV(term=slots.has_evidence_of_type, value=_evidence_types),
        # 1 to 3, from the assay types and the number of laboratories.
        CV(
            term=slots.has_confidence_score,
            value=lambda row: measurement(row.get(SCORE), source_field=SCORE),
        ),
    )


enzyme_substrate_schema = RelationBuilder(
    subject=_protein(ENZYME, 'Phosphatase entry names'),
    predicate=slots.affects,
    # The dephosphorylated residue is in the substrate's molecular form.
    object=_protein(SUBSTRATE, 'Substrate entry names', molecular_form=_substrate_form),
    annotations=_dephosphorylation_annotations(),
)

non_protein_substrates_schema = RelationBuilder(
    subject=_protein(ENZYME, 'Phosphatase entry name'),
    predicate=slots.affects,
    object=EntityBuilder(
        entity_type=model.ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.CHEBI, value=f('Substrate ChEBI ID', extract='chebi')),
            CV(
                term=Namespace.KEGG,
                value=f('Substrate KEGG ID', extract=r'^(C\d{5})$'),
            ),
            CV(term=Namespace.NAME, value=f('Substrate ChEBI name')),
        ),
    ),
    annotations=_dephosphorylation_annotations(),
)


def _download(filename: str) -> Download:
    return Download(
        url=_BASE_URL + filename,
        filename=filename,
        subfolder='depod',
        large=True,
        encoding=None,
        default_mode='rb',
    )


resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=_download('PPase_protSubtrates_201903.xls'),
        mapper=enzyme_substrate_schema,
        raw_parser=iter_protein_substrates,
    ),
    non_protein_substrates=Dataset(
        download=_download('PPase_nonprotSubtrates_201903.xls'),
        mapper=non_protein_substrates_schema,
        raw_parser=iter_non_protein_substrates,
    ),
)
