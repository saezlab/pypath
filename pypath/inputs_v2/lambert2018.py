"""
Parse the human transcription factors of Lambert et al. 2018 and emit the list and its families.

Lambert, Jolma, Campitelli et al. reviewed 2,765 candidate genes by hand and
decided for each whether it is a sequence-specific DNA-binding transcription
factor (TF). Table S1 gives the verdict ("Is TF?"), the DNA-binding domain
(DBD) family, the likelihood that the gene is a TF ("TF assessment") and the
binding mode, with the assessments of the two curators. The module emits

- ``terms``: the list of human TFs as one term;
- ``transcription_factors``: each gene with the verdict ``Yes`` as a member of
  the list, with the assessment as confidence level and the binding mode,
  motif status, notes and the curators' assessments as attributes;
- ``families``: the DBD families as gene families with their TF genes.

The genes with the verdict ``No`` are statements that the gene is not a TF:
they are not members of the list and are not emitted.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, GeneFamily, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.lambert2018 import (
    iter_dbd_families,
    iter_list,
    iter_transcription_factors,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    Membership,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.LAMBERT2018,
    name='Lambert2018',
    full='The Human Transcription Factors',
    url='https://doi.org/10.1016/j.cell.2018.01.029',
    # License as recorded in the legacy registry (resources.json): the
    # Elsevier user license. The article page was not readable to confirm
    # what it says today.
    license=LicenseCV.ELSEVIER_USER_LICENSE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='29425488',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human transcription factor genes with the curators\' assessment, '
        'binding mode and DNA-binding domain family.'
    ),
    description=(
        'The Human Transcription Factors is a manually curated list of the '
        'human genes that encode sequence-specific DNA-binding transcription '
        'factors. Expert curators reviewed 2,765 candidate genes, assessed '
        'how likely each is to be a TF, with emphasis on whether it has a '
        'known motif, and how it binds DNA.'
    ),
)

# The genes are human.
TAXON = 'NCBITaxon:9606'
# The list and the DBD families have no Namespace member.
LIST_NAMESPACE = 'lambert2018_list'
LIST_ID = 'human_transcription_factors'
FAMILY_NAMESPACE = 'lambert2018_dbd_family'
# Attributes of a gene's statement, with the source's words as the values.
BINDING_MODE = 'lambert2018:binding_mode'
MOTIF_STATUS = 'lambert2018:motif_status'
CURATOR_ASSESSMENT = 'lambert2018:curator_assessment'
CURATOR_BINDING_MODE = 'lambert2018:curator_binding_mode'
ASSESSMENT_DISAGREEMENT = 'lambert2018:assessment_disagreement'
BINDING_DISAGREEMENT = 'lambert2018:binding_disagreement'
# Where the source gives the ID of a gene as a UniProt entry name.
ID_NAMESPACES = {
    'ensg': Namespace.ENSG,
    'uniprot_entry': Namespace.UNIPROT_ENTRY,
}
# "TF assessment" of the genes with the verdict Yes, the confidence level of
# the statement. The other two values are of genes with the verdict No.
ASSESSMENTS = frozenset(
    {'Known motif', 'Inferred motif', 'Likely to be sequence specific TF'}
)


def map_list(row: dict[str, Any]) -> Entity:
    """The list of human TFs as a term."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=LIST_NAMESPACE, value=LIST_ID),
            Identifier(
                type=Namespace.NAME,
                value='Human transcription factors (Lambert et al. 2018)',
            ),
        ],
        annotations=[
            Annotation(
                term=slots.description,
                value=(
                    'Human genes that the curators of Lambert et al. 2018 '
                    'judged to encode sequence-specific DNA-binding '
                    'transcription factors ("Is TF?" Yes: '
                    f"{row['verdicts']['Yes']} of the "
                    f"{sum(row['verdicts'].values())} candidates)."
                ),
            )
        ],
    )


def _gene(row: dict[str, Any]) -> Entity:
    return Entity(
        type=Gene,
        identifiers=[
            Identifier(type=ID_NAMESPACES[row['id_type']], value=row['id']),
            Identifier(type=Namespace.GENESYMBOL, value=row['symbol']),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def map_transcription_factor(row: dict[str, Any]) -> Relation | None:
    """A gene with the verdict Yes as a member of the list."""
    if row['is_tf'] != 'Yes':
        return None
    if row['assessment'] not in ASSESSMENTS:
        raise ValueError(f'Unknown Lambert 2018 assessment of a TF: {row!r}')
    annotations = [Annotation(term=slots.has_confidence_level, value=row['assessment'])]
    for term, field in (
        (BINDING_MODE, 'binding_mode'),
        (MOTIF_STATUS, 'motif_status'),
        (slots.description, 'notes'),
        (CURATOR_ASSESSMENT, 'assessment_1'),
        (CURATOR_ASSESSMENT, 'assessment_2'),
        (CURATOR_BINDING_MODE, 'binding_1'),
        (CURATOR_BINDING_MODE, 'binding_2'),
        (ASSESSMENT_DISAGREEMENT, 'assessment_disagreement'),
        (BINDING_DISAGREEMENT, 'binding_disagreement'),
    ):
        if row[field] is not None:
            annotations.append(Annotation(term=term, value=row[field]))
    return Relation(
        subject=_gene(row),
        predicate=slots.associated_with,
        object=EntityRef(OntologyClass, LIST_NAMESPACE, LIST_ID),
        annotations=annotations,
    )


def map_family(row: dict[str, Any]) -> Entity:
    """A DNA-binding domain family with its TF genes."""
    return Entity(
        type=GeneFamily,
        identifiers=[
            Identifier(type=FAMILY_NAMESPACE, value=row['dbd']),
            Identifier(type=Namespace.NAME, value=row['dbd']),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        membership=[Membership(member=_gene(member)) for member in row['members']],
    )


# The rescued copy of Table S1 (mmc2.xlsx of the article); cell.com serves the
# supplement only to a browser. The legacy module used this copy.
download = Download(
    url='https://rescued.omnipathdb.org/lambert2018-supplementary-S1.xlsx',
    filename='lambert2018-supplementary-S1.xlsx',
    subfolder='lambert2018',
    large=True,
    encoding=None,
    default_mode='rb',
)

resource = Resource(
    config,
    terms=Dataset(download=download, mapper=map_list, raw_parser=iter_list),
    transcription_factors=Dataset(
        download=download,
        mapper=map_transcription_factor,
        raw_parser=iter_transcription_factors,
    ),
    families=Dataset(
        download=download, mapper=map_family, raw_parser=iter_dbd_families
    ),
)
