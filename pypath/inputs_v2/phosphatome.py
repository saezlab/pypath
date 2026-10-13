"""
Parse the human phosphatases of Chen et al. 2017 and emit folds, families and substrate classes.

Table S1 of "Genomics and evolution of protein phosphatases" lists the 264
human phosphatase genes with their fold, family and subfamily (the
classification of Phosphatome.Net), and says for each whether it is a protein
phosphatase, has protein or non-protein substrates and has catalytic
activity. The module emits

- ``families``: each fold as a gene family with its families as members, each
  family with its subfamilies and the genes that have no subfamily, and each
  subfamily with its genes. A fold has the SCOP name of Table S2 as its
  description.
- ``terms`` and ``class_annotations``: the classes of the substrate-type
  columns as terms of the resource, and each gene associated with the classes
  that apply to it.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, GeneFamily, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phosphatome import (
    CATALYTIC_ACTIVITY,
    NON_PROTEIN_SUBSTRATES,
    PROTEIN_PHOSPHATOME,
    PROTEIN_SUBSTRATES,
    iter_class_values,
    iter_phosphatase_classes,
    iter_phosphatase_tree,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.PHOSPHATOME,
    name='Phosphatome',
    url='http://phosphatome.net/',
    # The site and the supplement state no license. License as recorded in the
    # legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='28400531',
    topics=(Topic.SIGNALING,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    # The article (Science Signaling 2017), abstract: "we cataloged the human
    # protein phosphatome, composed of 189 known and predicted human protein
    # phosphatase genes", "classified protein phosphatases [...] into a hierarchy of
    # 10 protein folds, 21 families, and 178 subfamilies" and "Sequence analysis
    # predicts that 26 human protein phosphatase domains are catalytically
    # disabled".
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Folds, families and subfamilies of the human phosphatases, and '
        'their substrate types and catalytic activity.'
    ),
    description=(
        'Phosphatome.Net and the study of Chen, Dixon and Manning classify '
        'the protein phosphatases and other phosphatases of the human genome '
        'by fold, family and subfamily, and record which of them act on '
        'protein or on non-protein substrates and whether they are catalytic.'
    ),
)

# The genes are those of the human genome.
TAXON = 'NCBITaxon:9606'
# The names of the classification levels have no Namespace member. A name is
# only unique within its level (the fold AP has the family AP).
NAMESPACES = {
    'fold': 'phosphatome_fold',
    'family': 'phosphatome_family',
    'subfamily': 'phosphatome_subfamily',
}

SET_NAMESPACE = 'phosphatome_set'
SUBSTRATE_NAMESPACE = 'phosphatome_substrate_type'
ACTIVITY_NAMESPACE = 'phosphatome_catalytic_activity'
# Each column value of the sheet -> the term it states, as (namespace,
# identifier, name), or None where the value states no class: the "No" of a
# yes/no column, and the two catalytic-activity values "Unknown" and "N/A".
# The identifier is the source's value (the column's name without "?" for the
# yes/no columns, where the value "Yes" means that the gene has it).
TERMS: dict[tuple[str, str], tuple[str, str, str] | None] = {
    (PROTEIN_PHOSPHATOME, 'Yes'): (
        SET_NAMESPACE, 'Protein phosphatome', 'Protein phosphatome',
    ),
    (PROTEIN_PHOSPHATOME, 'No'): None,
    (PROTEIN_SUBSTRATES, 'Yes'): (
        SUBSTRATE_NAMESPACE, 'Protein substrates', 'Protein substrates',
    ),
    (PROTEIN_SUBSTRATES, 'No'): None,
    (NON_PROTEIN_SUBSTRATES, 'Yes'): (
        SUBSTRATE_NAMESPACE, 'Non-protein substrates', 'Non-protein substrates',
    ),
    (NON_PROTEIN_SUBSTRATES, 'No'): None,
    (CATALYTIC_ACTIVITY, 'Yes'): (
        ACTIVITY_NAMESPACE, 'Yes', 'Catalytically active',
    ),
    (CATALYTIC_ACTIVITY, 'Verified pseudophosphatase'): (
        ACTIVITY_NAMESPACE, 'Verified pseudophosphatase', 'Verified pseudophosphatase',
    ),
    (CATALYTIC_ACTIVITY, 'Predicted pseudophosphatase'): (
        ACTIVITY_NAMESPACE, 'Predicted pseudophosphatase', 'Predicted pseudophosphatase',
    ),
    (CATALYTIC_ACTIVITY, 'Unknown'): None,
    (CATALYTIC_ACTIVITY, 'N/A'): None,
}
# The columns of a gene row of the ``class_annotations`` dataset and the sheet column
# each one holds.
CLASS_FIELDS = {
    'protein_phosphatome': PROTEIN_PHOSPHATOME,
    'protein_substrates': PROTEIN_SUBSTRATES,
    'non_protein_substrates': NON_PROTEIN_SUBSTRATES,
    'catalytic_activity': CATALYTIC_ACTIVITY,
}


def _term(column: str, value: Any) -> tuple[str, str, str] | None:
    """The term that a value of a class column states; unknown values raise."""
    if (column, value) not in TERMS:
        raise ValueError(f'Unknown Phosphatome value {value!r} in {column!r}')
    return TERMS[(column, value)]


def _gene(member: dict[str, Any]) -> Entity:
    return Entity(
        type=Gene,
        identifiers=[
            Identifier(type=Namespace.ENTREZ, value=member['entrez']),
            Identifier(type=Namespace.GENESYMBOL, value=member['symbol']),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def map_family(node: dict[str, Any]) -> Entity:
    """A fold, family or subfamily with the next level and its genes."""
    annotations = [Annotation(term=slots.in_taxon, value=TAXON)]
    if node['description']:
        # The SCOP identifier and name of the fold (Table S2).
        annotations.append(
            Annotation(term=slots.description, value=f"SCOP {node['description']}")
        )
    return Entity(
        type=GeneFamily,
        identifiers=[
            Identifier(type=NAMESPACES[node['level']], value=node['name']),
            Identifier(type=Namespace.NAME, value=node['name']),
        ],
        annotations=annotations,
        membership=[
            *(Membership(member=map_family(child)) for child in node['children']),
            *(Membership(member=_gene(member)) for member in node['members']),
        ],
    )


def map_term(row: dict[str, Any]) -> Entity | None:
    """A class of the substrate-type columns as a term of the resource."""
    term = _term(row['column'], row['value'])
    if term is None:
        return None
    namespace, identifier, name = term
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=namespace, value=identifier),
            Identifier(type=Namespace.NAME, value=name),
        ],
    )


def map_classes(row: dict[str, Any]) -> Entity | None:
    """A gene associated with each class that the sheet gives it.

    A gene that the sheet puts in no class (25 genes that are no protein
    phosphatase and have no substrate type and catalytic activity "N/A") has
    no statement and gives no record.
    """
    terms = [
        term
        for term in (
            _term(column, row[field]) for field, column in CLASS_FIELDS.items()
        )
        if term is not None
    ]
    if not terms:
        return None
    return Entity(
        type=Gene,
        identifiers=[
            Identifier(type=Namespace.ENTREZ, value=row['entrez']),
            Identifier(type=Namespace.GENESYMBOL, value=row['symbol']),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        associations=[
            Association(
                object=EntityRef(OntologyClass, namespace, identifier),
                predicate=slots.associated_with,
            )
            for namespace, identifier, _ in terms
        ],
    )


# The rescued copy of the supplement of Science Signaling (the journal serves
# it only to a browser; the legacy module used this copy).
download = Download(
    url='https://rescued.omnipathdb.org/phosphatome-supplementary-S1-S23.zip',
    filename='phosphatome-supplementary-S1-S23.zip',
    subfolder='phosphatome',
    large=True,
    ext='zip',
    needed=['aag1796_Tables S1 to S23.xlsx'],
    encoding=None,
    default_mode='rb',
)

resource = Resource(
    config,
    families=Dataset(
        download=download, mapper=map_family, raw_parser=iter_phosphatase_tree
    ),
    terms=Dataset(
        download=download, mapper=map_term, raw_parser=iter_class_values
    ),
    class_annotations=Dataset(
        download=download, mapper=map_classes, raw_parser=iter_phosphatase_classes
    ),
)
