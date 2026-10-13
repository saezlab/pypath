"""
Parse the in silico human surfaceome (SURFY) and emit protein-to-list and protein-to-class records.

Bausch-Fluck et al. 2018 trained a random forest classifier, SURFY, on the
cell surface proteins of the Cell Surface Protein Atlas and predicted a human
surfaceome of 2,886 proteins. Table S3 of the paper has these proteins with the
prediction (score, false positive rate class, how the label was reached, the
membrane topology) and, for most, the class that Almen et al. 2009 gave the
protein.

The list is one term, the surfaceome, and each protein is related to it with the
fields of the prediction on the statement. The Almen classes are the terms of
``almen2009.py`` (namespace ``almen2009_class``, EC classes for enzymes; their
mapping is shared in ``_almen_classes.py``), and each protein is related to its
class as in that module; the classes are emitted here again because the table is
their source for these proteins.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._almen_classes import (
    flag,
    map_class,
    protein_entity,
    surfaceome_table,
    term_entity,
)
from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Resource, ResourceConfig
from pypath.inputs_v2.parsers.almen2009 import CLASS
from pypath.inputs_v2.parsers.surfaceome import SURFACEOME, iter_proteins, iter_terms
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation

config = ResourceConfig(
    id=ResourceCv.SURFACEOME,
    name='Surfaceome',
    full='In silico human surfaceome',
    url='https://wollscheidlab.org/SURFY/',
    # The paper (PNAS 2018): "This open access article is distributed under Creative
    # Commons Attribution-NonCommercial-NoDerivatives License 4.0 (CC BY-NC-ND)."
    # The website states no license of its own. The legacy registry
    # (resources.json) says CC BY 4.0.
    license=LicenseCV.CC_BY_NC_ND_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='30373828',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Human proteins predicted to be at the cell surface, with the '
        'prediction score, topology and membrane protein class.'
    ),
    description=(
        'The in silico human surfaceome is the set of 2,886 human proteins '
        'that SURFY, a random forest classifier trained on the experimentally '
        'verified cell surface proteins of the Cell Surface Protein Atlas, '
        'predicts to be at the cell surface. The paper gives for each protein '
        'the prediction, its membrane topology, and the class of the human '
        'membrane proteome classification of Almen et al.'
    ),
)

# Table S3 of the paper (``almen2009.py`` reads its other sheet).
download = surfaceome_table

# The list of predicted surface proteins: a plain list, so one term.
LIST_NAMESPACE = 'surfaceome_list'
LIST_ID = 'in_silico_surfaceome'

# Source attributes of the statement that puts a protein in the list.
# The classifier's score is a Measurement without a unit: higher for a protein
# more likely to be at the surface, between 0 and 1 in the sheet. The number of
# transmembrane domains is a Measurement too (a count, no unit).
_SCORE = 'surfaceome:score'
_LABEL_SOURCE = 'surfaceome:label_source'
_TRAINING_SET = 'surfaceome:training_set'
_TM_DOMAINS = 'surfaceome:transmembrane_domains'
_SIGNAL_PEPTIDE = 'surfaceome:signal_peptide'
_TOPOLOGY = 'surfaceome:topology'
_TOPOLOGY_SOURCE = 'surfaceome:topology_source'

# "MachineLearning FPR class (1=1%, 2=5%, 3=15%)": the false positive rate the
# score is cut at. Some values carry ``(nonTM)``, which the table does not explain.
_FPR = {'1': 1, '2': 5, '3': 15}


def _confidence_level(value: str) -> str:
    code, _, note = value.partition(' ')
    if code not in _FPR or note not in ('', '(nonTM)'):
        raise ValueError(f'Unknown Surfaceome FPR class: {value!r}')
    return f'class {value}: false positive rate {_FPR[code]}%'


def _number(row: dict[str, Any], field: str, column: str):
    value = measurement(row[field], source_field=column)
    if value is None:
        raise ValueError(f'Surfaceome {row["uniprot"]}: {column} is {row[field]!r}')
    return value


def _list_entity() -> Entity:
    return Entity(
        type=OntologyClass,
        identifiers=[Identifier(type=LIST_NAMESPACE, value=LIST_ID)],
    )


def _prediction(row: dict[str, Any]) -> list[Annotation]:
    annotations = [
        Annotation(term=_SCORE, value=_number(row, 'score', 'MachineLearning score')),
        Annotation(term=_LABEL_SOURCE, value=row['label_source']),
        Annotation(
            term=_TM_DOMAINS, value=_number(row, 'tm_domains', 'TM domains')
        ),
        Annotation(term=_SIGNAL_PEPTIDE, value=flag(row['signal_peptide'])),
    ]
    for term, field in (
        (_TRAINING_SET, 'trainingset'),
        (_TOPOLOGY, 'topology'),
        (_TOPOLOGY_SOURCE, 'topology_source'),
    ):
        if row[field]:
            annotations.append(Annotation(term=term, value=row[field]))
    if row['fpr_class']:
        annotations.append(
            Annotation(
                term=slots.has_confidence_level,
                value=_confidence_level(row['fpr_class']),
            )
        )
    return annotations


def map_protein(row: dict[str, Any]) -> Relation:
    """A protein in the list, or in a class of Almen et al.: one record each."""
    if row['kind'] == SURFACEOME:
        object_, annotations = _list_entity(), _prediction(row)
    else:
        object_, annotations = term_entity(row['kind'], row['term']), []
    return Relation(
        subject=protein_entity(row),
        # A prediction and a classification place the protein in a group; they
        # state no cause or function. The prediction is of a protein at the
        # cell surface, which would be `located_in`, but a predicted membership
        # of a list is what the table gives.
        predicate=slots.associated_with,
        object=object_,
        annotations=annotations,
    )


def map_term(row: dict[str, Any]) -> Entity:
    """The list of the surfaceome, or a class of Almen et al."""
    if row['kind'] == SURFACEOME:
        return Entity(
            type=OntologyClass,
            identifiers=[
                Identifier(type=LIST_NAMESPACE, value=LIST_ID),
                Identifier(type=Namespace.NAME, value='In silico human surfaceome'),
            ],
            annotations=[
                Annotation(
                    term=slots.description,
                    value=(
                        'Human proteins that the random forest classifier SURFY '
                        'predicts to be at the cell surface.'
                    ),
                )
            ],
        )
    if row['kind'] == CLASS:
        return map_class(row)
    raise ValueError(f'Unknown kind of Surfaceome term: {row["kind"]!r}')


resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=map_term,
        raw_parser=iter_terms,
    ),
    proteins=Dataset(
        download=download,
        mapper=map_protein,
        raw_parser=iter_proteins,
    ),
)
