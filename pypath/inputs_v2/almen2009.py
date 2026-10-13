"""
Parse the classification of the human membrane proteome of Almen et al. 2009.

The paper mined the human proteome with three predictors of transmembrane
helices (Phobius, TMHMM and SOSUI), clustered the 6,718 membrane proteins it
found and classified them by curation into groups under five main classes:
enzymes, receptors, transporters, miscellaneous and unclassified. The
additional file lists each protein by its IPI accession with the main class and
the path of classes below it (``GPCR;Rhodopsin``, ``SLC;APC;SLC12``).

The classes are terms of their own, in the namespace ``almen2009_class``, read
from the additional file (``terms``). A protein is related to the most specific
class of its path; the classes above it follow from the ``subclass_of``
relations of the terms. The classes of enzymes are EC numbers, so an enzyme is
related to its EC classes in the EC namespace and to the main class ``Enzymes``.

The proteins (``proteins``) are not read from the additional file: its
accessions are those of the International Protein Index, which no longer
exists, and nothing resolves them. They are read from the table of the in
silico surfaceome (Bausch-Fluck et al. 2018, ``surfaceome.py``), whose first
sheet gives the class of Almen et al. to 6,264 human proteins by UniProt
accession. Compared with the additional file (6,717 proteins, 311 classes
with a protein), that sheet lacks 453 proteins and has no protein for six
classes: four groups of ``Unknown_function`` (37 proteins in the additional
file), EC 2.3.3 (1) and EC 1.10 (3). The 14 proteins whose EC class Excel
stored as the number ``1.1`` (EC 1.1 or EC 1.10) are in the main class
``Enzymes`` only (``parsers/surfaceome.py``). The predictions of the three
tools for each protein, which only the additional file has, are not carried.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import slots
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._almen_classes import (
    map_class,
    protein_entity,
    surfaceome_table,
    term_entity,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.almen2009 import iter_classes
from pypath.inputs_v2.parsers.surfaceome import iter_almen_classes
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Relation

config = ResourceConfig(
    id=ResourceCv.ALMEN2009,
    name='Almen2009',
    full='Classification of the human membrane proteome',
    url='https://doi.org/10.1186/1741-7007-7-50',
    # "This is an Open Access article distributed under the terms of the Creative
    # Commons Attribution License (http://creativecommons.org/licenses/by/2.0)",
    # BMC Biology 2009, 7:50 (the legacy registry says 4.0). The classification
    # is the data of Almen et al. and has their license, also where it is read
    # from the table of Bausch-Fluck et al. 2018 (whose article is CC BY-NC-ND
    # 4.0); decided by the project owner, 2026-10-13.
    license=LicenseCV.CC_BY_2_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='19678920',
    topics=(Topic.SUBCELLULAR_LOCALIZATION, Topic.MOLECULAR_CLASSIFICATION),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Classification of 6,264 human membrane proteins into receptor, '
        'transporter, enzyme and other families.'
    ),
    description=(
        'Almen et al. identified the human membrane proteome with three '
        'predictors of transmembrane alpha-helices (Phobius, TMHMM, SOSUI), '
        'reduced it to a non-redundant set by aligning it to the genome, '
        'clustered it, and curated the clusters against a range of databases. '
        'They classified the 6,718 proteins into families under the main '
        'classes enzymes, receptors, transporters, miscellaneous and '
        'unclassified. The paper identifies the proteins by IPI accessions, '
        'which no longer exist; the proteins are imported from the table of '
        'the in silico human surfaceome (Bausch-Fluck et al. 2018), which '
        'gives the class of Almen et al. to 6,264 proteins by UniProt accession.'
    ),
)

# The additional file of the paper, from the publisher: the classes.
download = Download(
    url=(
        'https://static-content.springer.com/esm/'
        'art%3A10.1186%2F1741-7007-7-50/MediaObjects/12915_2009_258_MOESM1_ESM.xls'
    ),
    filename='almen2009_12915_2009_258_MOESM1_ESM.xls',
    subfolder='almen2009',
    large=True,
    default_mode='rb',
    encoding=None,
)

def map_protein_class(row: dict[str, Any]) -> Relation:
    """A protein is associated with a class: one record per class."""
    return Relation(
        subject=protein_entity(row),
        # The classification makes no claim of function or cause: it places
        # the protein in a group.
        predicate=slots.associated_with,
        object=term_entity(row['kind'], row['term']),
        annotations=[],
    )


resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=map_class,
        raw_parser=iter_classes,
    ),
    proteins=Dataset(
        download=surfaceome_table,
        mapper=map_protein_class,
        raw_parser=iter_almen_classes,
    ),
)
