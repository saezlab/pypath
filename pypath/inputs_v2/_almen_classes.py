"""The classes of Almen et al. 2009 as entities, shared by ``almen2009.py`` and ``surfaceome.py``.

Both resources carry the classification of the human membrane proteome (the
surfaceome table repeats the class of each protein), so a class must be the same
entity whichever resource emits it: same type, namespace, identifier and name.
The mapping lives here, once; the reading of the classification is in
``parsers/almen2009.py``.

Both also read the proteins from the same file, the table of the in silico
surfaceome (``parsers/surfaceome.py``): its download and the protein of one of
its rows are here too.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.almen2009 import CLASS, EC
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    OntologyRelation,
)

# Table S3 of Bausch-Fluck et al. 2018, from the web site of the Wollscheid lab
# (the old address, wlab.ethz.ch/surfaceome, no longer serves the file).
surfaceome_table = Download(
    url='https://wollscheidlab.org/SURFY/table_S3_surfaceome.xlsx',
    filename='table_S3_surfaceome.xlsx',
    subfolder='surfaceome',
    large=True,
    default_mode='rb',
    encoding=None,
)

# The table is of the human proteome.
TAXON = 'NCBITaxon:9606'


def protein_entity(row: dict[str, Any]) -> Entity:
    """The protein of a row of the surfaceome table, with the identifiers the
    table gives it: UniProt accession and entry name, gene symbol, Entrez Gene IDs."""
    identifiers = [
        Identifier(type=Namespace.UNIPROT, value=row['uniprot']),
        Identifier(type=Namespace.UNIPROT_ENTRY, value=row['entry_name']),
    ]
    if row['gene']:
        identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=row['gene']))
    identifiers.extend(
        Identifier(type=Namespace.ENTREZ, value=entrez) for entrez in row['entrez']
    )
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )

# The classes of the source: the main classes and the paths below them, each
# identified by the path as the source writes it (``GPCR;Rhodopsin``).
CLASS_NAMESPACE = 'almen2009_class'


def flag(value: bool) -> str:
    """A yes/no flag of the source as the string ``'true'`` or ``'false'``.

    The parser has read the file's own words (``Y`` or ``0`` of Phobius SP, ``1``
    or ``0`` of the surfaceome sheet) into a boolean; any other value raises.
    """
    if value is True:
        return 'true'
    if value is False:
        return 'false'
    raise ValueError(f'Not a flag: {value!r}')


def class_entity(identifier: str) -> Entity:
    """A class of the source's hierarchy; the same entity as in ``classes``."""
    return Entity(
        type=OntologyClass,
        identifiers=[Identifier(type=CLASS_NAMESPACE, value=identifier)],
    )


def ec_entity(number: str) -> Entity:
    """An EC class, written as the enzyme classification writes it."""
    return Entity(
        type=OntologyClass,
        identifiers=[Identifier(type=Namespace.EC, value=f'EC:{number}')],
    )


def term_entity(kind: str, term: str) -> Entity:
    """The class object of a statement, from the parser's ``kind`` and ``term``."""
    if kind == CLASS:
        return class_entity(term)
    if kind == EC:
        return ec_entity(term)
    raise ValueError(f'Unknown kind of Almen2009 class: {kind!r}')


def map_class(row: dict[str, Any]) -> Entity:
    """A class of the hierarchy, below its parent."""
    identifiers = [
        Identifier(type=CLASS_NAMESPACE, value=row['id']),
        Identifier(type=Namespace.NAME, value=row['name']),
    ]
    if row['label'] != row['name']:
        identifiers.append(Identifier(type=Namespace.SYNONYM, value=row['label']))
    parent = row['parent']
    return Entity(
        type=OntologyClass,
        identifiers=identifiers,
        ontology_relations=[
            OntologyRelation(
                predicate=slots.subclass_of,
                object=EntityRef(
                    type=OntologyClass,
                    identifier_type=CLASS_NAMESPACE,
                    identifier=parent,
                ),
                ontology_id=CLASS_NAMESPACE,
            )
        ]
        if parent
        else None,
    )
