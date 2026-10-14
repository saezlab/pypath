"""The mapping shared by ExoCarta and Vesiclepedia.

Both databases publish the same two tables (see ``parsers/exocarta.py``), and
both are mapped the same way: a gene of the protein/RNA table, ``located_in``
the vesicle type of the study, with the evidence of the study. The modules of
the two resources keep what differs (the species of their tables, the
namespace and the description of their vesicle types) and pass it in.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2._publications import pmid_curie
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation


def map_detection(
    row: dict[str, Any], *, prefix: str, namespace: str, taxa: dict[str, str]
) -> Relation:
    """A molecule found in the vesicles of a type, in an experiment.

    ``prefix`` names the database in the source attributes and ``namespace``
    is the namespace of its vesicle types. A species that is not in ``taxa``
    raises.
    """
    if row['species'] not in taxa:
        raise ValueError(f'Unknown {prefix} species: {row["species"]!r}')
    identifiers = [Identifier(type=Namespace.ENTREZ, value=row['entrez'])]
    if row['symbol']:
        identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=row['symbol']))
    gene = Entity(
        type=Gene,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=taxa[row['species']])],
    )
    # The content type says whether the table lists a protein or an RNA of the gene.
    annotations = [
        Annotation(term=f'{prefix}:content_type', value=row['content_type']),
        Annotation(term=f'{prefix}:experiment', value=row['experiment']),
    ]
    for key in ('sample', 'sample_name'):
        if row[key]:
            annotations.append(Annotation(term=f'{prefix}:{key}', value=row[key]))
    for key, term in (
        ('sample_source', 'sample_source'),
        ('isolation_methods', 'isolation_method'),
        ('methods', 'detection_method'),
    ):
        annotations.extend(Annotation(term=f'{prefix}:{term}', value=v) for v in row[key])
    if row['pmid']:
        curie = pmid_curie(row['pmid'])
        if curie is None:
            raise ValueError(f'Not a PubMed identifier: {row["pmid"]!r}')
        annotations.append(Annotation(term=slots.publications, value=curie))
    return Relation(
        subject=gene,
        predicate=slots.located_in,
        object=Entity(
            type=OntologyClass,
            identifiers=[Identifier(type=namespace, value=row['vesicle_type'])],
        ),
        annotations=annotations,
    )


def map_vesicle_type(row: dict[str, Any], *, namespace: str, description: str) -> Entity:
    """A vesicle type of the database."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=namespace, value=row['name']),
            Identifier(type=Namespace.NAME, value=row['name']),
        ],
        annotations=[Annotation(term=slots.description, value=description)],
    )
