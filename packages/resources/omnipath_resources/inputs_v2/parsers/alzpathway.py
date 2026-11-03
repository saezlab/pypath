"""
Reader for the complexes of a CellDesigner map (SBML with the CellDesigner extension).

Nothing in this file knows about AlzPathway: :func:`iter_celldesigner_complexes`
reads any CellDesigner 4 file. The AlzPathway-specific interpretation is in
``inputs_v2/alzpathway.py``.

In a CellDesigner file a complex is a species of class ``COMPLEX``. Its
components are the "included species" that name it in ``complexSpecies``; a
component can itself be a complex. A protein, gene or RNA species points to an
entry of the model's protein, gene or RNA list (``proteinReference`` and so
on), which is what makes two species the same molecule; other classes are
identified by their class and name. Cross-references are MIRIAM URNs in the
RDF annotation of the species and of the reactions.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any
import xml.etree.ElementTree as ET

from pypath.inputs_v2.parsers.base import read_xml_root

__all__ = ['iter_celldesigner_complexes', 'iter_complexes']

_CD = '{http://www.sbml.org/2001/ns/celldesigner}'
_RDF = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
_MIRIAM = re.compile(r'urn:miriam:([^:]+):(.+)')
# ``proteinReference``, ``geneReference``, ``rnaReference``, ``antisensernaReference``
_REFERENCE_SUFFIX = 'Reference'


def _miriam(element: ET.Element) -> dict[str, list[str]]:
    """MIRIAM cross-references of an element: collection -> identifiers."""
    found: dict[str, list[str]] = {}
    annotation = element.find('{*}annotation')
    if annotation is None:
        return found
    for resource in annotation.iter(_RDF + 'li'):
        match = _MIRIAM.fullmatch(resource.get(_RDF + 'resource') or '')
        if match and match[2] not in found.setdefault(match[1], []):
            found[match[1]].append(match[2])
    return found


def _identity(species: ET.Element) -> dict[str, Any]:
    """Class, molecule reference and state of a species or included species."""
    identity = species.find(f'.//{_CD}speciesIdentity')
    if identity is None:
        raise ValueError(
            f'CellDesigner species {species.get("id")!r} has no speciesIdentity'
        )
    reference = next(
        (
            (child.text or '').strip()
            for child in identity
            if child.tag.endswith(_REFERENCE_SUFFIX)
        ),
        None,
    )
    homodimer = identity.findtext(f'{_CD}state/{_CD}homodimer')
    return {
        'class': (identity.findtext(_CD + 'class') or '').strip(),
        'reference': reference or None,
        'homodimer': int(homodimer) if homodimer and homodimer.isdigit() else None,
        'modifications': sorted(
            {
                modification.get('state')
                for modification in identity.iter(_CD + 'modification')
                if modification.get('state')
            }
        ),
    }


def iter_celldesigner_complexes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per complex species of a CellDesigner file.

    Output fields:

    - ``species_id``, ``name``: of the complex species.
    - ``compartment``: name of the compartment of the species (its id if the
      compartment has no name).
    - ``producing_reactions``: the reactions that have the complex as a
      product, each with ``id``, ``type`` (CellDesigner ``reactionType``) and
      ``xrefs`` (MIRIAM collection -> identifiers, e.g. ``pubmed``).
    - ``members``: the components, in file order, each with ``species_id``,
      ``name``, ``class`` (``PROTEIN``, ``SIMPLE_MOLECULE``, ``COMPLEX`` ...),
      ``reference`` (id of the protein, gene or RNA entry), ``homodimer``
      (number of copies, None if not drawn as a multimer), ``modifications``
      (residue states such as ``phosphorylated``), ``xrefs`` and ``members``
      (the components of a component that is itself a complex).

    Components carry no annotation of their own in a CellDesigner file. The
    ``xrefs`` of a component are those of the species of the model that are
    the same molecule: the same protein, gene or RNA entry or, for the other
    classes, the same class and name. Nothing else is inferred.

    A document that is not a CellDesigner model raises ``ValueError``.
    """
    root = read_xml_root(opener)
    if root is None:
        return
    model = root.find('{*}model')
    extension = (
        model.find(f'{{*}}annotation/{_CD}extension') if model is not None else None
    )
    if not root.tag.endswith('}sbml') or extension is None:
        raise ValueError('Expected a CellDesigner SBML document')

    compartments = {
        compartment.get('id'): compartment.get('name') or compartment.get('id')
        for compartment in model.iterfind('{*}listOfCompartments/{*}compartment')
    }
    species = list(model.iterfind('{*}listOfSpecies/{*}species'))
    # Cross-references of each molecule, from every species that shows it.
    xrefs: dict[tuple[str, str], dict[str, list[str]]] = {}
    for element in species:
        identity = _identity(element)
        key = (identity['class'], identity['reference'] or element.get('name') or '')
        collected = xrefs.setdefault(key, {})
        for collection, identifiers in _miriam(element).items():
            known = collected.setdefault(collection, [])
            known.extend(i for i in identifiers if i not in known)

    components: dict[str, list[ET.Element]] = {}
    for element in extension.iterfind(f'{_CD}listOfIncludedSpecies/{_CD}species'):
        parent = (element.findtext(f'.//{_CD}complexSpecies') or '').strip()
        components.setdefault(parent, []).append(element)

    producing: dict[str, list[dict[str, Any]]] = {}
    for reaction in model.iterfind('{*}listOfReactions/{*}reaction'):
        record = {
            'id': reaction.get('id'),
            'type': (reaction.findtext(f'.//{_CD}reactionType') or '').strip(),
            'xrefs': _miriam(reaction),
        }
        for product in reaction.iterfind('{*}listOfProducts/{*}speciesReference'):
            producing.setdefault(product.get('species'), []).append(record)

    def members(species_id: str, seen: tuple[str, ...]) -> list[dict[str, Any]]:
        if species_id in seen:
            raise ValueError(f'CellDesigner complex {species_id!r} contains itself')
        result = []
        for element in components.get(species_id, ()):
            identity = _identity(element)
            name = element.get('name') or ''
            result.append(
                {
                    'species_id': element.get('id'),
                    'name': name,
                    **identity,
                    'xrefs': {
                        collection: list(identifiers)
                        for collection, identifiers in xrefs.get(
                            (identity['class'], identity['reference'] or name), {}
                        ).items()
                    },
                    'members': members(element.get('id'), (*seen, species_id)),
                }
            )
        return result

    for element in species:
        if _identity(element)['class'] != 'COMPLEX':
            continue
        yield {
            'species_id': element.get('id'),
            'name': element.get('name') or '',
            'compartment': compartments.get(element.get('compartment')),
            'producing_reactions': producing.get(element.get('id'), []),
            'members': members(element.get('id'), ()),
        }


def iter_complexes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the complex species of the map; a map without complexes is an error."""
    count = 0
    for row in iter_celldesigner_complexes(opener):
        count += 1
        yield row
    if not count:
        raise ValueError('AlzPathway map has no complex species')
