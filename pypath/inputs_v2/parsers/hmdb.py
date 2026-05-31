"""
HMDB (Human Metabolome Database) XML parser.

Parses HMDB's metabolites XML format and yields flat dictionaries.
"""

from __future__ import annotations

from collections.abc import Generator
import re

import lxml.etree as etree

from pypath.inputs_v2.parsers.obo import iter_obo


def _raw(
    opener,
    max_records: int | None = None,
    chemont_name_map: dict[str, list[str]] | None = None,
    **_kwargs: object,
) -> Generator[dict[str, str], None, None]:
    """
    Parse HMDB metabolites XML file.

    Args:
        opener: File opener from download_and_open
        max_records: Optional limit on number of records to process

    Yields:
        Dictionary for each metabolite with flattened fields
    """
    if not opener or not opener.result:
        return

    xml_file = None
    for filename, file_handle in opener.result.items():
        if 'metabolites.xml' in filename.lower():
            xml_file = file_handle
            break

    if not xml_file:
        return

    xmlns = '{http://www.hmdb.ca}'

    # XXX: Testing if get_from_path can substitute all use-cases

    #def get_text(elem: etree._Element, tag: str, parent: etree._Element = None) -> str:
    #    parent = parent if parent is not None else elem
    #    child = parent.find(f'{xmlns}{tag}')
    #    return child.text.strip() if child is not None and child.text else ''

    #def get_list(elem: etree._Element, container_tag: str, item_tag: str) -> str:
    #    container = elem.find(f'{xmlns}{container_tag}')
    #    if container is None:
    #        return ''
    #    items = [
    #        item.text.strip() for item in container.findall(f'{xmlns}{item_tag}')
    #        if item.text and item.text.strip()
    #    ]
    #    return ';'.join(items)

    def get_from_path(elem, *args):

        path = '/'.join([f'{xmlns}{a}' for a in args])
        res = [e.text.strip() for e in elem.findall(path) if e.text is not None]

        return ';'.join(res)

    context = etree.iterparse(xml_file, events=('end',), tag=f'{xmlns}metabolite')

    count = 0
    for _, elem in context:
        if max_records is not None and count >= max_records:
            break

        #pubmed_ids = []
        #general_refs = elem.find(f'{xmlns}general_references')
        #if general_refs is not None:
        #    for ref in general_refs.findall(f'{xmlns}reference'):
        #        pubmed_id = ref.find(f'{xmlns}pubmed_id')
        #        if pubmed_id is not None and pubmed_id.text:
        #            pubmed_ids.append(pubmed_id.text.strip())

        taxonomy_terms = _taxonomy_terms(elem, xmlns)
        taxonomy_association_terms = _taxonomy_association_terms(elem, xmlns)

        yield {
            #'accession': get_text(elem, 'accession'),
            #'name': get_text(elem, 'name'),
            #'traditional_iupac': get_text(elem, 'traditional_iupac'),
            #'iupac_name': get_text(elem, 'iupac_name'),
            #'synonyms': get_list(elem, 'synonyms', 'synonym'),
            #'inchikey': get_text(elem, 'inchikey'),
            #'inchi': get_text(elem, 'inchi'),
            #'smiles': get_text(elem, 'smiles'),
            #'chebi_id': get_text(elem, 'chebi_id'),
            #'pubchem_compound_id': get_text(elem, 'pubchem_compound_id'),
            #'kegg_id': get_text(elem, 'kegg_id'),
            #'drugbank_id': get_text(elem, 'drugbank_id'),
            #'cas_registry_number': get_text(elem, 'cas_registry_number'),
            #'description': get_text(elem, 'description'),
            #'pubmed_ids': ';'.join(pubmed_ids) if pubmed_ids else '',
            'accession': get_from_path(elem, 'accession'),
            'name': get_from_path(elem, 'name'),
            'traditional_iupac': get_from_path(elem, 'traditional_iupac'),
            'iupac_name': get_from_path(elem, 'iupac_name'),
            'synonyms': get_from_path(elem, 'synonyms', 'synonym'),
            'inchikey': get_from_path(elem, 'inchikey'),
            'inchi': get_from_path(elem, 'inchi'),
            'smiles': get_from_path(elem, 'smiles'),
            'chebi_id': get_from_path(elem, 'chebi_id'),
            'pubchem_compound_id': get_from_path(elem, 'pubchem_compound_id'),
            'kegg_id': get_from_path(elem, 'kegg_id'),
            'drugbank_id': get_from_path(elem, 'drugbank_id'),
            'cas_registry_number': get_from_path(elem, 'cas_registry_number'),
            'description': get_from_path(elem, 'description'),
            'pubmed_ids': get_from_path(elem, 'general_references', 'reference', 'pubmed_id'),
            'diseases': get_from_path(elem, 'diseases', 'disease', 'name'),
            'cellular_locations': get_from_path(elem, 'biological_properties', 'cellular_locations', 'cellular'),
            'tissue_locations': get_from_path(elem, 'biological_properties', 'tissue_locations', 'tissue'),
            'biospecimen_locations': get_from_path(elem, 'biological_properties', 'biospecimen_locations', 'biospecimen'),
            'chemont_ids': _mapped_ids(taxonomy_association_terms, chemont_name_map),
            'chemont_annotation_ids': _mapped_ids(taxonomy_terms, chemont_name_map),
        }

        count += 1

        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    del context


def chemont_name_map(opener) -> dict[str, list[str]]:
    """Build a normalized ChemOnt label-to-accession map from the OBO file."""
    mapping: dict[str, list[str]] = {}
    for row in iter_obo(opener):
        _add_mapping(mapping, row.get('name'), row.get('id'))
    return mapping


def _taxonomy_terms(elem: etree._Element, xmlns: str) -> list[str]:
    taxonomy = elem.find(f'{xmlns}taxonomy')
    if taxonomy is None:
        return []

    values = []
    for tag in (
        'direct_parent',
        'kingdom',
        'super_class',
        'class',
        'sub_class',
        'molecular_framework',
    ):
        value = _child_text(taxonomy, xmlns, tag)
        if value:
            values.append(value)

    for container_tag, item_tag in (
        ('alternative_parents', 'alternative_parent'),
        ('substituents', 'substituent'),
    ):
        container = taxonomy.find(f'{xmlns}{container_tag}')
        if container is None:
            continue
        for item in container.findall(f'{xmlns}{item_tag}'):
            value = _clean_text(item.text)
            if value:
                values.append(value)

    return _dedupe(values)


def _taxonomy_association_terms(elem: etree._Element, xmlns: str) -> list[str]:
    taxonomy = elem.find(f'{xmlns}taxonomy')
    if taxonomy is None:
        return []

    values = []
    for tag in ('direct_parent', 'sub_class'):
        value = _child_text(taxonomy, xmlns, tag)
        if value:
            values.append(value)

    alternative_parents = taxonomy.find(f'{xmlns}alternative_parents')
    if alternative_parents is not None:
        for item in alternative_parents.findall(f'{xmlns}alternative_parent'):
            value = _clean_text(item.text)
            if value:
                values.append(value)

    return _dedupe(values)


def _mapped_ids(
    labels: list[str],
    mapping: dict[str, list[str]] | None,
) -> str:
    if not mapping:
        return ''

    ids = []
    for label in labels:
        ids.extend(mapping.get(_normalize_label(label), []))

    return ';'.join(_dedupe(ids))


def _add_mapping(
    mapping: dict[str, list[str]],
    label: object,
    accession: object,
) -> None:
    label_text = _clean_text(label)
    accession_text = _clean_text(accession)
    if not label_text or not accession_text:
        return
    key = _normalize_label(label_text)
    values = mapping.setdefault(key, [])
    if accession_text not in values:
        values.append(accession_text)


def _normalize_label(value: str) -> str:
    return re.sub(r'\s+', ' ', value.casefold()).strip()


def _child_text(elem: etree._Element, xmlns: str, tag: str) -> str:
    child = elem.find(f'{xmlns}{tag}')
    return _clean_text(child.text if child is not None else None)


def _clean_text(value: object) -> str:
    return ' '.join(str(value or '').split())


def _dedupe(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out
