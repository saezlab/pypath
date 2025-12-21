"""
Parse HMDB (Human Metabolome Database) data and emit Entity records.

This module downloads HMDB metabolite XML data and converts it into Entity
records using the declarative schema pattern from tabular_builder.
"""

from __future__ import annotations

from collections.abc import Generator

import lxml.etree as etree

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig


def _iterate_metabolites(opener, max_records: int | None = None) -> Generator[dict[str, str], None, None]:
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

    def get_text(elem: etree._Element, tag: str, parent: etree._Element = None) -> str:
        parent = parent if parent is not None else elem
        child = parent.find(f'{xmlns}{tag}')
        return child.text.strip() if child is not None and child.text else ''

    def get_list(elem: etree._Element, container_tag: str, item_tag: str) -> str:
        container = elem.find(f'{xmlns}{container_tag}')
        if container is None:
            return ''
        items = [
            item.text.strip() for item in container.findall(f'{xmlns}{item_tag}')
            if item.text and item.text.strip()
        ]
        return ';'.join(items)

    context = etree.iterparse(xml_file, events=('end',), tag=f'{xmlns}metabolite')

    count = 0
    for _, elem in context:
        if max_records is not None and count >= max_records:
            break

        pubmed_ids = []
        general_refs = elem.find(f'{xmlns}general_references')
        if general_refs is not None:
            for ref in general_refs.findall(f'{xmlns}reference'):
                pubmed_id = ref.find(f'{xmlns}pubmed_id')
                if pubmed_id is not None and pubmed_id.text:
                    pubmed_ids.append(pubmed_id.text.strip())

        yield {
            'accession': get_text(elem, 'accession'),
            'name': get_text(elem, 'name'),
            'traditional_iupac': get_text(elem, 'traditional_iupac'),
            'iupac_name': get_text(elem, 'iupac_name'),
            'synonyms': get_list(elem, 'synonyms', 'synonym'),
            'inchikey': get_text(elem, 'inchikey'),
            'inchi': get_text(elem, 'inchi'),
            'smiles': get_text(elem, 'smiles'),
            'chebi_id': get_text(elem, 'chebi_id'),
            'pubchem_compound_id': get_text(elem, 'pubchem_compound_id'),
            'kegg_id': get_text(elem, 'kegg_id'),
            'drugbank_id': get_text(elem, 'drugbank_id'),
            'cas_registry_number': get_text(elem, 'cas_registry_number'),
            'description': get_text(elem, 'description'),
            'pubmed_ids': ';'.join(pubmed_ids) if pubmed_ids else '',
        }

        count += 1

        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    del context


def _hmdb_raw(opener, max_records: int | None = None, **_kwargs: object):
    yield from _iterate_metabolites(opener, max_records=max_records)


config = ResourceConfig(
    id=ResourceCv.HMDB,
    name='Human Metabolome Database',
    url='https://hmdb.ca/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37953221',
    description=(
        'The Human Metabolome Database (HMDB) is a comprehensive database '
        'containing detailed information about small molecule metabolites '
        'found in the human body. It includes chemical, clinical, and '
        'biochemical/molecular biology data, with over 220,000 metabolite '
        'entries including both water-soluble and lipid-soluble metabolites.'
    ),
)

f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
    transform={
        'chebi': lambda v: f'CHEBI:{v}' if v else None,
    },
)

metabolites_schema = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.HMDB, value=f('accession')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.IUPAC_TRADITIONAL_NAME, value=f('traditional_iupac')),
        CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=f('iupac_name')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms', delimiter=';')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('inchikey')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('inchi')),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('smiles')),
        CV(
            term=IdentifierNamespaceCv.CHEBI,
            value=f('chebi_id', extract='chebi'),
        ),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('pubchem_compound_id')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('kegg_id')),
        CV(term=IdentifierNamespaceCv.DRUGBANK, value=f('drugbank_id')),
        CV(term=IdentifierNamespaceCv.CAS, value=f('cas_registry_number')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('description')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_ids', delimiter=';')),
    ),
)

resource = Resource(
    config,
    metabolites=Dataset(
        download=Download(
            url='https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
            filename='hmdb_metabolites.zip',
            subfolder='hmdb',
            large=True,
            ext='zip',
            default_mode='rb',
        ),
        mapper=metabolites_schema,
        raw_parser=_hmdb_raw,
    ),
)

metabolites = resource.metabolites
