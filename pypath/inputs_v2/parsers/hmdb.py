"""
HMDB (Human Metabolome Database) XML parser.

Parses HMDB's metabolites XML format and yields flat dictionaries.
"""

from __future__ import annotations

from collections.abc import Generator

import lxml.etree as etree


def _raw(opener, max_records: int | None = None, **_kwargs: object) -> Generator[dict[str, str], None, None]:
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
