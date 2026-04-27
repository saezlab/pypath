"""
WikiPathways RDF parser.
"""

from __future__ import annotations

import re
import urllib.request
from collections.abc import Generator, Iterable
from typing import Any

from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import DC, DCTERMS, FOAF, OWL, RDF, RDFS

from pypath.internals.cv_terms import EntityTypeCv


WP = Namespace('http://vocabularies.wikipathways.org/wp#')
PAV = Namespace('http://purl.org/pav/')
CITO = Namespace('http://purl.org/spar/cito/')

_CURRENT_RDF_URL = 'https://data.wikipathways.org/current/rdf/'
_RDF_FILENAME_RE = re.compile(r'href="(wikipathways-\d+-rdf-wp\.zip)"')
_WP_ID_RE = re.compile(r'(WP\d+(?:_r\d+)?)')
_NCBI_TAXON_RE = re.compile(r'NCBITaxon_(\d+)$')
_PUBMED_RE = re.compile(r'(\d+)$')
_INTERACTION_URI_RE = re.compile(r'/([^/]+)$')

_DATA_CACHE: dict[str, list[dict[str, str]]] = {}

_XREF_FIELDS = {
    'bdbUniprot': 'uniprot',
    'bdbEntrezGene': 'entrez',
    'bdbEnsembl': 'ensembl',
    'bdbChEBI': 'chebi',
    'bdbHmdb': 'hmdb',
    'bdbKeggCompound': 'kegg_compound',
    'bdbPubChem': 'pubchem_compound',
    'bdbWikidata': 'wikidata',
    'bdbHgncSymbol': 'hgnc',
}

_IDENTIFIER_PATTERNS = {
    'uniprot': re.compile(r'/uniprot/([^/?#]+)$'),
    'entrez': re.compile(r'/ncbigene/([^/?#]+)$'),
    'ensembl': re.compile(r'/ensembl/([^/?#]+)$'),
    'chebi': re.compile(r'/chebi/([^/?#]+)$'),
    'hmdb': re.compile(r'/hmdb/([^/?#]+)$'),
    'kegg_compound': re.compile(r'/kegg\.compound/([^/?#]+)$'),
    'pubchem_compound': re.compile(r'/CID(\d+)$'),
    'wikidata': re.compile(r'/entity/([^/?#]+)$'),
    'hgnc': re.compile(r'/hgnc\.symbol/([^/?#]+)$'),
}

_INTERACTION_TYPES = {
    'Interaction',
    'DirectedInteraction',
    'Stimulation',
    'Inhibition',
    'Conversion',
    'Catalysis',
    'Binding',
    'ComplexBinding',
    'TranscriptionTranslation',
}


def current_rdf_url(**_kwargs: Any) -> str:
    """
    Resolve the current WikiPathways RDF pathway zip URL.
    """

    with urllib.request.urlopen(_CURRENT_RDF_URL, timeout=30) as response:
        html = response.read().decode('utf-8', 'ignore')

    match = _RDF_FILENAME_RE.search(html)

    if not match:
        raise RuntimeError('Could not resolve current WikiPathways RDF archive URL.')

    return f'{_CURRENT_RDF_URL}{match.group(1)}'


def _raw(
        opener,
        data_type: str,
        force_refresh: bool = False,
        **_kwargs: Any,
    ) -> Generator[dict[str, str], None, None]:
    """
    Parse WikiPathways RDF and emit normalized raw records.
    """

    records = _load_records(opener, force_refresh=force_refresh)
    yield from records.get(data_type, [])


def _load_records(opener, force_refresh: bool = False) -> dict[str, list[dict[str, str]]]:

    if _DATA_CACHE and not force_refresh:
        return _DATA_CACHE

    records = {
        'pathways': [],
        'pathway_terms': [],
        'interactions': [],
    }

    if not opener or not opener.result:
        return records

    for name, handle in sorted(opener.result.items()):

        if not name.endswith('.ttl'):
            continue

        graph = Graph()
        data = handle.read()
        if isinstance(data, bytes):
            data = data.decode('utf-8', 'ignore')
        graph.parse(data=data, format='ttl')

        pathway = _extract_pathway_record(graph)

        if not pathway:
            continue

        records['pathways'].append(pathway)
        records['pathway_terms'].append(_extract_pathway_term_record(pathway))
        records['interactions'].extend(_extract_interaction_records(graph, pathway))

    _DATA_CACHE.clear()
    _DATA_CACHE.update(records)

    return records


def _extract_pathway_record(graph: Graph) -> dict[str, str] | None:

    pathway_uri = next(graph.subjects(RDF.type, WP.Pathway), None)

    if not pathway_uri:
        return None

    pathway_uri_str = str(pathway_uri)
    pathway_version_id = _extract_wp_id(pathway_uri_str)
    pathway_id = _extract_pathway_id(pathway_version_id)

    return {
        'pathway_uri': pathway_uri_str,
        'pathway_id': pathway_id,
        'pathway_version_id': pathway_version_id,
        'title': _first_literal(graph, pathway_uri, DC.title),
        'description': _first_literal(graph, pathway_uri, DCTERMS.description),
        'organism_name': _first_literal(graph, pathway_uri, WP.organismName),
        'taxon_id': _extract_taxon_id(_first_uri(graph, pathway_uri, WP.organism)),
        'pubmed_ids': _join_unique(_pubmed_ids(graph, pathway_uri)),
        'ontology_terms': _join_unique(_ontology_terms(graph, pathway_uri)),
    }


def _extract_pathway_term_record(pathway: dict[str, str]) -> dict[str, str]:
    xrefs = []
    if pathway.get('pathway_id'):
        xrefs.append(f"WikiPathways:{pathway['pathway_id']}")
    if pathway.get('pathway_version_id'):
        xrefs.append(f"WikiPathwaysVersion:{pathway['pathway_version_id']}")
    if pathway.get('taxon_id'):
        xrefs.append(f"NCBITaxon:{pathway['taxon_id']}")
    for pubmed_id in pathway.get('pubmed_ids', '').split(';'):
        if pubmed_id:
            xrefs.append(f"PMID:{pubmed_id}")

    return {
        'id': pathway.get('pathway_id', ''),
        'name': pathway.get('title', ''),
        'definition': pathway.get('description', ''),
        'synonyms': pathway.get('pathway_version_id', ''),
        'comments': '',
        'xrefs': _join_unique(xrefs),
    }



def _extract_interaction_records(
        graph: Graph,
        pathway: dict[str, str],
    ) -> list[dict[str, str]]:

    result = []
    seen_records = set()

    for interaction_uri in sorted(set(graph.subjects(RDF.type, WP.Interaction)), key=str):

        source_uri = _first_uri(graph, interaction_uri, WP.source)
        target_uri = _first_uri(graph, interaction_uri, WP.target)

        if not source_uri or not target_uri:
            continue

        sources = _terminal_entity_records(graph, source_uri)
        targets = _terminal_entity_records(graph, target_uri)

        if not sources or not targets:
            continue

        interaction_types = _interaction_types(graph, interaction_uri)
        interaction_uri_str = str(interaction_uri)

        for source in sources:
            for target in targets:

                if (
                    source['uri'] == target['uri']
                    or source['entity_type'] == EntityTypeCv.PATHWAY.value
                    or target['entity_type'] == EntityTypeCv.PATHWAY.value
                ):
                    continue

                record_key = (interaction_uri_str, source['uri'], target['uri'])
                if record_key in seen_records:
                    continue
                seen_records.add(record_key)

                result.append(
                    {
                        'interaction_uri': interaction_uri_str,
                        'interaction_local_id': _interaction_local_id(interaction_uri_str),
                        'interaction_types': _join_unique(interaction_types),
                        'pathway_id': pathway['pathway_id'],
                        'pathway_version_id': pathway['pathway_version_id'],
                        'pathway_term_accession': pathway['pathway_id'],
                        'pathway_ontology_terms': pathway['ontology_terms'],
                        'organism_name': pathway['organism_name'],
                        'taxon_id': pathway['taxon_id'],
                        **_prefix_record('source', source),
                        **_prefix_record('target', target),
                    }
                )

    return result


def _extract_entity_record(graph: Graph, entity_uri: str) -> dict[str, str]:

    subject = URIRef(entity_uri)
    rdf_types = {
        _local_name(obj)
        for obj in graph.objects(subject, RDF.type)
    }

    xrefs = {
        field: []
        for field in _XREF_FIELDS.values()
    }

    for predicate_suffix, field in _XREF_FIELDS.items():
        values = [
            _extract_identifier_value(field, str(obj))
            for obj in graph.objects(subject, WP[predicate_suffix])
        ]
        xrefs[field].extend(value for value in values if value)

    for field, value in _infer_identifier_from_uri(entity_uri).items():
        xrefs[field].append(value)

    identifier_uri = _first_uri(graph, subject, DC.identifier)
    if identifier_uri:
        for field, value in _infer_identifier_from_uri(identifier_uri).items():
            xrefs[field].append(value)

    return {
        'uri': entity_uri,
        'entity_type': _entity_type(rdf_types).value,
        'label': _first_literal(graph, subject, RDFS.label),
        'uniprot': _join_unique(xrefs['uniprot']),
        'entrez': _join_unique(xrefs['entrez']),
        'ensembl': _join_unique(xrefs['ensembl']),
        'chebi': _join_unique(xrefs['chebi']),
        'hmdb': _join_unique(xrefs['hmdb']),
        'kegg_compound': _join_unique(xrefs['kegg_compound']),
        'pubchem_compound': _join_unique(xrefs['pubchem_compound']),
        'wikidata': _join_unique(xrefs['wikidata']),
        'hgnc': _join_unique(xrefs['hgnc']),
    }


def _terminal_entity_records(
        graph: Graph,
        entity_uri: str,
        visited: set[str] | None = None,
    ) -> list[dict[str, str]]:
    """
    Resolve WikiPathways interaction references to biological participants.

    WikiPathways RDF may use another ``wp:Interaction`` node as the source or
    target of an interaction (line segments chained through anchors in GPML).
    Those nodes are graphical/relational constructs, not biological entities, so
    we recursively expand them to their terminal non-interaction endpoints.
    """

    visited = visited or set()

    if entity_uri in visited:
        return []

    visited.add(entity_uri)

    if not _is_interaction_node(graph, entity_uri):
        return [_extract_entity_record(graph, entity_uri)]

    subject = URIRef(entity_uri)
    records = []

    for predicate in (WP.source, WP.target):
        nested_uri = _first_uri(graph, subject, predicate)
        if nested_uri:
            records.extend(_terminal_entity_records(graph, nested_uri, visited.copy()))

    return _unique_entity_records(records)


def _unique_entity_records(records: Iterable[dict[str, str]]) -> list[dict[str, str]]:

    result = []
    seen = set()

    for record in records:
        uri = record.get('uri')
        if uri and uri not in seen:
            seen.add(uri)
            result.append(record)

    return result


def _is_interaction_node(graph: Graph, entity_uri: str) -> bool:

    subject = URIRef(entity_uri)
    rdf_types = {
        _local_name(obj)
        for obj in graph.objects(subject, RDF.type)
    }

    return bool(rdf_types & _INTERACTION_TYPES)


def _interaction_types(graph: Graph, interaction_uri: URIRef) -> list[str]:
    return [
        rdf_type
        for rdf_type in sorted({
            _local_name(obj)
            for obj in graph.objects(interaction_uri, RDF.type)
        })
        if rdf_type in _INTERACTION_TYPES
    ]


def _entity_type(rdf_types: set[str]) -> EntityTypeCv:

    if 'Complex' in rdf_types:
        return EntityTypeCv.COMPLEX

    if 'Metabolite' in rdf_types:
        return EntityTypeCv.SMALL_MOLECULE

    if rdf_types & {'Protein', 'GeneProduct'}:
        return EntityTypeCv.PROTEIN

    if 'Rna' in rdf_types:
        return EntityTypeCv.RNA

    if 'Dna' in rdf_types:
        return EntityTypeCv.DNA

    if 'Pathway' in rdf_types:
        return EntityTypeCv.PATHWAY

    return EntityTypeCv.PHYSICAL_ENTITY


def _ontology_terms(graph: Graph, pathway_uri: URIRef) -> list[str]:

    predicates = (WP.ontologyTag, WP.pathwayOntologyTag)
    values = []

    for predicate in predicates:
        values.extend(
            _ontology_term_from_uri(str(obj))
            for obj in graph.objects(pathway_uri, predicate)
        )

    return [value for value in values if value]


def _pubmed_ids(graph: Graph, pathway_uri: URIRef) -> list[str]:

    predicates = (DCTERMS.references, CITO.cites)
    values = []

    for predicate in predicates:
        values.extend(
            _extract_pubmed_id(str(obj))
            for obj in graph.objects(pathway_uri, predicate)
        )

    return [value for value in values if value]


def _prefix_record(prefix: str, record: dict[str, str]) -> dict[str, str]:
    return {f'{prefix}_{key}': value for key, value in record.items()}


def _first_literal(graph: Graph, subject: URIRef, predicate: URIRef) -> str:

    value = next(graph.objects(subject, predicate), None)
    return str(value) if value is not None else ''


def _first_uri(graph: Graph, subject: URIRef, predicate: URIRef) -> str:

    value = next(graph.objects(subject, predicate), None)
    return str(value) if value is not None else ''


def _extract_wp_id(value: str) -> str:

    match = _WP_ID_RE.search(value)
    return match.group(1) if match else ''


def _extract_pathway_id(pathway_version_id: str) -> str:
    return pathway_version_id.split('_', 1)[0] if pathway_version_id else ''


def _extract_taxon_id(value: str) -> str:

    match = _NCBI_TAXON_RE.search(value)
    return match.group(1) if match else ''


def _extract_pubmed_id(value: str) -> str:

    match = _PUBMED_RE.search(value)
    return match.group(1) if match else ''


def _ontology_term_from_uri(value: str) -> str:

    if value.startswith('http://purl.obolibrary.org/obo/'):
        return value.rsplit('/', maxsplit=1)[-1].replace('_', ':', 1)

    if '#' in value:
        return value.rsplit('#', maxsplit=1)[-1]

    return value.rsplit('/', maxsplit=1)[-1]


def _interaction_local_id(value: str) -> str:

    match = _INTERACTION_URI_RE.search(value)
    return match.group(1) if match else value


def _infer_identifier_from_uri(value: str) -> dict[str, str]:

    for field, pattern in _IDENTIFIER_PATTERNS.items():
        match = pattern.search(value)
        if match:
            return {field: match.group(1)}

    return {}


def _extract_identifier_value(field: str, value: str) -> str:

    match = _IDENTIFIER_PATTERNS[field].search(value)
    return match.group(1) if match else ''


def _local_name(value: URIRef) -> str:

    value_str = str(value)

    if '#' in value_str:
        return value_str.rsplit('#', maxsplit=1)[-1]

    return value_str.rsplit('/', maxsplit=1)[-1]


def _join_unique(values: Iterable[str]) -> str:

    seen = []

    for value in values:
        if value and value not in seen:
            seen.append(value)

    return ';'.join(seen)
