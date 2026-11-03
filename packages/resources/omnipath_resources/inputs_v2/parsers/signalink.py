"""
Parsers of the SignaLink 3 files.

SignaLink is a resource of signaling pathways. Its pathway layer (layer 0) comes as
two tab separated files with a header line.

The nodes file has a row for each protein or RNA:

    id  name  alt_accession  tax_id  pathways  aliases  topology
    3   Uniprot:P49765      taxid:9606  Receptor tyrosine kinase      Ligand

``name`` is the database and accession (``Uniprot:P49765``, ``RNACentral:URS...``),
``pathways`` the signaling pathways of the protein and ``topology`` its functions in
them (``Ligand``, ``Receptor``, ``Mediator``, ``Scaffold``, ``Co-factor``,
``Transcription factor``, ``Endocytosis related``); both separate several values by
``|``. ``alt_accession`` and ``aliases`` are empty.

The edges file has a row for each directed pair of nodes. The pair merges the
statements of the databases listed in ``source_db``:

    interactor_a_node_name  interactor_b_node_name  interaction_detection_method
    first_author  publication_ids  interaction_types  source_db
    interaction_identifiers  confidence_scores  layer

``interaction_types`` mixes two flags (``is_direct:true``, ``is_directed:true``) with
PSI-MI terms written ``MI:0624(stimulation)``. A pair can have terms of both signs.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.parsers.base import iter_tsv

NODE_COLUMNS = ('id', 'name', 'tax_id', 'pathways', 'topology')
EDGE_COLUMNS = (
    'interactor_a_node_name',
    'interactor_b_node_name',
    'interaction_detection_method',
    'publication_ids',
    'interaction_types',
    'source_db',
    'layer',
)

# The PSI-MI terms of the interaction types that state an effect of the first
# node on the second, as (direction, aspect) of the effect on the second node.
# Directions are ``increased`` and ``decreased``, None where the term has none;
# the aspect is the property that changes (``activity``, ``phosphorylation``...),
# None where the term names none. The names are the PSI-MI names; a term stated
# with a modification means that modification is made (the enzyme adds the group)
# or undone (``de-`` terms).
EFFECT_TERMS = {
    'MI:0624': ('increased', None),  # stimulation
    'MI:0623': ('decreased', None),  # inhibition
    'MI:2235': ('increased', None),  # up-regulates
    'MI:2236': ('increased', 'activity'),  # up-regulates activity
    'MI:2238': ('increased', 'expression'),  # up-regulates quantity by expression
    'MI:2239': ('increased', 'stability'),  # up-regulates quantity by stabilization
    'MI:2240': ('decreased', None),  # down-regulates
    'MI:2241': ('decreased', 'activity'),  # down-regulates activity
    'MI:2242': ('decreased', 'abundance'),  # down-regulates quantity
    'MI:2243': ('decreased', 'expression'),  # down-regulates quantity by repression
    'MI:2244': ('decreased', 'stability'),  # down-regulates quantity by destabilization
    'MI:2254': ('increased', None),  # chemical activation reaction
    'MI:2247': (None, 'expression'),  # transcriptional regulation
    'MI:0414': (None, None),  # enzymatic reaction
    'MI:0217': ('increased', 'phosphorylation'),  # phosphorylation reaction
    'MI:0203': ('decreased', 'phosphorylation'),  # dephosphorylation reaction
    'MI:0220': ('increased', 'ubiquitination'),  # ubiquitination reaction
    'MI:0204': ('decreased', 'ubiquitination'),  # deubiquitination reaction
    'MI:0192': ('increased', 'acetylation'),  # acetylation reaction
    'MI:0197': ('decreased', 'acetylation'),  # deacetylation reaction
    'MI:0213': ('increased', 'methylation'),  # methylation reaction
    'MI:0566': ('increased', 'sumoylation'),  # sumoylation reaction
    'MI:0557': ('increased', 'ADP-ribosylation'),  # adp ribosylation reaction
    # Cleavage is a reaction that happens, like the phosphorylation above and as in
    # the proteolytic processing of HPRD, SPIKE and the macrophage resource.
    'MI:0194': ('increased', 'cleavage'),  # cleavage reaction
    'MI:0570': ('increased', 'cleavage'),  # protein cleavage
    'MI:2256': (None, 'localization'),  # relocalization reaction
    'MI:2252': (None, 'guanyl_nucleotide_exchange'),  # guanine nucleotide exchange factor
}
# Interaction types that state no effect: binding or the interaction itself.
BINDING_TERMS = frozenset(
    {
        'MI:0190',  # interaction type
        'MI:0914',  # association
        'MI:0915',  # physical association
        'MI:0462',  # bind
        'MI:0408',  # disulfide bond
        'MI:0195',  # covalent binding
        'MI:1126',  # self interaction
        'MI:0401',  # biochemical
        # An obsolete residue modification, which InnateDB writes as the interaction
        # type of some records; it names no effect and no aspect. innatedb.py takes
        # it as an untyped interaction, too.
        'MI:0179',  # other modification
    }
)
_FLAGS = ('is_direct', 'is_directed')
_NODE_NAME = re.compile(r'(Uniprot|RNACentral):(\S+)')
_TAXON = re.compile(r'taxid:(\d+)')


def _values(cell: str | None) -> list[str]:
    return list(dict.fromkeys(v.strip() for v in (cell or '').split('|') if v.strip()))


def _nodes(opener) -> dict[str, dict[str, Any]]:
    nodes = {}
    for row in iter_tsv(opener, quotechar=None, strict=True, required=NODE_COLUMNS):
        name = (row['name'] or '').strip()
        match = _NODE_NAME.fullmatch(name)
        taxon = _TAXON.fullmatch((row['tax_id'] or '').strip())
        if not match or not taxon or name in nodes:
            raise ValueError(f'SignaLink node not understood: {row!r}')
        nodes[name] = {
            'name': name,
            'database': match.group(1),
            'accession': match.group(2),
            'taxon': taxon.group(1),
            'pathways': _values(row['pathways']),
            'functions': _values(row['topology']),
        }
    if not nodes:
        raise ValueError('SignaLink nodes file has no rows')
    return nodes


def iter_nodes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    The nodes that have a pathway or a function.

    Each dictionary has ``name`` (``Uniprot:P49765``), ``database`` (``Uniprot`` or
    ``RNACentral``), ``accession``, ``taxon`` (digits), ``pathways`` and
    ``functions`` (lists, as written). The other nodes of the file, most of them
    members of no pathway, carry no annotation and are left out.
    """
    for node in _nodes(opener).values():
        if node['pathways'] or node['functions']:
            yield node


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The distinct pathways, then the distinct functions of the nodes file."""
    pathways, functions = {}, {}
    for node in _nodes(opener).values():
        pathways.update(dict.fromkeys(node['pathways']))
        functions.update(dict.fromkeys(node['functions']))
    for term in pathways:
        yield {'kind': 'pathway', 'term': term}
    for term in functions:
        yield {'kind': 'function', 'term': term}


def _pmids(cell: str | None) -> list[str]:
    """
    The PubMed CURIEs of a ``pubmed:123|pubmed:456`` cell.

    A number that is not a PubMed identifier (the file has some of more than eight
    digits, which are fused numbers) is left out. Another database is an error.
    """
    result = []
    for token in _values(cell):
        database, _, number = token.partition(':')
        if database != 'pubmed':
            raise ValueError(f'SignaLink publication not understood: {token!r}')
        if curie := pmid_curie(number):
            result.append(curie)
    return list(dict.fromkeys(result))


def _effects(terms: list[str]) -> list[tuple[str | None, str | None]]:
    """
    The (direction, aspect) pairs of the effect terms of a pair of nodes.

    One pair per direction that a term states: terms of opposite directions are
    separate statements. A term without a direction goes with the pair of no
    direction when no term states one. The aspect is set when the terms of the
    direction (and those without a direction) agree on one; terms that name
    different aspects leave it unset.
    """
    terms = [EFFECT_TERMS[t] for t in terms if t in EFFECT_TERMS]
    if not terms:
        return []
    directions = list(dict.fromkeys(d for d, _ in terms if d))
    general = {a for d, a in terms if d is None and a}
    result = []
    for direction in directions or [None]:
        aspects = {a for d, a in terms if d == direction and a} | general
        result.append((direction, next(iter(aspects)) if len(aspects) == 1 else None))
    return result


def iter_interactions(
    opener, *, nodes_download, force_refresh: bool = False, **kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """
    One dictionary per effect of each directed pair of nodes in the edges file.

    A pair is one record, or several when its effect terms have opposite directions
    (``direction`` is then ``increased`` in one and ``decreased`` in the other).

    - ``a``, ``b``: the nodes, as dictionaries of ``iter_nodes`` (all the nodes).
    - ``direction``, ``aspect``: the effect of ``a`` on ``b``, or None.
    - ``affects``: whether the pair has an effect term; ``direction`` and ``aspect``
      may still be None.
    - ``is_direct``: whether any of the merged statements is direct.
    - ``interaction_types``: the PSI-MI terms of the pair. ``methods``: the detection
      methods. ``pmids``: the PubMed identifiers. ``sources``: the labels of
      ``source_db``.

    A node the nodes file lacks, a PSI-MI term this parser has no meaning for, a
    flag it does not know or a layer other than 0 is an error.
    """
    nodes = _nodes(nodes_download.open(force_refresh=force_refresh, **kwargs))
    count = 0
    for row in iter_tsv(opener, quotechar=None, strict=True, required=EDGE_COLUMNS):
        a = nodes.get((row['interactor_a_node_name'] or '').strip())
        b = nodes.get((row['interactor_b_node_name'] or '').strip())
        if a is None or b is None or row['layer'] != '0':
            raise ValueError(f'SignaLink edge not understood: {row!r}')
        flags, terms = {}, []
        for item in _values(row['interaction_types']):
            name, _, value = item.partition(':')
            if name in _FLAGS:
                if value not in ('true', 'false'):
                    raise ValueError(f'SignaLink flag not understood: {item!r}')
                flags.setdefault(name, set()).add(value == 'true')
                continue
            term = _mitab.mi_term(_mitab.parse_mitab_item(item)) if name == 'MI' else None
            if term not in EFFECT_TERMS and term not in BINDING_TERMS:
                raise ValueError(f'SignaLink interaction type not understood: {item!r}')
            terms.append(term)
        if flags.get('is_directed') != {True}:
            raise ValueError(f'SignaLink edge is not directed: {row!r}')
        effects = _effects(terms)
        for direction, aspect in effects or [(None, None)]:
            count += 1
            yield {
                'a': a,
                'b': b,
                'affects': bool(effects),
                'direction': direction,
                'aspect': aspect,
                'is_direct': True in flags.get('is_direct', set()),
                'interaction_types': list(dict.fromkeys(terms)),
                'methods': _mitab.mi_terms(row['interaction_detection_method']),
                'pmids': _pmids(row['publication_ids']),
                'sources': _values(row['source_db']),
            }
    if not count:
        raise ValueError('SignaLink edges file has no rows')
