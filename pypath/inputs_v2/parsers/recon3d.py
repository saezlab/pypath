"""Raw parser for Recon3D BiGG JSON.

Each data type is parsed independently from the downloaded JSON file via a
dedicated ``_parse_*`` function.  ``_raw()`` loads the file and dispatches
to the appropriate function based on ``data_type``.
"""

from __future__ import annotations

import json
import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.base import _first_handle

_ISOFORM_RE = re.compile(r'_[A-Z]+\d*$')


def _annotation_list(annotation_dict: dict, key: str) -> list[str] | None:
    """Extract annotation values for a given key.

    Looks up ``key`` in ``annotation_dict`` and strips any leading URL path
    components (e.g. ``https://hmdb.ca/metabolites/HMDB0000122`` →
    ``HMDB0000122``).

    Args:
        annotation_dict: The ``annotation`` sub-dict from a BiGG metabolite
            or reaction record.
        key: The annotation key to look up (e.g. ``'hmdb'``, ``'chebi'``).

    Returns:
        A list of cleaned annotation value strings, or ``None`` if the key
        is absent or the value list is empty.
    """
    values = annotation_dict.get(key, [])
    if isinstance(values, str):
        values = [values]
    cleaned = [v.split('/')[-1] if '/' in v else v for v in values]
    return cleaned or None


def _strip_compartment(met_id: str) -> tuple[str, str]:
    """Split a BiGG metabolite ID into base ID and compartment code.

    BiGG metabolite IDs encode the compartment as a single-letter suffix
    separated by an underscore (e.g. ``glc__D_c`` → cytosol ``'c'``,
    ``glc__D_e`` → extracellular ``'e'``).

    Args:
        met_id: A BiGG metabolite ID, optionally with a compartment suffix
            (e.g. ``'atp_c'``, ``'glc__D_e'``).

    Returns:
        A tuple of ``(base_id, compartment)`` where ``compartment`` is the
        single-letter code, or an empty string if no compartment suffix is
        present.
    """
    parts = met_id.rsplit('_', 1)
    if len(parts) == 2 and len(parts[1]) == 1 and parts[1].isalpha():
        return parts[0], parts[1]
    return met_id, ''


def _strip_isoform(gene_id: str) -> str:
    """Strip the Recon3D isoform suffix from a gene ID.

    Recon3D appends isoform tags of the form ``_AT1``, ``_AT2``, etc. to
    gene identifiers.  These are model-internal labels with no external
    database equivalent, so they are removed to recover the base Entrez ID.

    Args:
        gene_id: A raw gene identifier from the BiGG JSON, e.g.
            ``'1591_AT1'`` or ``'2645'``.

    Returns:
        The Entrez gene ID with any isoform suffix removed.
    """
    return _ISOFORM_RE.sub('', gene_id)


def _parse_gene_rule(rule: str) -> list[list[str]]:
    """Parse a Boolean gene reaction rule into structured enzyme groups.

    Recon3D gene reaction rules use Boolean logic where ``or`` separates
    alternative isoenzymes and ``and`` joins the subunits of a multi-protein
    complex.  Parentheses are removed before splitting.  Word-boundary
    matching is used for ``or``/``and`` to avoid false splits on gene IDs
    containing those substrings.

    Recon3D-internal isoform suffixes of the form ``_AT1``, ``_AT2``, etc.
    are stripped so that gene identifiers match standard Entrez IDs.

    Args:
        rule: The raw gene reaction rule string from the BiGG JSON, e.g.
            ``'2645_AT1 or (3098_AT1 and 3099_AT2)'``.  An empty string or
            ``'nan'`` is treated as no association.

    Returns:
        A list of enzyme groups.  Each group is a list of Entrez gene ID
        strings: a single-element list represents a monomeric enzyme, and a
        multi-element list represents a protein complex (AND group).
        Returns an empty list if the rule is absent or empty.
    """
    if not rule or not rule.strip() or rule.strip() == '0':
        return []
    rule = re.sub(r'[()]', '', rule).strip()
    result = []
    for or_part in re.split(r'\bor\b', rule, flags=re.IGNORECASE):
        subunits = [g.strip() for g in re.split(r'\band\b', or_part, flags=re.IGNORECASE)]
        cleaned = [_strip_isoform(g) for g in subunits if g]
        cleaned = [g for g in cleaned if g]
        if cleaned:
            result.append(cleaned)
    return result


def _parse_metabolites(data: dict) -> Generator[dict, None, None]:
    """Yield one record per unique metabolite base ID.

    Iterates over the ``metabolites`` array in the BiGG JSON, strips
    compartment suffixes, and deduplicates so that the same molecule in
    different compartments produces only one record.

    Args:
        data: The parsed top-level BiGG JSON dict.

    Yields:
        dict: One record per unique base metabolite ID with the keys
            ``bigg_metabolite_id``, ``name``, ``formula``, ``charge``,
            ``hmdb``, ``chebi``, ``kegg_compound``, and ``metanetx``.
            Absent annotation values are ``None`` and will be skipped by
            the framework.
    """
    seen: set[str] = set()
    for m in data.get('metabolites', []):
        base_id, _ = _strip_compartment(m['id'])
        if base_id in seen:
            continue
        seen.add(base_id)
        ann = m.get('annotation', {})
        yield {
            'bigg_metabolite_id': base_id,
            'name': m.get('name'),
            'formula': m.get('formula'),
            'charge': m.get('charge'),
            'hmdb': _annotation_list(ann, 'hmdb'),
            'chebi': _annotation_list(ann, 'chebi'),
            'kegg_compound': _annotation_list(ann, 'kegg.compound'),
            'metanetx': _annotation_list(ann, 'metanetx.chemical'),
        }


def _parse_reactions(data: dict) -> Generator[dict, None, None]:
    """Yield one record per reaction.

    Iterates over the ``reactions`` array in the BiGG JSON and encodes
    stoichiometry as ``||``-delimited ``base_id:compartment:coefficient``
    strings for reactants and products separately.  Direction is ``'reversible'`` when
    ``lower_bound < 0 < upper_bound`` (consistent with the legacy module),
    and ``'left_to_right'`` otherwise.

    Args:
        data: The parsed top-level BiGG JSON dict.

    Yields:
        dict: One record per reaction with the keys ``bigg_reaction_id``,
            ``name``, ``subsystem``, ``direction``, ``ec``, ``reactants``,
            and ``products``.  Absent values are ``None`` and will be
            skipped by the framework.
    """
    for r in data.get('reactions', []):
        ann = r.get('annotation', {})
        reactants = []
        products = []
        for met_id, stoich in r.get('metabolites', {}).items():
            base_id, compartment = _strip_compartment(met_id)
            if stoich < 0:
                reactants.append(f'{base_id}:{compartment}:{abs(stoich)}')
            else:
                products.append(f'{base_id}:{compartment}:{stoich}')
        lb = r.get('lower_bound', 0)
        ub = r.get('upper_bound', 0)
        direction = 'reversible' if lb < 0 < ub else 'left_to_right'
        yield {
            'bigg_reaction_id': r['id'],
            'name': r.get('name'),
            'subsystem': r.get('subsystem'),
            'direction': direction,
            'ec': _annotation_list(ann, 'ec-code'),
            'metanetx_reaction': _annotation_list(ann, 'metanetx.reaction'),
            'reactants': '||'.join(reactants),
            'products': '||'.join(products),
        }


def _parse_genes(data: dict) -> Generator[dict, None, None]:
    """Yield one record per unique Entrez gene ID.

    Strips Recon3D isoform suffixes (e.g. ``_AT1``) from all gene IDs,
    discards the placeholder gene ``'0'``, and deduplicates so that
    isoforms of the same gene produce only one record.

    Args:
        data: The parsed top-level BiGG JSON dict.

    Yields:
        dict: One record per unique gene with the keys ``entrez_id`` and
            ``name``.  Absent names are ``None`` and will be skipped by
            the framework.
    """
    seen: set[str] = set()
    for g in data.get('genes', []):
        entrez = _strip_isoform(g['id'])
        if not entrez or entrez == '0' or entrez in seen:
            continue
        seen.add(entrez)
        yield {
            'entrez_id': entrez,
            'name': g.get('name'),
        }


def _parse_catalysis(data: dict) -> Generator[dict, None, None]:
    """Yield one record per unique (enzyme, reaction) catalysis association.

    Parses the ``gene_reaction_rule`` for every reaction and emits one row
    per OR group (alternative isoenzyme or complex).  Single-gene groups
    produce a ``protein`` row; multi-gene AND groups produce a ``complex``
    row.  Isoform suffixes are stripped and duplicates arising from
    multiple isoforms of the same gene are suppressed.

    Args:
        data: The parsed top-level BiGG JSON dict.

    Yields:
        dict: One record per unique (enzyme, reaction) pair with the keys
            ``enzyme_type`` (``'protein'`` or ``'complex'``),
            ``enzyme_entrez`` (Entrez ID for proteins, ``None`` for
            complexes), ``reaction_bigg_id``, ``reaction_name``, and
            ``subsystem``.
    """
    seen: set[tuple] = set()
    for r in data.get('reactions', []):
        for subunit_list in _parse_gene_rule(r.get('gene_reaction_rule', '')):
            subunit_list = list(dict.fromkeys(
                _strip_isoform(g) for g in subunit_list
                if _strip_isoform(g) not in ('', '0') # what does it do?
            ))
            if not subunit_list:
                continue
            if len(subunit_list) == 1:
                key = ('protein', subunit_list[0], r['id'])
                if key in seen:
                    continue
                seen.add(key)
                yield {
                    'enzyme_type': 'protein',
                    'enzyme_entrez': subunit_list[0],
                    'reaction_bigg_id': r['id'],
                    'reaction_name': r.get('name'),
                    'subsystem': r.get('subsystem'),
                }
            else:
                key = ('complex', tuple(sorted(subunit_list)), r['id'])
                if key in seen:
                    continue
                seen.add(key)
                yield {
                    'enzyme_type': 'complex',
                    'enzyme_entrez': None,
                    'reaction_bigg_id': r['id'],
                    'reaction_name': r.get('name'),
                    'subsystem': r.get('subsystem'),
                }


def _parse_enzyme_complexes(data: dict) -> Generator[dict, None, None]:
    """Yield one record per unique multi-subunit enzyme complex.

    Collects all AND groups from gene reaction rules across all reactions,
    deduplicates by sorted subunit composition, and emits one record per
    unique complex.  Isoform suffixes are stripped before deduplication.

    Args:
        data: The parsed top-level BiGG JSON dict.

    Yields:
        dict: One record per unique complex with the key
            ``complex_subunits``, a ``||``-delimited string of Entrez IDs
            in original (non-sorted) order.
    """
    seen: set[tuple] = set()
    for r in data.get('reactions', []):
        for subunit_list in _parse_gene_rule(r.get('gene_reaction_rule', '')):
            subunit_list = list(dict.fromkeys(
                _strip_isoform(g) for g in subunit_list if g
            ))
            if len(subunit_list) < 2:
                continue
            key = tuple(sorted(subunit_list))
            if key in seen:
                continue
            seen.add(key)
            yield {
                'complex_subunits': '||'.join(subunit_list),
            }


_PARSERS = {
    'metabolites': _parse_metabolites,
    'reactions': _parse_reactions,
    'genes': _parse_genes,
    'catalysis': _parse_catalysis,
    'enzyme_complexes': _parse_enzyme_complexes,
}


def _raw(opener: Any, data_type: str, **kwargs) -> Generator[dict, None, None]:
    """Yield raw record dicts for the requested data type.

    Entry point called by each ``Dataset.raw_parser`` lambda in
    ``pypath/inputs_v2/recon3d.py``.  Loads the BiGG JSON from the cached
    file and dispatches to the appropriate ``_parse_*`` function.

    Args:
        opener: An ``Opener`` instance providing access to the downloaded
            Recon3D JSON file.
        data_type: One of ``'metabolites'``, ``'reactions'``, ``'genes'``,
            ``'catalysis'``, or ``'enzyme_complexes'``.
        **kwargs: Accepted but unused; present for compatibility with the
            inputs_v2 raw parser signature.

    Yields:
        dict: One record dict per entity.  The keys differ by data type;
            see the corresponding ``_parse_*`` function for field definitions.
    """
    data = json.load(_first_handle(opener))
    yield from _PARSERS[data_type](data)
