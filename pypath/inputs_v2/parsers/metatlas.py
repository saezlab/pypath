"""Raw parser for the Human-GEM YAML model.

The Human-GEM YAML (``model/Human-GEM.yml``) follows the standard-GEM format
from SysBioChalmers.  Its top-level keys are ``metabolites``, ``reactions``,
``genes``, ``compartments``, and ``metaData``.  Each section may be serialised
as a YAML ordered-map (``!!omap``), which PyYAML parses as a list of
``(key, value)`` pairs rather than a plain dict; every section is normalised
with ``dict()`` on first access.

Differences from the Recon3D BiGG JSON that require a separate parser:

- File format: YAML (with a broken-quote pre-processing step) vs JSON.
- Compartment encoding: Human-GEM metabolite IDs append the compartment letter
  directly (``MAM01039c``); the compartment is also stored in a dedicated
  ``compartment`` field, so stripping is done by reading that field.
- Gene identifiers: Ensembl ENSG IDs with no isoform suffix — no stripping.
- Cross-references: absent from the YAML; only available in separate TSV files.
- EC codes: stored in an ``eccodes`` field (list) rather than in an annotation
  sub-dict.
- Subsystem: may be a list in some reactions; joined to a ``'; '``-delimited
  string.
- The stoichiometry ``metabolites`` dict inside each reaction may be an omap.

The gene-rule boolean logic (``or`` / ``and``) is identical to Recon3D, so
``_parse_gene_rule`` is imported from the recon3d parser rather than
duplicated.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

import yaml

from pypath.inputs_v2.base import read_opener_text
from pypath.inputs_v2.parsers.recon3d import _parse_gene_rule

_RE_BROKEN_QUOTES = re.compile(
    r'^(\s*- \w+: )'
    r'("(?:[^"\\]|\\.)*")'
    r'(\S.*)$',
    re.MULTILINE,
)


def _fix_yaml_quoting(content: str) -> str:
    """
    Fix broken double-quote usage in the Human-GEM YAML.

    Some reactions contain metabolite names with unescaped double quotes
    that confuse the YAML parser.  This wraps such values in single quotes.

    Args:
        content: Raw YAML file content as a string.

    Returns:
        The content string with broken double-quoted scalars repaired.
    """

    def _fix_match(m):

        prefix = m.group(1)
        value = m.group(2) + m.group(3)
        escaped = value.replace("'", "''")

        return f"{prefix}'{escaped}'"

    return _RE_BROKEN_QUOTES.sub(_fix_match, content)


def _load(opener: Any) -> dict:
    """
    Parse the Human-GEM YAML from the opener.

    The opener already points at a locally cached file managed by the
    framework's cachedir; no additional caching is needed here.

    Args:
        opener: An opener object as provided by the inputs_v2 framework,
            pointing at the locally cached Human-GEM YAML file.

    Returns:
        The top-level YAML content as a plain dict with keys
        ``metabolites``, ``reactions``, ``genes``, ``compartments``,
        and ``metaData``.
    """

    content = read_opener_text(opener)
    content = _fix_yaml_quoting(content)
    data = yaml.safe_load(content)

    if isinstance(data, list):
        data = dict(data)

    return data


def _as_dict(obj) -> dict:
    """
    Normalise a YAML omap (list of pairs) or plain dict to a dict.

    Args:
        obj: Either a list of ``(key, value)`` pairs as produced by PyYAML
            for ``!!omap`` nodes, or an existing dict, or ``None``.

    Returns:
        A plain dict.  Returns an empty dict for falsy input.
    """

    if isinstance(obj, list):
        return dict(obj)

    return obj or {}


def _strip_compartment(met_id: str, compartment: str) -> str:
    """
    Remove the compartment letter from the end of a MAM metabolite ID.

    Human-GEM IDs append the compartment letter directly (e.g. ``MAM01039c``).
    The compartment code is stored in the ``compartment`` field so stripping
    is reliable.  All compartment codes are single letters.

    Args:
        met_id: Full metabolite ID including compartment suffix,
            e.g. ``'MAM01039c'``.
        compartment: Single-letter compartment code, e.g. ``'c'``.

    Returns:
        Base metabolite ID with the compartment letter removed,
        e.g. ``'MAM01039'``.
    """

    if compartment and met_id.endswith(compartment):
        return met_id[:-1]

    return met_id


def _metabolites(data: dict) -> Generator[dict, None, None]:
    """
    Yield one record per unique base metabolite ID.

    Compartment suffixes are stripped using the explicit ``compartment`` field
    and the results are deduplicated so that the same molecule in different
    compartments produces a single record.

    Args:
        data: The parsed top-level Human-GEM YAML dict.

    Yields:
        One dict per unique metabolite with keys ``human_gem_metabolite_id``,
        ``name``, ``formula``, and ``charge``.
    """

    seen: set[str] = set()

    for m in data.get('metabolites', []):

        m = _as_dict(m)
        compartment = m.get('compartment', '')
        raw_id = m.get('id', '')
        base_id = _strip_compartment(raw_id, compartment)

        if base_id in seen:
            continue

        seen.add(base_id)

        yield {
            'human_gem_metabolite_id': base_id,
            'name': m.get('name'),
            'formula': m.get('formula'),
            'charge': m.get('charge'),
        }


def _reactions(data: dict) -> Generator[dict, None, None]:
    """
    Yield one record per reaction with stoichiometry encoded as strings.

    Stoichiometry is packed as ``||``-delimited ``base_id:compartment:coef``
    strings for reactants and products — the same encoding used in Recon3D so
    the ``metatlas.py`` schema can reuse the same ``FieldConfig`` map lambdas.

    Direction follows the same convention as Recon3D: ``'reversible'`` when
    ``lower_bound < 0 < upper_bound``, ``'left_to_right'`` otherwise.

    Args:
        data: The parsed top-level Human-GEM YAML dict.

    Yields:
        One dict per reaction with keys ``human_gem_reaction_id``, ``name``,
        ``subsystem``, ``direction``, ``eccodes``, ``reactants``,
        and ``products``.
    """

    for r in data.get('reactions', []):

        r = _as_dict(r)
        mets = _as_dict(r.get('metabolites', {}))

        subsystem = r.get('subsystem')

        if isinstance(subsystem, list):
            subsystem = '; '.join(str(s) for s in subsystem)

        eccodes = r.get('eccodes')

        if isinstance(eccodes, list):
            eccodes = '; '.join(str(e) for e in eccodes)

        reactants = []
        products = []

        for met_id, stoich in mets.items():

            base_id, compartment = met_id[:-1], met_id[-1]

            if stoich < 0:
                reactants.append(f'{base_id}:{compartment}:{abs(stoich)}')
            else:
                products.append(f'{base_id}:{compartment}:{stoich}')

        lb = float(r['lower_bound'])
        ub = float(r['upper_bound'])
        direction = 'reversible' if lb < 0 < ub else 'left_to_right'

        yield {
            'human_gem_reaction_id': r.get('id', ''),
            'name': r.get('name'),
            'subsystem': subsystem,
            'direction': direction,
            'eccodes': eccodes,
            'reactants': '||'.join(reactants),
            'products': '||'.join(products),
        }


_TRANSPORT_SUBSYSTEM = 'Transport reactions'


def _transport_reactions(data: dict) -> Generator[dict, None, None]:
    """
    Yield reaction records for transport reactions only.

    Transport reactions are identified by ``subsystem == 'Transport reactions'``
    in the Human-GEM YAML, consistent with the expert-curated annotation.

    Args:
        data: The parsed top-level Human-GEM YAML dict.

    Yields:
        One dict per transport reaction; same keys as :func:`_reactions`.
    """

    for record in _reactions(data):

        if record.get('subsystem') == _TRANSPORT_SUBSYSTEM:
            yield record


def _metabolic_reactions(data: dict) -> Generator[dict, None, None]:
    """
    Yield reaction records for metabolic (non-transport) reactions only.

    Args:
        data: The parsed top-level Human-GEM YAML dict.

    Yields:
        One dict per metabolic reaction; same keys as :func:`_reactions`.
    """

    for record in _reactions(data):

        if record.get('subsystem') != _TRANSPORT_SUBSYSTEM:
            yield record


def _catalysis(data: dict) -> Generator[dict, None, None]:
    """
    Yield one record per unique (enzyme, reaction) catalysis association.

    Gene rules use Ensembl ENSG IDs with no isoform suffix.  OR groups give
    independent isoenzyme rows; AND groups give complex rows.  Gene symbols
    are looked up from the genes section and added as ``enzyme_name`` on
    protein rows (complexes have no single name).

    Args:
        data: The parsed top-level Human-GEM YAML dict.

    Yields:
        One dict per unique (enzyme, reaction) pair with keys
        ``enzyme_type``, ``enzyme_ensembl``, ``enzyme_name``,
        ``complex_subunits``, ``reaction_id``, ``reaction_name``,
        and ``subsystem``.
    """

    gene_names = {
        dict(g)['id']: dict(g).get('name')
        for g in data.get('genes', [])
    }

    seen: set[tuple] = set()

    for r in data.get('reactions', []):

        r = _as_dict(r)
        subsystem = (
            '; '.join(r['subsystem'])
            if isinstance(r.get('subsystem'), list)
            else r.get('subsystem')
        )

        for subunit_list in _parse_gene_rule(r.get('gene_reaction_rule', '')):

            if len(subunit_list) == 1:

                ensembl_id = subunit_list[0]
                key = ('protein', ensembl_id, r.get('id', ''))

                if key in seen:
                    continue

                seen.add(key)

                yield {
                    'enzyme_type': 'protein',
                    'enzyme_ensembl': ensembl_id,
                    'enzyme_name': gene_names.get(ensembl_id),
                    'complex_subunits': None,
                    'reaction_id': r.get('id', ''),
                    'reaction_name': r.get('name'),
                    'subsystem': subsystem,
                }

            else:

                key = ('complex', tuple(sorted(subunit_list)), r.get('id', ''))

                if key in seen:
                    continue

                seen.add(key)

                yield {
                    'enzyme_type': 'complex',
                    'enzyme_ensembl': None,
                    'enzyme_name': None,
                    'complex_subunits': '||'.join(subunit_list),
                    'reaction_id': r.get('id', ''),
                    'reaction_name': r.get('name'),
                    'subsystem': subsystem,
                }


def _enzyme_complexes(data: dict) -> Generator[dict, None, None]:
    """
    Yield one record per unique multi-subunit enzyme complex.

    Deduplicates by sorted subunit composition across all reactions.

    Args:
        data: The parsed top-level Human-GEM YAML dict.

    Yields:
        One dict per unique complex with key ``complex_subunits``, a
        ``||``-delimited string of Ensembl IDs.
    """

    seen: set[tuple] = set()

    for r in data.get('reactions', []):

        r = _as_dict(r)

        for subunit_list in _parse_gene_rule(r.get('gene_reaction_rule', '')):

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
    'metabolites': _metabolites,
    'reactions': _reactions,
    'transport_reactions': _transport_reactions,
    'metabolic_reactions': _metabolic_reactions,
    'catalysis': _catalysis,
    'enzyme_complexes': _enzyme_complexes,
}


def _raw(opener: Any, data_type: str, **kwargs) -> Generator[dict, None, None]:
    """
    Entry point called by each Dataset.raw_parser lambda in metatlas.py.

    Loads the Human-GEM YAML and dispatches to the appropriate parser
    function based on ``data_type``.

    Args:
        opener: An opener object as provided by the inputs_v2 framework.
        data_type: One of ``'metabolites'``, ``'reactions'``,
            ``'catalysis'``, or ``'enzyme_complexes'``.
        **kwargs: Accepted but unused; present for compatibility with the
            inputs_v2 raw parser signature.

    Yields:
        One record dict per entity; keys differ by ``data_type``.
    """

    data = _load(opener)

    yield from _PARSERS[data_type](data)
