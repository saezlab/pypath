#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
Standard GEM file parsing functions.

Standard GEMs may contain annotation files in TSV format:
- reactions.tsv - Reaction ID mappings
- metabolites.tsv - Metabolite ID mappings
- genes.tsv - Gene annotations

Note: Column names vary between GEMs. SysBioChalmers GEMs use names like
'rxnKEGGID', 'rxnBiGGID', while others use 'kegg.reaction', 'bigg.reaction'.
The TSV functions return raw dictionaries to accommodate this variation.

The YAML model file (model/{GEM}.yml) is present in all standard-GEM repos
and contains stoichiometry, flux bounds, gene rules, compartments, and
cross-references. The YAML functions return typed GemReaction and
GemMetabolite named tuples.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import io
import re

import yaml

from ._common import _log
from ._git import git_raw_file, gem_file_path, _parse_gem_index
from ._records import GemReaction, GemMetabolite, GemInteraction

__all__ = [
    'metatlas_gem_reactions',
    'metatlas_gem_metabolites',
    'metatlas_gem_genes',
    'metatlas_gem_tsv',
    'metatlas_gem_yaml',
    'metatlas_gem_yaml_reactions',
    'metatlas_gem_yaml_metabolites',
    'metatlas_gem_network',
    'metatlas_gem_transport_ids',
    'metatlas_gem_detect_gene_id_type',
    'metatlas_gem_transport_network',
]


def _get_gem_info(gem: str) -> tuple[str, str] | None:
    """
    Retrieves git host and repo for a GEM from the index.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').

    Returns:
        Tuple of (host, repo) or None if not found.
    """

    index = _parse_gem_index()

    for host, repos in index.items():
        for repo in repos:
            if repo.split('/')[-1] == gem:
                return host, repo

    _log(f'GEM {gem} not found in index.')
    return None


def metatlas_gem_tsv(
        gem: str,
        file: str,
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Downloads and parses a TSV file from a standard GEM repository.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        file: Relative path to the TSV file (e.g., 'model/reactions.tsv').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries for each row in the TSV file.
        Column names vary between GEMs.
    """

    gem_info = _get_gem_info(gem)

    if gem_info is None:
        return

    host, repo = gem_info

    _log(f'Downloading {file} from {gem}.')

    content = git_raw_file(host, repo, ref, file)

    if content is None:
        _log(f'File `{file}` not available in {gem}; '
             'not all GEMs provide TSV annotations.')
        return

    reader = csv.DictReader(io.StringIO(content), delimiter='\t')

    for row in reader:
        yield row


def metatlas_gem_reactions(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Parses reaction annotations from a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries with reaction ID mappings.
        Column names vary between GEMs (e.g., 'rxnKEGGID' vs 'kegg.reaction').
    """

    _log(f'Parsing reactions from {gem}.')

    yield from metatlas_gem_tsv(gem, 'model/reactions.tsv', ref)


def metatlas_gem_metabolites(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Parses metabolite annotations from a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries with metabolite ID mappings.
        Column names vary between GEMs.
    """

    _log(f'Parsing metabolites from {gem}.')

    yield from metatlas_gem_tsv(gem, 'model/metabolites.tsv', ref)


def metatlas_gem_genes(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Parses gene annotations from a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries with gene annotations.
        Column names vary between GEMs.
    """

    _log(f'Parsing genes from {gem}.')

    yield from metatlas_gem_tsv(gem, 'model/genes.tsv', ref)


_RE_BROKEN_QUOTES = re.compile(
    r'^(\s*- \w+: )'        # YAML omap key prefix, e.g. "      - name: "
    r'("(?:[^"\\]|\\.)*")'  # a double-quoted scalar, e.g. '"5"'
    r'(\S.*)$',             # trailing text that shouldn't be there
    re.MULTILINE,
)


def _fix_yaml_quoting(content: str) -> str:
    """
    Fix broken double-quote usage in GEM YAML files.

    Some GEM YAML files contain unescaped double quotes in string values,
    e.g. ``- name: "5"-deoxyadenosine...`` where ``"5"`` is parsed by YAML
    as a complete quoted scalar.  This wraps such values in single quotes.
    """

    def _fix_match(m):
        prefix = m.group(1)
        value = m.group(2) + m.group(3)
        escaped = value.replace("'", "''")
        return f"{prefix}'{escaped}'"

    return _RE_BROKEN_QUOTES.sub(_fix_match, content)


def metatlas_gem_yaml(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> dict | None:
    """
    Downloads and parses the YAML model file from a standard GEM repository.

    The YAML file (model/{gem}.yml) is present in all standard-GEM repos and
    contains the full model: reactions with stoichiometry and gene rules,
    metabolites with compartment annotations, and compartment definitions.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Returns:
        Dictionary with keys 'reactions', 'metabolites', 'compartments',
        'genes', and 'metaData'; or None if download failed.
    """

    gem_info = _get_gem_info(gem)

    if gem_info is None:
        return None

    host, repo = gem_info
    path = gem_file_path(gem, 'yml')

    _log(f'Downloading YAML model from {gem}.')

    content = git_raw_file(host, repo, ref, path)

    if content is None:
        _log(f'YAML model file not found for {gem} at {path}.')
        return None

    _log(f'Parsing YAML model for {gem}.')

    content = _fix_yaml_quoting(content)
    data = yaml.safe_load(content)

    # yaml.safe_load with !!omap returns list of (key, value) tuples
    if isinstance(data, list):
        data = dict(data)

    return data


def metatlas_gem_yaml_reactions(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> Generator[GemReaction, None, None]:
    """
    Parses reactions from the YAML model of a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        GemReaction named tuples with stoichiometry, bounds, and gene rules.
    """

    data = metatlas_gem_yaml(gem, ref)

    if data is None:
        return

    reactions = data.get('reactions', [])

    _log(f'Processing {len(reactions)} reactions from {gem} YAML.')

    for rxn in reactions:

        if isinstance(rxn, list):
            rxn = dict(rxn)

        mets = rxn.get('metabolites', {})

        if isinstance(mets, list):
            mets = dict(mets)

        subsystem = rxn.get('subsystem', '')

        if isinstance(subsystem, list):
            subsystem = '; '.join(str(s) for s in subsystem)

        eccodes = rxn.get('eccodes', '')

        if isinstance(eccodes, list):
            eccodes = '; '.join(str(e) for e in eccodes)

        yield GemReaction(
            id=rxn.get('id', ''),
            name=rxn.get('name', ''),
            metabolites=mets,
            lower_bound=float(rxn.get('lower_bound', 0)),
            upper_bound=float(rxn.get('upper_bound', 0)),
            gene_reaction_rule=rxn.get('gene_reaction_rule', ''),
            subsystem=subsystem,
            eccodes=eccodes,
        )


def metatlas_gem_yaml_metabolites(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> Generator[GemMetabolite, None, None]:
    """
    Parses metabolites from the YAML model of a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        GemMetabolite named tuples with compartment, formula, and charge.
    """

    data = metatlas_gem_yaml(gem, ref)

    if data is None:
        return

    metabolites = data.get('metabolites', [])

    _log(f'Processing {len(metabolites)} metabolites from {gem} YAML.')

    for met in metabolites:

        if isinstance(met, list):
            met = dict(met)

        charge = met.get('charge')

        if charge is not None:
            charge = int(charge)

        yield GemMetabolite(
            id=met.get('id', ''),
            name=met.get('name', ''),
            compartment=met.get('compartment', ''),
            formula=met.get('formula', ''),
            charge=charge,
        )


def _parse_gene_rule(rule: str) -> list[str]:
    """
    Parse a gene_reaction_rule string into enzyme identifiers.

    The rule is a Boolean expression with ``or`` (isoenzymes) and ``and``
    (complex subunits).  Parentheses are stripped — nested expressions are
    flattened, matching the OmnipathR approach.

    AND-connected genes (complex subunits) are joined into a single
    identifier by ``_``, with subunits sorted for consistency.

    Returns:
        List of enzyme identifiers (one per OR-alternative).  Single genes
        remain as-is; complexes become ``GENE1_GENE2`` strings.
        Empty list for orphan reactions (no gene rule).
    """

    if not rule or not rule.strip():
        return []

    rule = re.sub(r'[()]', '', rule)

    enzymes = []

    for group in rule.split(' or '):

        subunits = sorted(g.strip() for g in group.split(' and ') if g.strip())

        if subunits:
            enzymes.append('_'.join(subunits))

    return enzymes


def metatlas_gem_network(
        gem: str = 'Human-GEM',
        ref: str | None = None,
        include_orphans: bool = False,
) -> Generator[GemInteraction, None, None]:
    """
    Build a binary interaction network from a standard GEM.

    Each metabolic reaction is decomposed into metabolite-enzyme edges:
    reactant → enzyme and enzyme → product.  Reversible reactions produce
    additional edges with ``reverse=True`` where products become sources
    and reactants become targets.

    Gene rules are parsed: ``or`` gives independent isoenzyme edges,
    ``and`` gives per-subunit edges for enzyme complexes.

    Orphan reactions (no gene rule) are skipped by default.  When
    ``include_orphans`` is ``True``, orphan reactions are yielded using
    the reaction ID as a pseudo-enzyme node (``target_type='reaction'``
    or ``source_type='reaction'``).

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.
        include_orphans: If True, yield edges for reactions that have no
            gene rule, using the reaction ID as a pseudo-enzyme node.
            Default: False.

    Yields:
        GemInteraction named tuples — binary edges between metabolites
        and enzymes with compartment annotations.
    """

    data = metatlas_gem_yaml(gem, ref)

    if data is None:
        return

    # Build metabolite ID → compartment lookup
    met_comp = {}

    for met in data.get('metabolites', []):

        if isinstance(met, list):
            met = dict(met)

        met_comp[met.get('id', '')] = met.get('compartment', '')

    reactions = data.get('reactions', [])
    n_orphan = 0

    _log(f'Building network from {len(reactions)} reactions in {gem}.')

    for rxn in reactions:

        if isinstance(rxn, list):
            rxn = dict(rxn)

        enzymes = _parse_gene_rule(rxn.get('gene_reaction_rule', ''))

        if not enzymes:
            n_orphan += 1

            if not include_orphans:
                continue

            rxn_id = rxn.get('id', '')
            mets = rxn.get('metabolites', {})

            if isinstance(mets, list):
                mets = dict(mets)

            lb = float(rxn.get('lower_bound', 0))
            ub = float(rxn.get('upper_bound', 0))
            reversible = lb < 0 < ub
            direction = 1 if lb + ub >= 0 else -1

            reactants = [m for m, c in mets.items() if c * direction < 0]
            products = [m for m, c in mets.items() if c * direction > 0]

            for met_id in reactants:
                yield GemInteraction(
                    source=met_id, target=rxn_id,
                    source_type='metabolite', target_type='reaction',
                    source_compartment=met_comp.get(met_id, ''),
                    target_compartment='',
                    reaction_id=rxn_id, reverse=False,
                )

            for met_id in products:
                yield GemInteraction(
                    source=rxn_id, target=met_id,
                    source_type='reaction', target_type='metabolite',
                    source_compartment='',
                    target_compartment=met_comp.get(met_id, ''),
                    reaction_id=rxn_id, reverse=False,
                )

            if reversible:

                for met_id in products:
                    yield GemInteraction(
                        source=met_id, target=rxn_id,
                        source_type='metabolite', target_type='reaction',
                        source_compartment=met_comp.get(met_id, ''),
                        target_compartment='',
                        reaction_id=rxn_id, reverse=True,
                    )

                for met_id in reactants:
                    yield GemInteraction(
                        source=rxn_id, target=met_id,
                        source_type='reaction', target_type='metabolite',
                        source_compartment='',
                        target_compartment=met_comp.get(met_id, ''),
                        reaction_id=rxn_id, reverse=True,
                    )

            continue

        rxn_id = rxn.get('id', '')

        mets = rxn.get('metabolites', {})

        if isinstance(mets, list):
            mets = dict(mets)

        lb = float(rxn.get('lower_bound', 0))
        ub = float(rxn.get('upper_bound', 0))
        reversible = lb < 0 < ub
        direction = 1 if lb + ub >= 0 else -1

        reactants = [m for m, coef in mets.items() if coef * direction < 0]
        products = [m for m, coef in mets.items() if coef * direction > 0]

        for enzyme in enzymes:

            for met_id in reactants:
                yield GemInteraction(
                    source=met_id,
                    target=enzyme,
                    source_type='metabolite',
                    target_type='protein',
                    source_compartment=met_comp.get(met_id, ''),
                    target_compartment='',
                    reaction_id=rxn_id,
                    reverse=False,
                )

            for met_id in products:
                yield GemInteraction(
                    source=enzyme,
                    target=met_id,
                    source_type='protein',
                    target_type='metabolite',
                    source_compartment='',
                    target_compartment=met_comp.get(met_id, ''),
                    reaction_id=rxn_id,
                    reverse=False,
                )

            if reversible:

                for met_id in products:
                    yield GemInteraction(
                        source=met_id,
                        target=enzyme,
                        source_type='metabolite',
                        target_type='protein',
                        source_compartment=met_comp.get(met_id, ''),
                        target_compartment='',
                        reaction_id=rxn_id,
                        reverse=True,
                    )

                for met_id in reactants:
                    yield GemInteraction(
                        source=enzyme,
                        target=met_id,
                        source_type='protein',
                        target_type='metabolite',
                        source_compartment='',
                        target_compartment=met_comp.get(met_id, ''),
                        reaction_id=rxn_id,
                        reverse=True,
                    )

    if n_orphan:
        _log(f'Skipped {n_orphan} orphan reactions (no gene rule) in {gem}.')


def _strip_compartment(
        met_id: str,
        compartment: str,
) -> str:
    """
    Remove single-letter compartment suffix from a metabolite ID.

    GEM metabolite IDs encode the compartment as a trailing letter
    (e.g. ``MAM01039c`` → ``MAM01039``).  The suffix is stripped only
    when it matches ``compartment``.

    Args:
        met_id: Full metabolite ID including compartment suffix.
        compartment: Compartment code to strip (e.g. ``'c'``).

    Returns:
        Base metabolite ID without compartment suffix, or the original
        ID if the suffix does not match.
    """

    if compartment and met_id.endswith(compartment):
        return met_id[: -len(compartment)]

    return met_id


def _crossing_metabolites(
        edges: list,
) -> set[str]:
    """
    Return base metabolite IDs that cross a compartment boundary.

    Given all GemInteraction edges for a single transport reaction, finds
    metabolites whose base ID (compartment suffix stripped) appears on
    both the reactant side (``source_type == 'metabolite'``) and the
    product side (``target_type == 'metabolite'``) with different
    compartment codes.

    Args:
        edges: All GemInteraction records for a single reaction.

    Returns:
        Set of base metabolite IDs that physically cross a membrane.
    """

    reactant_comps: dict[str, str] = {}
    product_comps: dict[str, str] = {}

    for rec in edges:

        if rec.source_type == 'metabolite':
            base = _strip_compartment(rec.source, rec.source_compartment)
            reactant_comps[base] = rec.source_compartment

        if rec.target_type == 'metabolite':
            base = _strip_compartment(rec.target, rec.target_compartment)
            product_comps[base] = rec.target_compartment

    return {
        base
        for base, comp in reactant_comps.items()
        if base in product_comps and product_comps[base] != comp
    }


def metatlas_gem_transport_ids(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> frozenset[str]:
    """
    Return the set of transport reaction IDs for a GEM.

    Reads the YAML model and collects every reaction whose ``subsystem``
    field equals ``'Transport reactions'``.

    Args:
        gem: GEM name (e.g. ``'Human-GEM'``).
        ref: Git reference; uses the repository's default branch if ``None``.

    Returns:
        Frozenset of reaction ID strings.
    """

    return frozenset(
        rxn.id
        for rxn in metatlas_gem_yaml_reactions(gem=gem, ref=ref)
        if rxn.subsystem == 'Transport reactions'
    )


_ENSG_RE = re.compile(r'^ENSG\d{11}$')


def metatlas_gem_detect_gene_id_type(
        gem: str = 'Human-GEM',
        ref: str | None = None,
) -> str:
    """
    Detect the gene identifier type used in a GEM.

    Human-GEM uses Ensembl gene IDs (``ENSG00000000000``).  Other
    MetAtlas GEMs (Mouse-GEM, Rat-GEM, etc.) use gene symbols.  The type
    is determined by inspecting the first non-empty gene token in the
    gene-reaction rules.

    Args:
        gem: GEM name (e.g. ``'Human-GEM'``, ``'Mouse-GEM'``).
        ref: Git reference; uses the repository's default branch if ``None``.

    Returns:
        ``'ensembl'`` if genes match the ENSG pattern,
        ``'genesymbol'`` otherwise.
    """

    for rxn in metatlas_gem_yaml_reactions(gem=gem, ref=ref):
        rule = rxn.gene_reaction_rule or ''
        token = re.split(r'[\s()]+', rule.strip())[0]

        if token and token.lower() not in ('or', 'and'):
            return 'ensembl' if _ENSG_RE.match(token) else 'genesymbol'

    return 'ensembl'


def metatlas_gem_transport_network(
        gem: str = 'Human-GEM',
        ref: str | None = None,
        include_reverse: bool = True,
        include_orphans: bool = False,
) -> Generator[GemInteraction, None, None]:
    """
    Yield GemInteraction records for transport reactions in a GEM.

    Transport reactions are identified by the subsystem annotation
    ``'Transport reactions'``.  For each such reaction the
    compartment-crossing filter is applied: only metabolites whose base
    ID appears on both the reactant and product sides with *different*
    compartment codes generate edges.  This mechanically excludes
    cofactors (e.g. ATP hydrolysed and regenerated in the same
    compartment) while preserving genuine substrates that physically
    cross a membrane.

    Orphan reactions (no gene rule) are handled as pseudo-enzyme nodes
    using the reaction ID when ``include_orphans`` is ``True``.

    Args:
        gem:
            GEM name (e.g. ``'Human-GEM'``, ``'Mouse-GEM'``).
        ref:
            Git reference (branch, tag, or commit).  ``None`` uses the
            repository's default branch.
        include_reverse:
            If ``True``, yield reversed edges for reversible reactions
            (``reverse=True``).  Default: ``True``.
        include_orphans:
            If ``True``, yield edges for transport reactions that have no
            gene rule, using the reaction ID as a pseudo-enzyme node with
            ``target_type='reaction'`` or ``source_type='reaction'``.
            Default: ``False``.

    Yields:
        :class:`GemInteraction` records.  For metabolite → enzyme edges,
        ``source_type='metabolite'`` and ``target_type='protein'`` (or
        ``'reaction'`` for orphan pseudo-enzyme nodes).  Roles are
        swapped for enzyme → metabolite edges.
    """

    data = metatlas_gem_yaml(gem, ref)

    if data is None:
        return

    # Build metabolite ID → compartment lookup.
    met_comp: dict[str, str] = {}

    for met in data.get('metabolites', []):

        if isinstance(met, list):
            met = dict(met)

        met_comp[met.get('id', '')] = met.get('compartment', '')

    reactions = data.get('reactions', [])
    n_transport = 0

    _log(f'Building transport network from {len(reactions)} reactions in {gem}.')

    for rxn in reactions:

        if isinstance(rxn, list):
            rxn = dict(rxn)

        subsystem = rxn.get('subsystem', '')

        if isinstance(subsystem, list):
            subsystem = '; '.join(str(s) for s in subsystem)

        if subsystem != 'Transport reactions':
            continue

        rxn_id = rxn.get('id', '')
        enzymes = _parse_gene_rule(rxn.get('gene_reaction_rule', ''))

        if not enzymes and not include_orphans:
            continue

        mets = rxn.get('metabolites', {})

        if isinstance(mets, list):
            mets = dict(mets)

        lb = float(rxn.get('lower_bound', 0))
        ub = float(rxn.get('upper_bound', 0))
        reversible = lb < 0 < ub
        direction = 1 if lb + ub >= 0 else -1

        reactants = [m for m, coef in mets.items() if coef * direction < 0]
        products = [m for m, coef in mets.items() if coef * direction > 0]

        node_type = 'protein' if enzymes else 'reaction'
        nodes = enzymes if enzymes else [rxn_id]

        # Collect candidate edges; crossing-metabolite filter applied after.
        raw: list[GemInteraction] = []

        for node in nodes:

            for met_id in reactants:
                raw.append(GemInteraction(
                    source=met_id, target=node,
                    source_type='metabolite', target_type=node_type,
                    source_compartment=met_comp.get(met_id, ''),
                    target_compartment='',
                    reaction_id=rxn_id, reverse=False,
                ))

            for met_id in products:
                raw.append(GemInteraction(
                    source=node, target=met_id,
                    source_type=node_type, target_type='metabolite',
                    source_compartment='',
                    target_compartment=met_comp.get(met_id, ''),
                    reaction_id=rxn_id, reverse=False,
                ))

            if reversible and include_reverse:

                for met_id in products:
                    raw.append(GemInteraction(
                        source=met_id, target=node,
                        source_type='metabolite', target_type=node_type,
                        source_compartment=met_comp.get(met_id, ''),
                        target_compartment='',
                        reaction_id=rxn_id, reverse=True,
                    ))

                for met_id in reactants:
                    raw.append(GemInteraction(
                        source=node, target=met_id,
                        source_type=node_type, target_type='metabolite',
                        source_compartment='',
                        target_compartment=met_comp.get(met_id, ''),
                        reaction_id=rxn_id, reverse=True,
                    ))

        # Apply compartment-crossing filter.
        crossing = _crossing_metabolites(raw)

        if not crossing:
            continue

        n_transport += 1

        for rec in raw:
            met = rec.source if rec.source_type == 'metabolite' else rec.target
            comp = (
                rec.source_compartment
                if rec.source_type == 'metabolite'
                else rec.target_compartment
            )
            base = _strip_compartment(met, comp)

            if base in crossing:
                yield rec

    _log(f'{gem}: {n_transport} transport reactions yielded.')
