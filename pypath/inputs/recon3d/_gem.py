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
Recon3D parser: metabolites, reactions, genes, and interaction network.

Primary data source is the BiGG JSON model.  Metabolite annotations
(HMDB, ChEBI, KEGG) are extracted directly from the JSON ``annotation``
field; an optional VMH MATLAB supplement can add HMDB IDs not present in
the BiGG file.

Gene identifiers in Recon3D are NCBI Entrez Gene IDs (integers stored as
strings in the ``gene_reaction_rule`` expressions).

The :func:`recon3d_network` generator yields :class:`GemInteraction`
records (same namedtuple used by ``pypath.inputs.metatlas``), making
Recon3D compatible with the ``omnipath-metabo`` GEM processor.
"""

from __future__ import annotations

__all__ = [
    'recon3d_metabolites',
    'recon3d_reactions',
    'recon3d_genes',
    'recon3d_network',
    'recon3d_transporter_network',
]

import re
from collections import defaultdict
from collections.abc import Generator

from pypath.inputs.metatlas._records import GemInteraction

from ._common import _log
from ._raw import gem_matlab_extract, recon3d_raw, recon3d_raw_vmh


# ── helpers ──────────────────────────────────────────────────────────────────


def _strip_compartment(met_id: str) -> tuple[str, str]:
    """
    Split a BiGG metabolite ID into (base_id, compartment).

    BiGG encodes compartment as a single-letter suffix after the last
    underscore: ``'10fthf_c'`` → ``('10fthf', 'c')``.  If the last segment
    is not a single letter, the full ID is returned with an empty compartment.

    Args:
        met_id: BiGG metabolite ID with compartment suffix.

    Returns:
        Tuple of (base_id, compartment_code).
    """

    parts = met_id.rsplit('_', 1)

    if len(parts) == 2 and len(parts[1]) == 1 and parts[1].isalpha():
        return parts[0], parts[1]

    return met_id, ''


def _parse_gene_rule(
        rule: str,
        strip_isoforms: bool = True,
) -> list[str]:
    """
    Parse a Recon3D ``gene_reaction_rule`` into enzyme identifiers.

    ``or`` → isoenzymes (each yielded separately), ``and`` → complex
    subunits (joined with ``_``).  Parentheses are stripped.

    Recon3D gene IDs carry ``_ATN`` isoform suffixes (e.g. ``1234_AT1``).
    When ``strip_isoforms`` is ``True`` (default) these are removed before
    building complex subunit strings, yielding clean Entrez IDs.

    Args:
        rule: Gene-reaction rule string (e.g. ``'1234_AT1 or (5678_AT1 and 9012_AT2)'``).
        strip_isoforms: Remove ``_ATN`` isoform suffixes.  Default: ``True``.

    Returns:
        List of enzyme identifier strings.  Returns empty list for orphan
        reactions (empty or whitespace-only rule).
    """

    if not rule or not rule.strip():
        return []

    rule = re.sub(r'[()]', '', rule).strip()
    enzymes = []

    for or_part in re.split(r'\bor\b', rule, flags=re.IGNORECASE):
        genes = [g.strip() for g in re.split(r'\band\b', or_part, flags=re.IGNORECASE)]
        genes = [g for g in genes if g]

        if not genes:
            continue

        if strip_isoforms:
            stripped = []

            for g in genes:
                m = re.search(r'_AT\d+$', g)
                stripped.append(g[: m.start()] if m else g)

            genes = [g for g in stripped if g]

        if not genes:
            continue

        enzymes.append('_'.join(sorted(genes)) if len(genes) > 1 else genes[0])

    return enzymes


def _annotation_list(annotation: dict | list | None, key: str) -> list[str]:
    """
    Extract a flat list of annotation values from a BiGG annotation dict.

    BiGG annotation values can be strings or lists of strings.

    Args:
        annotation: The ``annotation`` field from a BiGG metabolite or gene.
        key: Annotation key (e.g. ``'hmdb'``, ``'chebi'``, ``'kegg.compound'``).

    Returns:
        List of strings.  Empty list if key is absent.
    """

    if not isinstance(annotation, dict):
        return []

    val = annotation.get(key, [])

    if isinstance(val, str):
        return [val]

    return list(val)


# ── public functions ──────────────────────────────────────────────────────────


def recon3d_metabolites(extra_hmdb: bool = True) -> list[dict]:
    """
    Return Recon3D metabolites with cross-reference annotations.

    Parses the BiGG JSON model.  When ``extra_hmdb`` is ``True``, the VMH
    MATLAB file is downloaded and used to supplement HMDB IDs for metabolites
    that have none or fewer IDs in the BiGG JSON (mirrors OmnipathR's
    ``recon3d_metabolites(extra_hmdb = TRUE)``).

    Each returned dict has keys:

    - ``id`` — full BiGG ID including compartment (e.g. ``'10fthf_c'``)
    - ``base_id`` — BiGG ID without compartment (e.g. ``'10fthf'``)
    - ``compartment`` — single-letter code (e.g. ``'c'``)
    - ``name`` — metabolite name
    - ``formula`` — molecular formula
    - ``charge`` — formal charge (int or None)
    - ``hmdb`` — list of HMDB IDs
    - ``chebi`` — list of ChEBI IDs
    - ``kegg`` — list of KEGG compound IDs
    - ``bigg`` — list of BiGG metabolite IDs (cross-references)

    Args:
        extra_hmdb:
            If ``True``, enrich HMDB IDs from the VMH MATLAB file.
            Requires ``scipy``.

    Returns:
        List of metabolite dicts.
    """

    data = recon3d_raw()
    metabolites = data.get('metabolites', [])
    _log(f'Recon3D: parsing {len(metabolites)} metabolites from BiGG JSON.')

    # Build VMH HMDB supplement: base_id → list[hmdb_id]
    vmh_hmdb: dict[str, list[str]] = {}

    if extra_hmdb:

        try:
            mat = recon3d_raw_vmh()
            extracted = gem_matlab_extract(mat, 'mets', 'metHMDBID')
            mets_raw = extracted.get('mets', [])
            hmdb_raw = extracted.get('metHMDBID', [])

            for met_full, hmdb_val in zip(mets_raw, hmdb_raw):

                if isinstance(met_full, str):
                    base, _ = _strip_compartment(met_full)
                else:
                    continue

                if isinstance(hmdb_val, str) and hmdb_val:
                    vmh_hmdb.setdefault(base, []).append(hmdb_val)

                elif isinstance(hmdb_val, list):
                    for v in hmdb_val:
                        if v:
                            vmh_hmdb.setdefault(base, []).append(str(v))

            _log(f'Recon3D: VMH MATLAB supplement loaded for {len(vmh_hmdb)} base metabolites.')

        except Exception as exc:
            _log(f'Recon3D: VMH MATLAB supplement failed ({exc}); proceeding without it.')

    result = []

    for met in metabolites:

        full_id = met.get('id', '')
        base_id, compartment = _strip_compartment(full_id)
        annotation = met.get('annotation', {})

        hmdb = _annotation_list(annotation, 'hmdb')
        chebi = _annotation_list(annotation, 'chebi')
        kegg = _annotation_list(annotation, 'kegg.compound')
        bigg = _annotation_list(annotation, 'bigg.metabolite')
        metanetx = _annotation_list(annotation, 'metanetx.chemical')

        if extra_hmdb and base_id in vmh_hmdb:
            hmdb = list(dict.fromkeys(hmdb + vmh_hmdb[base_id]))

        charge = met.get('charge')

        if charge is not None:
            try:
                charge = int(charge)
            except (TypeError, ValueError):
                charge = None

        result.append({
            'id': full_id,
            'base_id': base_id,
            'compartment': compartment,
            'name': met.get('name', ''),
            'formula': met.get('formula', ''),
            'charge': charge,
            'hmdb': hmdb,
            'chebi': chebi,
            'kegg': kegg,
            'bigg': bigg,
            'metanetx': metanetx,
        })

    return result


def recon3d_reactions() -> list[dict]:
    """
    Return Recon3D reactions with stoichiometry, bounds, and gene rules.

    Parses the BiGG JSON model.  Each returned dict has keys:

    - ``id`` — reaction BiGG ID
    - ``name`` — reaction name
    - ``metabolites`` — dict mapping metabolite full ID → stoich coefficient
    - ``lower_bound`` — flux lower bound (float)
    - ``upper_bound`` — flux upper bound (float)
    - ``reversible`` — ``True`` if ``lower_bound < 0 < upper_bound``
    - ``gene_reaction_rule`` — raw Boolean gene rule string
    - ``subsystem`` — subsystem/pathway name (str or ``''``)
    - ``annotation`` — dict of cross-reference annotations

    Returns:
        List of reaction dicts.
    """

    data = recon3d_raw()
    reactions = data.get('reactions', [])
    _log(f'Recon3D: parsing {len(reactions)} reactions from BiGG JSON.')

    result = []

    for rxn in reactions:

        lb = float(rxn.get('lower_bound', 0))
        ub = float(rxn.get('upper_bound', 0))

        result.append({
            'id': rxn.get('id', ''),
            'name': rxn.get('name', ''),
            'metabolites': rxn.get('metabolites', {}),
            'lower_bound': lb,
            'upper_bound': ub,
            'reversible': lb < 0 < ub,
            'gene_reaction_rule': rxn.get('gene_reaction_rule', ''),
            'subsystem': rxn.get('subsystem', ''),
            'annotation': rxn.get('annotation', {}),
        })

    return result


def recon3d_genes() -> list[dict]:
    """
    Return Recon3D genes with cross-reference annotations.

    Parses the BiGG JSON model.  Gene IDs are NCBI Entrez Gene IDs stored
    as strings.  Each returned dict has keys:

    - ``id`` — Entrez Gene ID string (e.g. ``'1'``)
    - ``name`` — gene name / symbol
    - ``annotation`` — dict of cross-reference annotations (may include
      ``'ncbigene'``, ``'uniprot'``, ``'ensembl.gene'``, etc.)

    Returns:
        List of gene dicts.
    """

    data = recon3d_raw()
    genes = data.get('genes', [])
    _log(f'Recon3D: parsed {len(genes)} genes from BiGG JSON.')

    return [
        {
            'id': g.get('id', ''),
            'name': g.get('name', ''),
            'annotation': g.get('annotation', {}),
        }
        for g in genes
    ]


def recon3d_network(
        metab_max_degree: int = 400,
        include_reverse: bool = True,
) -> Generator[GemInteraction, None, None]:
    """
    Yield Recon3D metabolite-enzyme interactions as :class:`GemInteraction`
    records.

    Each reaction is decomposed into binary directed edges using the same
    convention as :func:`pypath.inputs.metatlas.metatlas_gem_network`:

    - *metabolite → enzyme*: reactant consumed by the enzyme.
    - *enzyme → metabolite*: product released by the enzyme.

    Reversible reactions (``lower_bound < 0 < upper_bound``) additionally
    produce a reversed pair with ``reverse=True``.  Orphan reactions (no
    gene rule) are skipped.

    High-degree metabolites (likely cofactors) are filtered in a two-pass
    approach: all edges are collected first, metabolite degrees counted, then
    only edges involving metabolites with degree ≤ *metab_max_degree* are
    yielded.

    Metabolite IDs are BiGG base IDs (compartment suffix stripped); enzyme
    IDs are Entrez Gene ID strings, or ``'ENTREZ1_ENTREZ2'`` strings for
    AND-rule enzyme complexes.  The ``reaction_id`` field carries the BiGG
    reaction ID.

    Args:
        metab_max_degree:
            Cofactor filter threshold.  Default: 400.
        include_reverse:
            If ``True``, yield reversed edges for reversible reactions.

    Yields:
        :class:`~pypath.inputs.metatlas._records.GemInteraction` named tuples.
    """

    from collections import Counter

    reactions = recon3d_reactions()
    n_orphan = 0

    # Two-pass: collect raw edges, then filter by metabolite degree.
    raw: list[tuple] = []  # (source, target, source_type, target_type, src_comp, tgt_comp, rxn_id, reverse)

    for rxn in reactions:

        enzymes = _parse_gene_rule(rxn['gene_reaction_rule'])

        if not enzymes:
            n_orphan += 1
            continue

        mets: dict = rxn['metabolites']
        lb = rxn['lower_bound']
        ub = rxn['upper_bound']
        rxn_id = rxn['id']

        # Direction sign: +1 if forward-preferred, -1 if reverse-preferred.
        direction = 1 if lb + ub >= 0 else -1

        reactants = [m for m, coef in mets.items() if coef * direction < 0]
        products = [m for m, coef in mets.items() if coef * direction > 0]
        reversible = rxn['reversible']

        for enzyme in enzymes:

            for met_full in reactants:
                base, comp = _strip_compartment(met_full)
                raw.append((base, enzyme, 'metabolite', 'protein', comp, '', rxn_id, False))

                if reversible and include_reverse:
                    raw.append((enzyme, base, 'protein', 'metabolite', '', comp, rxn_id, True))

            for met_full in products:
                base, comp = _strip_compartment(met_full)
                raw.append((enzyme, base, 'protein', 'metabolite', '', comp, rxn_id, False))

                if reversible and include_reverse:
                    raw.append((base, enzyme, 'metabolite', 'protein', comp, '', rxn_id, True))

    _log(
        f'Recon3D: {len(reactions)} reactions, {n_orphan} orphan (skipped), '
        f'{len(raw)} raw edges before cofactor filtering.'
    )

    # Count metabolite degree across all edges.
    metab_degree: Counter = Counter()

    for src, tgt, src_type, tgt_type, *_ in raw:

        if src_type == 'metabolite':
            metab_degree[src] += 1

        if tgt_type == 'metabolite':
            metab_degree[tgt] += 1

    n_filtered = 0

    for src, tgt, src_type, tgt_type, src_comp, tgt_comp, rxn_id, reverse in raw:

        met = src if src_type == 'metabolite' else tgt

        if metab_degree[met] > metab_max_degree:
            n_filtered += 1
            continue

        yield GemInteraction(
            source=src,
            target=tgt,
            source_type=src_type,
            target_type=tgt_type,
            source_compartment=src_comp,
            target_compartment=tgt_comp,
            reaction_id=rxn_id,
            reverse=reverse,
        )

    _log(f'Recon3D: {n_filtered} edges removed by cofactor filter (metab_max_degree={metab_max_degree}).')


def recon3d_transporter_network(
        include_reverse: bool = True,
        include_orphans: bool = False,
) -> Generator[GemInteraction, None, None]:
    """
    Yield GemInteraction records for Recon3D transport reactions.

    Transport reactions are identified by detecting metabolites whose
    BiGG base ID appears on both the reactant and product sides of a
    reaction with *different* compartment codes — the molecular signature
    of membrane transport.

    Two directed edges are generated per transported metabolite per
    enzyme:

    - ``met[in_comp] → enzyme``
    - ``enzyme → met[out_comp]``

    Reversible reactions produce an additional reversed pair with
    ``reverse=True`` where ``in_comp`` and ``out_comp`` are swapped.

    Orphan reactions (no gene rule) generate edges using the reaction ID
    as a pseudo-enzyme node when ``include_orphans`` is ``True``.

    Args:
        include_reverse:
            If ``True``, yield reversed edges for reversible transport
            reactions (``reverse=True``).  Default: ``True``.
        include_orphans:
            If ``True``, yield edges for orphan transport reactions
            using the reaction ID as a pseudo-enzyme node with
            ``target_type='reaction'`` or ``source_type='reaction'``.
            Default: ``False``.

    Yields:
        :class:`GemInteraction` records.  Metabolite IDs are BiGG base
        IDs (compartment suffix stripped); enzyme IDs are Entrez Gene ID
        strings with ``_ATN`` isoform suffixes removed.
    """

    reactions = recon3d_reactions()
    n_transport = 0

    _log(f'Recon3D: building transporter network from {len(reactions)} reactions.')

    for rxn in reactions:

        enzymes = _parse_gene_rule(rxn['gene_reaction_rule'])

        if not enzymes and not include_orphans:
            continue

        mets: dict = rxn['metabolites']
        lb = rxn['lower_bound']
        ub = rxn['upper_bound']
        rxn_id = rxn['id']
        reversible = rxn['reversible']
        direction = 1 if lb + ub >= 0 else -1

        reactant_comps: defaultdict[str, list[str]] = defaultdict(list)
        product_comps: defaultdict[str, list[str]] = defaultdict(list)

        for met_full, coef in mets.items():
            base, comp = _strip_compartment(met_full)

            if coef * direction < 0:
                reactant_comps[base].append(comp)
            elif coef * direction > 0:
                product_comps[base].append(comp)

        transported = [
            (base, in_comp, out_comp)
            for base in set(reactant_comps) & set(product_comps)
            for in_comp in reactant_comps[base]
            for out_comp in product_comps[base]
            if in_comp != out_comp
        ]

        if not transported:
            continue

        n_transport += 1
        nodes = enzymes if enzymes else [rxn_id]
        node_type = 'protein' if enzymes else 'reaction'

        for base_id, in_comp, out_comp in transported:
            for node in nodes:

                yield GemInteraction(
                    source=base_id, target=node,
                    source_type='metabolite', target_type=node_type,
                    source_compartment=in_comp, target_compartment='',
                    reaction_id=rxn_id, reverse=False,
                )
                yield GemInteraction(
                    source=node, target=base_id,
                    source_type=node_type, target_type='metabolite',
                    source_compartment='', target_compartment=out_comp,
                    reaction_id=rxn_id, reverse=False,
                )

                if reversible and include_reverse:

                    yield GemInteraction(
                        source=base_id, target=node,
                        source_type='metabolite', target_type=node_type,
                        source_compartment=out_comp, target_compartment='',
                        reaction_id=rxn_id, reverse=True,
                    )
                    yield GemInteraction(
                        source=node, target=base_id,
                        source_type=node_type, target_type='metabolite',
                        source_compartment='', target_compartment=in_comp,
                        reaction_id=rxn_id, reverse=True,
                    )

    _log(f'Recon3D: {n_transport} transport reactions processed.')
