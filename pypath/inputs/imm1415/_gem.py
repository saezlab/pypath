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
iMM1415 parser: metabolites, reactions, genes, and interaction network.

iMM1415 (Sigurdsson et al. 2010) is a genome-scale metabolic model of Mus
musculus distributed via BiGG Models.  It uses the same JSON format as
Recon3D; the BiGG-generic parsing helpers are imported directly from
``pypath.inputs.recon3d._gem`` to avoid duplication.

Gene identifiers in iMM1415 are NCBI Entrez Gene IDs for mouse (integer
strings), analogous to Recon3D's human Entrez IDs.

The :func:`imm1415_network` and :func:`imm1415_transporter_network`
generators yield :class:`GemInteraction` records compatible with the
``omnipath-metabo`` GEM processor.

References:
    Sigurdsson MI, Jamshidi N, Steingrimsson E, et al.  A detailed
    genome-wide reconstruction of mouse metabolism based on human Recon 1.
    BMC Syst Biol. 2010;4:140. doi:10.1186/1752-0509-4-140
"""

from __future__ import annotations

__all__ = [
    'imm1415_metabolites',
    'imm1415_reactions',
    'imm1415_genes',
    'imm1415_network',
    'imm1415_transporter_network',
]

from collections import defaultdict
from collections.abc import Generator

from pypath.inputs.metatlas._records import GemInteraction

# Reuse the BiGG-generic parsing helpers from the Recon3D module.
from pypath.inputs.recon3d._gem import (
    _annotation_list,
    _parse_gene_rule,
    _strip_compartment,
)

from ._common import _log
from ._raw import imm1415_raw


def imm1415_metabolites() -> list[dict]:
    """
    Return iMM1415 metabolites with cross-reference annotations.

    Each returned dict has keys: ``id``, ``base_id``, ``compartment``,
    ``name``, ``formula``, ``charge``, ``hmdb``, ``chebi``, ``kegg``,
    ``bigg``, ``metanetx``.

    Returns:
        List of metabolite dicts.
    """

    data = imm1415_raw()
    metabolites = data.get('metabolites', [])
    _log(f'iMM1415: parsing {len(metabolites)} metabolites from BiGG JSON.')
    result = []

    for met in metabolites:

        full_id = met.get('id', '')
        base_id, compartment = _strip_compartment(full_id)
        annotation = met.get('annotation', {})

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
            'hmdb': _annotation_list(annotation, 'hmdb'),
            'chebi': _annotation_list(annotation, 'chebi'),
            'kegg': _annotation_list(annotation, 'kegg.compound'),
            'bigg': _annotation_list(annotation, 'bigg.metabolite'),
            'metanetx': _annotation_list(annotation, 'metanetx.chemical'),
        })

    return result


def imm1415_reactions() -> list[dict]:
    """
    Return iMM1415 reactions with stoichiometry, bounds, and gene rules.

    Each returned dict has keys: ``id``, ``name``, ``metabolites``,
    ``lower_bound``, ``upper_bound``, ``reversible``,
    ``gene_reaction_rule``, ``subsystem``, ``annotation``.

    Returns:
        List of reaction dicts.
    """

    data = imm1415_raw()
    reactions = data.get('reactions', [])
    _log(f'iMM1415: parsing {len(reactions)} reactions from BiGG JSON.')
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


def imm1415_genes() -> list[dict]:
    """
    Return iMM1415 genes with cross-reference annotations.

    Gene IDs are NCBI Entrez Gene IDs for mouse (integer strings).
    Each returned dict has keys: ``id``, ``name``, ``annotation``.

    Returns:
        List of gene dicts.
    """

    data = imm1415_raw()
    genes = data.get('genes', [])
    _log(f'iMM1415: parsed {len(genes)} genes from BiGG JSON.')

    return [
        {
            'id': g.get('id', ''),
            'name': g.get('name', ''),
            'annotation': g.get('annotation', {}),
        }
        for g in genes
    ]


def imm1415_network(
        metab_max_degree: int = 400,
        include_reverse: bool = True,
) -> Generator[GemInteraction, None, None]:
    """
    Yield iMM1415 metabolite-enzyme interactions as :class:`GemInteraction`
    records.

    Same edge-building convention as :func:`pypath.inputs.recon3d.recon3d_network`:
    reactants → enzyme and enzyme → products, with reversed pairs for
    reversible reactions.

    Args:
        metab_max_degree:
            Cofactor filter threshold.  Default: 400.
        include_reverse:
            If ``True``, yield reversed edges for reversible reactions.

    Yields:
        :class:`~pypath.inputs.metatlas._records.GemInteraction` named tuples.
        Gene IDs are mouse NCBI Entrez IDs (integer strings).
    """

    from collections import Counter

    reactions = imm1415_reactions()
    n_orphan = 0
    raw: list[tuple] = []

    for rxn in reactions:

        enzymes = _parse_gene_rule(rxn['gene_reaction_rule'])

        if not enzymes:
            n_orphan += 1
            continue

        mets: dict = rxn['metabolites']
        lb = rxn['lower_bound']
        ub = rxn['upper_bound']
        rxn_id = rxn['id']
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
        f'iMM1415: {len(reactions)} reactions, {n_orphan} orphan (skipped), '
        f'{len(raw)} raw edges before cofactor filtering.'
    )

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

    _log(f'iMM1415: {n_filtered} edges removed by cofactor filter (metab_max_degree={metab_max_degree}).')


def imm1415_transporter_network(
        include_reverse: bool = True,
        include_orphans: bool = False,
) -> Generator[GemInteraction, None, None]:
    """
    Yield GemInteraction records for iMM1415 transport reactions.

    Transport reactions are identified by detecting metabolites whose
    BiGG base ID appears on both the reactant and product sides with
    different compartment codes.

    Args:
        include_reverse:
            If ``True``, yield reversed edges for reversible transport
            reactions.  Default: ``True``.
        include_orphans:
            If ``True``, yield edges for orphan transport reactions using
            the reaction ID as a pseudo-enzyme node.  Default: ``False``.

    Yields:
        :class:`GemInteraction` records.  Gene IDs are mouse NCBI Entrez
        IDs (integer strings) with ``_ATN`` isoform suffixes removed.
    """

    reactions = imm1415_reactions()
    n_transport = 0

    _log(f'iMM1415: building transporter network from {len(reactions)} reactions.')

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

    _log(f'iMM1415: {n_transport} transport reactions processed.')
