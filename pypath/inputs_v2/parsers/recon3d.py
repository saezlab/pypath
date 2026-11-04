"""Recon3D's raw parser, over the generic BiGG JSON parsing in ``_bigg``.

Recon3D's JSON carries nothing the generic parser doesn't already handle —
every name below is defined in :mod:`pypath.inputs_v2.parsers._bigg` and
re-exported here unchanged, so existing imports of this module keep working.
"""

from pypath.inputs_v2.parsers._bigg import (
    _annotation_list,
    _build_metabolite_index,
    _merge_unique,
    _normalize_annotation_value,
    _parse_catalysis,
    _parse_enzyme_complexes,
    _parse_gene_rule,
    _parse_genes,
    _parse_metabolites,
    _parse_reactions,
    _raw,
    _reaction_enzyme_entrez_ids,
    _reaction_member_fields,
    _reaction_metabolite_ids,
    _serialize_member_value,
    _strip_compartment,
    _strip_isoform,
)

__all__ = [
    '_annotation_list',
    '_normalize_annotation_value',
    '_strip_compartment',
    '_strip_isoform',
    '_merge_unique',
    '_build_metabolite_index',
    '_reaction_metabolite_ids',
    '_serialize_member_value',
    '_reaction_member_fields',
    '_parse_gene_rule',
    '_parse_reactions',
    '_reaction_enzyme_entrez_ids',
    '_parse_metabolites',
    '_parse_genes',
    '_parse_catalysis',
    '_parse_enzyme_complexes',
    '_raw',
]
