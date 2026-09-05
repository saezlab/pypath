"""Recon3D metabolic activities, chemical participants and model genes.

Input/output edges retain stoichiometry; chemical formula and charge use typed
attributes. Compartments, numeric bounds and original Boolean gene rules remain source
payload fields. GPR alternatives link activities to source-scoped logical AND
groups of genes; these groups do not assert physical assemblies.
"""

from __future__ import annotations

import hashlib
import json

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.recon3d import _raw
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Membership
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
)


def _gene_rule_group(genes: list[str]) -> Entity:
    """Source-scoped logical AND clause, without a physical assembly claim."""
    genes = sorted(set(genes))
    digest = hashlib.sha256(json.dumps(genes, separators=(',', ':')).encode()).hexdigest()
    return Entity(
        type=model.NamedThing,
        identifiers=[Identifier(type='recon3d.gpr_and', value=digest)],
        membership=[
            Membership(
                member=Entity(
                    type=model.Gene,
                    identifiers=[Identifier(type=Namespace.ENTREZ, value=gene)],
                    annotations=[Annotation(term=slots.in_taxon, value='NCBITaxon:9606')],
                ),
                predicate=slots.has_member,
            )
            for gene in genes
        ],
    )


class _ReactionMembership(MembershipBuilder):
    """Attach alternative GPR clauses, retaining required genes within each."""

    def build(self, row, cache):
        members = super().build(row, cache)
        for genes in row.get('gene_rule_clauses', []):
            if not genes:
                continue
            group = _gene_rule_group(genes)
            member = group.membership[0].member if len(genes) == 1 else group
            members.append(Membership(member=member, predicate=slots.associated_with))
        return members


config = ResourceConfig(
    id=ResourceCv.RECON3D,
    name='Recon3D',
    url='https://www.vmh.life/',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='29457794',
    primary_category='metabolic_network',
    description=(
        'Recon3D is a comprehensive human genome-scale metabolic reconstruction '
        'that accounts for three-dimensional molecular structures.'
    ),
)

f = FieldConfig(
    delimiter=';',
    map={
        'entity_type': lambda value: model.MacromolecularComplex
        if value == 'complex'
        else model.Protein,
        'stoich_id': lambda value: value.split(':')[0] if value else None,
        'stoich_comp': lambda value: value.split(':')[1] if value else None,
        'stoich_val': lambda value: value.split(':')[2] if value else None,
        'split_member_values': lambda value: value.split(';;')
        if value
        else None,
        'complex_name': lambda value: 'complex:' + '-'.join(value.split('||'))
        if value
        else '',
    },
)


download = Download(
    url='http://bigg.ucsd.edu/static/models/Recon3D.json',
    filename='Recon3D.json',
    subfolder='recon3d',
    large=True,
    ext='json',
    default_mode='r',
)

HUMAN_TAXON_ID = '9606'


# ── metabolites ──────────────────────────────────────────────────────────────

metabolites_schema = EntityBuilder(
    entity_type=model.ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.BIGG_METABOLITE, value=f('bigg_metabolite_id')),
        CV(term=Namespace.HMDB, value=f('hmdb')),
        CV(term=Namespace.CHEBI, value=f('chebi')),
        CV(term=Namespace.KEGG, value=f('kegg_compound')),
        CV(term=Namespace.METANETX, value=f('metanetx')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term='chemrof:charge', value=f('charge')),
        CV(term=slots.has_chemical_formula, value=f('formula')),
    ),
)


# ── reactions ────────────────────────────────────────────────────────────────

reactions_schema = EntityBuilder(
    entity_type=model.MolecularActivity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.BIGG_REACTION, value=f('bigg_reaction_id')),
        CV(term=Namespace.METANETX_REACTION, value=f('metanetx_reaction')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_topic,
            value=f(
                'ec', transform=lambda v: 'EC:' + str(v).removeprefix('EC:')
            ),
        )
    ),
    membership=_ReactionMembership(
        MembersFromList(
            entity_type=model.ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.BIGG_METABOLITE,
                    value=f(
                        'reactants',
                        delimiter='||',
                        map='stoich_id',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'reactant_hmdb',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'reactant_chebi',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.KEGG,
                    value=f(
                        'reactant_kegg_compound',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.METANETX,
                    value=f(
                        'reactant_metanetx',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f(
                        'reactant_name', delimiter='||', preserve_indices=True
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f(
                        'reactants',
                        delimiter='||',
                        map='stoich_val',
                        preserve_indices=True,
                    ),
                )
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term='chemrof:charge',
                    value=f(
                        'reactant_charge', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=slots.has_chemical_formula,
                    value=f(
                        'reactant_formula',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
            ),
            predicate=slots.has_input,
        ),
        MembersFromList(
            entity_type=model.ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.BIGG_METABOLITE,
                    value=f(
                        'products',
                        delimiter='||',
                        map='stoich_id',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'product_hmdb',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'product_chebi',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.KEGG,
                    value=f(
                        'product_kegg_compound',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.METANETX,
                    value=f(
                        'product_metanetx',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f(
                        'product_name', delimiter='||', preserve_indices=True
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f(
                        'products',
                        delimiter='||',
                        map='stoich_val',
                        preserve_indices=True,
                    ),
                )
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term='chemrof:charge',
                    value=f(
                        'product_charge', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=slots.has_chemical_formula,
                    value=f(
                        'product_formula', delimiter='||', preserve_indices=True
                    ),
                ),
            ),
            predicate=slots.has_output,
        ),
    ),
)


# ── transport_reactions ──────────────────────────────────────────────────────

transport_reactions_schema = EntityBuilder(
    entity_type=model.MolecularActivity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.BIGG_REACTION, value=f('bigg_reaction_id')),
        CV(term=Namespace.METANETX_REACTION, value=f('metanetx_reaction')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_topic,
            value=f(
                'ec', transform=lambda v: 'EC:' + str(v).removeprefix('EC:')
            ),
        )
    ),
    membership=_ReactionMembership(
        MembersFromList(
            entity_type=model.ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.BIGG_METABOLITE,
                    value=f(
                        'reactants',
                        delimiter='||',
                        map='stoich_id',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'reactant_hmdb',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'reactant_chebi',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.KEGG,
                    value=f(
                        'reactant_kegg_compound',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.METANETX,
                    value=f(
                        'reactant_metanetx',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f(
                        'reactant_name', delimiter='||', preserve_indices=True
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f(
                        'reactants',
                        delimiter='||',
                        map='stoich_val',
                        preserve_indices=True,
                    ),
                )
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term='chemrof:charge',
                    value=f(
                        'reactant_charge', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=slots.has_chemical_formula,
                    value=f(
                        'reactant_formula',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
            ),
            predicate=slots.has_input,
        ),
        MembersFromList(
            entity_type=model.ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.BIGG_METABOLITE,
                    value=f(
                        'products',
                        delimiter='||',
                        map='stoich_id',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'product_hmdb',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'product_chebi',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.KEGG,
                    value=f(
                        'product_kegg_compound',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.METANETX,
                    value=f(
                        'product_metanetx',
                        delimiter='||',
                        map='split_member_values',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f(
                        'product_name', delimiter='||', preserve_indices=True
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f(
                        'products',
                        delimiter='||',
                        map='stoich_val',
                        preserve_indices=True,
                    ),
                )
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term='chemrof:charge',
                    value=f(
                        'product_charge', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=slots.has_chemical_formula,
                    value=f(
                        'product_formula', delimiter='||', preserve_indices=True
                    ),
                ),
            ),
            predicate=slots.has_output,
        ),
    ),
)


# ── catalysis ────────────────────────────────────────────────────────────────


# ── enzyme_complexes ─────────────────────────────────────────────────────────

def enzyme_complexes_schema(row: dict) -> Entity | None:
    """Retain a source model AND group with a deterministic internal ID."""
    genes = [gene for gene in str(row.get('complex_subunits') or '').split('||') if gene]
    return _gene_rule_group(genes) if genes else None


# ── genes ────────────────────────────────────────────────────────────────────

genes_schema = EntityBuilder(
    entity_type=model.Gene,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.ENTREZ, value=f('entrez_id')),
        CV(term=Namespace.GENESYMBOL, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value='NCBITaxon:9606')
    ),
)


# ── resource ─────────────────────────────────────────────────────────────────

resource = Resource(
    config,
    metabolites=Dataset(
        download=download,
        mapper=metabolites_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='metabolites', **kwargs
        ),
    ),
    reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='reactions', **kwargs
        ),
    ),
    transport_reactions=Dataset(
        download=download,
        mapper=transport_reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='transport_reactions', **kwargs
        ),
    ),
    metabolic_reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='metabolic_reactions', **kwargs
        ),
    ),
    enzyme_complexes=Dataset(
        download=download,
        mapper=enzyme_complexes_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='enzyme_complexes', **kwargs
        ),
    ),
    genes=Dataset(
        download=download,
        mapper=genes_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='genes', **kwargs
        ),
    ),
)
