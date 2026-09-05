"""Human-GEM chemicals and model reactions.

Input/output edges retain model orientation and stoichiometry. Numeric flux
bounds are source-defined quantities; compartments remain source context. Gene
rules are associations, not independent catalysis claims for every listed gene.
"""

from __future__ import annotations

import hashlib
import json

import math

from biolink_model.datamodel.model import (
    ChemicalEntity,
    Gene,
    MacromolecularComplex,
    MolecularActivity,
    NamedThing,
    QuantityValue,
    slots,
)
from omnipath_core.measurements import Measurement
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.metatlas import _raw
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
        type=NamedThing,
        identifiers=[Identifier(type='human_gem.gpr_and', value=digest)],
        membership=[
            Membership(
                member=Entity(
                    type=Gene,
                    identifiers=[Identifier(type=Namespace.ENSG, value=gene)],
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


def _flux_bound(row, field):
    """Preserve finite model bounds without inventing absent unit metadata."""
    value = row.get(field)
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return Measurement(
        quantity=QuantityValue(
            has_numeric_value=number,
            has_unit=row.get('flux_units') or None,
        ),
        source_field=field,
    )


config = ResourceConfig(
    id=ResourceCv.HUMAN_GEM,
    name='Human-GEM',
    url='https://github.com/SysBioChalmers/Human-GEM',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='41950094',
    primary_category='metabolic_network',
    description='Human-GEM is a genome-scale metabolic model of human metabolism maintained by the SysBioChalmers group.',
)
f = FieldConfig(
    delimiter=';',
    extract={'chebi': '^(?:CHEBI:)?(\\d+)$'},
    map={
        'entity_type': lambda value: MacromolecularComplex
        if value == 'complex'
        else Gene,
        'stoich_id': lambda value: value.split(':')[0] if value else None,
        'stoich_comp': lambda value: value.split(':')[1] if value else None,
        'stoich_val': lambda value: value.split(':')[2] if value else None,
        'complex_name': lambda value: 'complex:' + '-'.join(value.split('||'))
        if value
        else '',
    },
)
download = Download(
    url='https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/Human-GEM.yml',
    filename='Human-GEM.yml',
    subfolder='metatlas',
    large=True,
    ext='yml',
    default_mode='r',
)
metabolites_schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.HUMAN_GEM_METABOLITE,
            value=f('human_gem_metabolite_id'),
        ),
        CV(term=Namespace.HMDB, value=f('hmdb')),
        CV(term=Namespace.CHEBI, value=f('chebi', extract='chebi')),
        CV(term=Namespace.PUBCHEM, value=f('pubchem_compound')),
        CV(term=Namespace.LIPIDMAPS, value=f('lipidmaps')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('formula')),
        CV(term='chemrof:charge', value=f('charge')),
    ),
)
reactions_schema = EntityBuilder(
    entity_type=MolecularActivity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.HUMAN_GEM_REACTION, value=f('human_gem_reaction_id')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_quantitative_value,
            value=lambda row: _flux_bound(row, 'lower_bound'),
        ),
        CV(
            term=slots.has_quantitative_value,
            value=lambda row: _flux_bound(row, 'upper_bound'),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'eccodes',
                transform=lambda v: 'EC:' + str(v).removeprefix('EC:'),
            ),
        ),
    ),
    membership=_ReactionMembership(
        MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.HUMAN_GEM_METABOLITE,
                    value=f('reactants', delimiter='||', map='stoich_id'),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'reactant_hmdb', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'reactant_chebi',
                        delimiter='||',
                        extract='chebi',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.PUBCHEM,
                    value=f(
                        'reactant_pubchem_compound',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.LIPIDMAPS,
                    value=f(
                        'reactant_lipidmaps',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f('reactants', delimiter='||', map='stoich_val'),
                ),
                CV(
                    term='biopax:cellularLocation',
                    value=f('reactants', delimiter='||', map='stoich_comp'),
                )
            ),
            predicate=slots.has_input,
        ),
        MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.HUMAN_GEM_METABOLITE,
                    value=f('products', delimiter='||', map='stoich_id'),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'product_hmdb', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'product_chebi',
                        delimiter='||',
                        extract='chebi',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.PUBCHEM,
                    value=f(
                        'product_pubchem_compound',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.LIPIDMAPS,
                    value=f(
                        'product_lipidmaps',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f('products', delimiter='||', map='stoich_val'),
                ),
                CV(
                    term='biopax:cellularLocation',
                    value=f('products', delimiter='||', map='stoich_comp'),
                )
            ),
            predicate=slots.has_output,
        ),
    ),
)
transport_reactions_schema = EntityBuilder(
    entity_type=MolecularActivity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.HUMAN_GEM_REACTION, value=f('human_gem_reaction_id')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_quantitative_value,
            value=lambda row: _flux_bound(row, 'lower_bound'),
        ),
        CV(
            term=slots.has_quantitative_value,
            value=lambda row: _flux_bound(row, 'upper_bound'),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'eccodes',
                transform=lambda v: 'EC:' + str(v).removeprefix('EC:'),
            ),
        ),
    ),
    membership=_ReactionMembership(
        MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.HUMAN_GEM_METABOLITE,
                    value=f('reactants', delimiter='||', map='stoich_id'),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'reactant_hmdb', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'reactant_chebi',
                        delimiter='||',
                        extract='chebi',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.PUBCHEM,
                    value=f(
                        'reactant_pubchem_compound',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.LIPIDMAPS,
                    value=f(
                        'reactant_lipidmaps',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f('reactants', delimiter='||', map='stoich_val'),
                ),
                CV(
                    term='biopax:cellularLocation',
                    value=f('reactants', delimiter='||', map='stoich_comp'),
                )
            ),
            predicate=slots.has_input,
        ),
        MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.HUMAN_GEM_METABOLITE,
                    value=f('products', delimiter='||', map='stoich_id'),
                ),
                CV(
                    term=Namespace.HMDB,
                    value=f(
                        'product_hmdb', delimiter='||', preserve_indices=True
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'product_chebi',
                        delimiter='||',
                        extract='chebi',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.PUBCHEM,
                    value=f(
                        'product_pubchem_compound',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=Namespace.LIPIDMAPS,
                    value=f(
                        'product_lipidmaps',
                        delimiter='||',
                        preserve_indices=True,
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f('products', delimiter='||', map='stoich_val'),
                ),
                CV(
                    term='biopax:cellularLocation',
                    value=f('products', delimiter='||', map='stoich_comp'),
                )
            ),
            predicate=slots.has_output,
        ),
    ),
)


def enzyme_complexes_schema(row: dict) -> Entity | None:
    """Retain a source model AND group with a deterministic internal ID."""
    genes = [gene for gene in str(row.get('complex_subunits') or '').split('||') if gene]
    return _gene_rule_group(genes) if genes else None




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
)
