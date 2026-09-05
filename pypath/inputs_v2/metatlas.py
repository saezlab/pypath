"""Human-GEM chemicals and model reactions.

Input/output edges retain model orientation and stoichiometry. Numeric flux
bounds are source-defined quantities; compartments remain source context. Gene
rules are associations, not independent catalysis claims for every listed gene.
"""

from __future__ import annotations

import math

from biolink_model.datamodel.model import (
    BiologicalProcess,
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
from pypath.internals.silver_schema import Entity, Identifier, Membership
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
)


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
    membership=MembershipBuilder(
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
                )
            ),
            predicate=slots.has_output,
        ),
        MembersFromList(
            entity_type=Gene,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.ENSG, value=f('enzyme_ensembl')),
                CV(
                    term=Namespace.GENESYMBOL,
                    value=f('enzyme_name', preserve_indices=True),
                ),
            ),
            annotations=AnnotationsBuilder(),
            predicate=slots.associated_with,
        ),
    ),
)
transport_reactions_schema = EntityBuilder(
    entity_type=BiologicalProcess,
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
    membership=MembershipBuilder(
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
                )
            ),
            predicate=slots.has_output,
        ),
        MembersFromList(
            entity_type=Gene,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.ENSG, value=f('enzyme_ensembl')),
                CV(
                    term=Namespace.GENESYMBOL,
                    value=f('enzyme_name', preserve_indices=True),
                ),
            ),
            annotations=AnnotationsBuilder(),
            predicate=slots.associated_with,
        ),
    ),
)


def enzyme_complexes_schema(row):
    """An AND clause groups required genes without asserting a physical assembly."""
    members = [
        Entity(
            type=Gene, identifiers=[Identifier(type=Namespace.ENSG, value=gene)]
        )
        for gene in str(row.get('complex_subunits') or '').split('||')
        if gene
    ]
    if not members:
        return None
    return Entity(
        type=NamedThing,
        identifiers=[],
        membership=[
            Membership(member=member, predicate=slots.has_member)
            for member in members
        ],
    )


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
