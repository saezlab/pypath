"""Metabolic Atlas chemicals and model reactions: Human-GEM and the animal GEMs.

Human-GEM and the models derived from it for mouse, rat, zebrafish, fruit fly
and worm share Metabolic Atlas identifiers: ``MAM`` for metabolites and ``MAR``
for reactions name the same molecule and the same reaction in every model, so
each is one entity. What a model adds is its own: its genes (Ensembl for
Human-GEM, symbols for the others, each with the model's taxon) in the gene
rules of its reactions, its flux bounds, and its dataset. Human-GEM keeps the
dataset names ``metabolites``, ``reactions`` and ``enzyme_complexes``; the
other models prefix theirs (``mouse_reactions``, ...).


Each reaction is typed in SBO from its subsystem: "Transport reactions" are
transport, "Exchange/demand reactions" exchange, demand or sink, the rest
biochemical; the subsystem itself is kept beside the type, and a transport's
TCDB classes as topics. Input/output evidence retains model orientation,
stoichiometry, the GO compartment (as a Biolink anatomical context, with its
name) and reported conversion direction; the reaction states the compartments
it spans, which keeps its compartment variants apart. Numeric flux bounds are source-defined quantities.
Gene rules are associations, not independent catalysis claims for every listed gene.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from functools import partial

import math

from biolink_model.datamodel.model import (
    ChemicalEntity,
    Gene,
    MacromolecularComplex,
    MolecularActivity,
    QuantityValue,
    slots,
)
from omnipath_core.measurements import Measurement
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic
from omnipath_core.reactions import COMPARTMENT
from omnipath_core.source_attributes import CELLULAR_LOCATION

from omnipath_resources.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from omnipath_resources.inputs_v2._reactions import _gene_rule_group as _shared_gene_rule_group
from omnipath_resources.inputs_v2._reactions import (
    model_reaction_compartments,
    model_reaction_type,
    reaction_compartments_cv,
    reaction_type_cv,
    source_reaction_type_cv,
    stoichiometry_compartment_go,
)
from omnipath_resources.inputs_v2._source_context import compartment_name, conversion_direction_cv
from omnipath_resources.inputs_v2.parsers.metatlas import (
    HUMAN_GEM,
    _raw,
    model_metabolite_xrefs,
    model_url,
)
from omnipath_resources.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from omnipath_resources.internals.silver_schema import Annotation, Entity, Identifier, Membership
from omnipath_resources.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
)


@dataclass(frozen=True)
class GemModel:
    """A Metabolic Atlas model: its repository name, organism and gene identifiers."""

    name: str
    taxon: str
    gene_namespace: Namespace
    # Human-GEM's datasets keep their names; the other models' are prefixed.
    dataset_prefix: str

    @property
    def gene_rule_namespace(self) -> str:
        """The AND clauses of one model's gene rules; never merged across models."""
        return self.name.lower().replace('-', '_') + '.gpr_and'


MODELS = (
    GemModel(HUMAN_GEM, '9606', Namespace.ENSG, ''),
    GemModel('Mouse-GEM', '10090', Namespace.GENESYMBOL, 'mouse_'),
    GemModel('Rat-GEM', '10116', Namespace.GENESYMBOL, 'rat_'),
    GemModel('Zebrafish-GEM', '7955', Namespace.GENESYMBOL, 'zebrafish_'),
    GemModel('Fruitfly-GEM', '7227', Namespace.GENESYMBOL, 'fruitfly_'),
    GemModel('Worm-GEM', '6239', Namespace.GENESYMBOL, 'worm_'),
)
HUMAN = MODELS[0]


def _gene_rule_group(genes: list[str], model: GemModel = HUMAN) -> Entity:
    """Source-scoped logical AND clause, without a physical assembly claim."""
    return _shared_gene_rule_group(
        genes,
        identifier_type=model.gene_rule_namespace,
        taxon=f'NCBITaxon:{model.taxon}',
        gene_namespace=model.gene_namespace,
    )


class _ReactionMembership(MembershipBuilder):
    """Attach alternative GPR clauses, retaining required genes within each."""

    def __init__(self, model: GemModel, *members) -> None:
        super().__init__(*members)
        self.model = model

    def build(self, row, cache):
        members = super().build(row, cache)
        for genes in row.get('gene_rule_clauses', []):
            if not genes:
                continue
            group = _gene_rule_group(genes, self.model)
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
    # ±1000 is the model's default for an unconstrained flux, not a measured bound.
    if not math.isfinite(number) or abs(number) == 1000:
        return None
    return Measurement(
        quantity=QuantityValue(
            has_numeric_value=number,
            has_unit=row.get('flux_units') or None,
        ),
        source_field=field,
    )


_TRANSPORT_SUBSYSTEM = 'Transport reactions'
_PSEUDOREACTION_SUBSYSTEM = 'Exchange/demand reactions'


def _reaction_type(row) -> str:
    subsystems = {s.strip() for s in str(row.get('subsystem') or '').split(';')}
    return model_reaction_type(
        row,
        transport=_TRANSPORT_SUBSYSTEM in subsystems,
        pseudoreaction=_PSEUDOREACTION_SUBSYSTEM in subsystems,
    )


config = ResourceConfig(
    id=ResourceCv.HUMAN_GEM,
    name='Metabolic Atlas',
    full='Metabolic Atlas genome-scale metabolic models',
    url='https://metabolicatlas.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='41950094',
    topics=(Topic.METABOLISM,),
    content=(Content.BIOCHEMICAL_REACTIONS,),
    summary=(
        'Chemicals and metabolic reactions with stoichiometry and gene associations '
        'of Human-GEM and the mouse, rat, zebrafish, fruit fly and worm GEMs.'
    ),
    description=(
        'The genome-scale metabolic models of Metabolic Atlas maintained by the '
        'SysBioChalmers group: Human-GEM and the Mouse-, Rat-, Zebrafish-, '
        'Fruitfly- and Worm-GEM derived from it.'
    ),
)
f = FieldConfig(
    delimiter=';',
    extract={'chebi': '^(?:CHEBI:)?(\\d+)$'},
    map={
        'entity_type': lambda value: MacromolecularComplex
        if value == 'complex'
        else Gene,
        'stoich_id': lambda value: value.split(':')[0] if value else None,
        'stoich_comp': lambda value: compartment_name(value.split(':')[1])
        if value
        else None,
        'stoich_val': lambda value: value.split(':')[2] if value else None,
        'stoich_go': stoichiometry_compartment_go,
        'complex_name': lambda value: 'complex:' + '-'.join(value.split('||'))
        if value
        else '',
    },
)
def _download(model: GemModel) -> Download:
    return Download(
        url=model_url(model.name, f'model/{model.name}.yml'),
        filename=f'{model.name}.yml',
        subfolder='metatlas',
        large=True,
        ext='yml',
        default_mode='r',
    )


download = _download(HUMAN)
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
        CV(term=Namespace.KEGG, value=f('kegg_compound')),
        CV(term=Namespace.BIGG_METABOLITE, value=f('bigg')),
        CV(term=Namespace.METANETX, value=f('metanetx')),
        CV(term=Namespace.SMILES, value=f('smiles')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('formula')),
        CV(term='chemrof:charge', value=f('charge')),
    ),
)
def _reactions_schema(model: GemModel) -> EntityBuilder:
    return EntityBuilder(
        entity_type=MolecularActivity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.HUMAN_GEM_REACTION, value=f('human_gem_reaction_id')),
            CV(term=Namespace.RHEA, value=f('rhea')),
            CV(term=Namespace.KEGG_REACTION, value=f('kegg_reaction')),
            CV(term=Namespace.BIGG_REACTION, value=f('bigg_reaction')),
            CV(term=Namespace.METANETX_REACTION, value=f('metanetx_reaction')),
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
            CV(term=slots.has_topic, value=f('tcdb', transform=lambda v: 'TC:' + str(v))),
            reaction_type_cv(_reaction_type),
            source_reaction_type_cv(f('subsystem')),
            reaction_compartments_cv(model_reaction_compartments),
        ),
        membership=_ReactionMembership(
            model,
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
                    CV(
                        term=Namespace.KEGG,
                        value=f(
                            'reactant_kegg_compound',
                            delimiter='||',
                            preserve_indices=True,
                        ),
                    ),
                    CV(
                        term=Namespace.BIGG_METABOLITE,
                        value=f(
                            'reactant_bigg',
                            delimiter='||',
                            preserve_indices=True,
                        ),
                    ),
                    CV(
                        term=Namespace.METANETX,
                        value=f(
                            'reactant_metanetx',
                            delimiter='||',
                            preserve_indices=True,
                        ),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    conversion_direction_cv(),
                    CV(
                        term=slots.stoichiometry,
                        value=f('reactants', delimiter='||', map='stoich_val'),
                    ),
                    CV(
                        term=COMPARTMENT,
                        value=f('reactants', delimiter='||', map='stoich_go'),
                    ),
                    CV(
                        term=CELLULAR_LOCATION,
                        value=f('reactants', delimiter='||', map='stoich_comp'),
                    ),
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
                    CV(
                        term=Namespace.KEGG,
                        value=f(
                            'product_kegg_compound',
                            delimiter='||',
                            preserve_indices=True,
                        ),
                    ),
                    CV(
                        term=Namespace.BIGG_METABOLITE,
                        value=f(
                            'product_bigg',
                            delimiter='||',
                            preserve_indices=True,
                        ),
                    ),
                    CV(
                        term=Namespace.METANETX,
                        value=f(
                            'product_metanetx',
                            delimiter='||',
                            preserve_indices=True,
                        ),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    conversion_direction_cv(),
                    CV(
                        term=slots.stoichiometry,
                        value=f('products', delimiter='||', map='stoich_val'),
                    ),
                    CV(
                        term=COMPARTMENT,
                        value=f('products', delimiter='||', map='stoich_go'),
                    ),
                    CV(
                        term=CELLULAR_LOCATION,
                        value=f('products', delimiter='||', map='stoich_comp'),
                    ),
                ),
                predicate=slots.has_output,
            ),
        ),
    )


reactions_schema = _reactions_schema(HUMAN)


def _enzyme_complexes_schema(row: dict, model: GemModel = HUMAN) -> Entity | None:
    """Retain a source model AND group with a deterministic internal ID."""
    genes = [gene for gene in str(row.get('complex_subunits') or '').split('||') if gene]
    return _gene_rule_group(genes, model) if genes else None


def enzyme_complexes_schema(row: dict) -> Entity | None:
    """Human-GEM's AND groups."""
    return _enzyme_complexes_schema(row, HUMAN)




def _map_with_structures(schema, row, model: GemModel = HUMAN):
    """Join the model's native metabolite structures, including cached participants."""
    entity = schema(row)
    xrefs = model_metabolite_xrefs(model.name)
    structure_source = model_url(model.name, 'model/metabolites.tsv')

    def enrich(current):
        if current is None:
            return None
        identifiers = list(current.identifiers)
        annotations = list(current.annotations or [])
        native = next((i.value for i in identifiers if i.type == Namespace.HUMAN_GEM_METABOLITE), None)
        if native and not any(i.type == Namespace.SMILES for i in identifiers):
            smiles = xrefs.get(native, {}).get('smiles')
            if smiles:
                identifiers.append(Identifier(type=Namespace.SMILES, value=smiles))
                annotations.append(Annotation(
                    term='omnipath:structure_source', value=structure_source,
                ))
        members = [m._replace(member=enrich(m.member)) for m in current.membership or []]
        return current._replace(identifiers=identifiers, annotations=annotations, membership=members)

    return enrich(entity)


def preparation_inputs():
    """Expose the selected native auxiliary tables for preparation cache identity."""
    from omnipath_resources.share.downloads import _resolve_data_dir
    path = os.environ.get('OMNIPATH_HUMAN_GEM_XREFS')
    data = _resolve_data_dir() / 'metatlas'
    return [Path(path) if path else data / 'Human-GEM-metabolites.tsv'] + [
        data / f'{model.name}-metabolites.tsv' for model in MODELS if model is not HUMAN
    ]


def _datasets(model: GemModel) -> dict[str, Dataset]:
    model_download = download if model is HUMAN else _download(model)
    reactions = reactions_schema if model is HUMAN else _reactions_schema(model)
    complexes = enzyme_complexes_schema if model is HUMAN else partial(_enzyme_complexes_schema, model=model)
    return {
        f'{model.dataset_prefix}metabolites': Dataset(
            download=model_download,
            mapper=partial(_map_with_structures, metabolites_schema, model=model),
            raw_parser=partial(_raw, data_type='metabolites', model=model.name),
        ),
        f'{model.dataset_prefix}reactions': Dataset(
            download=model_download,
            mapper=partial(_map_with_structures, reactions, model=model),
            raw_parser=partial(_raw, data_type='reactions', model=model.name),
        ),
        f'{model.dataset_prefix}enzyme_complexes': Dataset(
            download=model_download,
            mapper=complexes,
            raw_parser=partial(_raw, data_type='enzyme_complexes', model=model.name),
        ),
    }


resource = Resource(
    config,
    **{name: dataset for model in MODELS for name, dataset in _datasets(model).items()},
)
