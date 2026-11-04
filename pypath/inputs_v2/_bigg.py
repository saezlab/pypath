"""The schema factory every BiGG genome-scale model builds its four datasets from.

A BiGG model gives: metabolites (chemical participants, each with cross-
references to HMDB, ChEBI, KEGG and MetaNetX, a formula and a charge),
reactions (an SBO-typed molecular activity, with directed, compartmented,
stoichiometric participants), enzyme complexes (the multi-subunit AND groups
a reaction's gene rule carries), and genes (one record per Entrez ID). All
four read the BiGG JSON format every BiGG model ships, parsed generically by
:mod:`pypath.inputs_v2.parsers._bigg`.

What is true of every BiGG model so far (SBO-based reaction typing, the
Entrez gene namespace, the nine-compartment map, the Recon3D-style isoform
suffix pattern) stays a literal here — a model that genuinely needs
something else is new scope for whoever onboards it, not a parameter this
factory pre-builds (see ``pypath/inputs_v2/parsers/_bigg.py`` and
``pypath/inputs_v2/_source_context.py``). What differs between models —
taxon, the GPR group identifier prefix, and whether the standalone
``metabolites`` dataset excludes metabolites already covered as reaction
participants — are the three parameters :func:`bigg_resource` takes.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.reactions import COMPARTMENT
from omnipath_core.source_attributes import CELLULAR_LOCATION

from pypath.inputs_v2._reactions import (
    _gene_rule_group,
    model_reaction_compartments,
    reaction_compartments_cv,
    reaction_type_cv,
    source_reaction_type_cv,
    stoichiometry_compartment_go,
)
from pypath.inputs_v2._source_context import (
    compartment_name,
    conversion_direction_cv,
)
from pypath.inputs_v2.base import Dataset, Download
from pypath.inputs_v2.parsers._bigg import _raw
from pypath.internals.silver_schema import Entity, Membership
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
)


def field_config() -> FieldConfig:
    """The ``FieldConfig`` every BiGG model's schemas read.

    Holds no model-specific value: the stoichiometry/compartment/complex-name
    maps operate on the generic ``base_id:compartment:coefficient`` and
    ``||``-joined member-list encodings every BiGG model's rows carry.
    """
    return FieldConfig(
        delimiter=';',
        map={
            'stoich_id': lambda value: value.split(':')[0] if value else None,
            'stoich_comp': lambda value: compartment_name(value.split(':')[1])
            if value
            else None,
            'stoich_val': lambda value: value.split(':')[2] if value else None,
            'stoich_go': stoichiometry_compartment_go,
            'split_member_values': lambda value: value.split(';;')
            if value
            else None,
            'complex_name': lambda value: 'complex:' + '-'.join(value.split('||'))
            if value
            else '',
        },
    )


def metabolites_schema() -> EntityBuilder:
    """The standalone ``metabolites`` dataset's entity shape.

    Takes no knob: every BiGG model's metabolite JSON carries the same six
    cross-reference namespaces, a formula and a charge.
    """
    f = field_config()
    return EntityBuilder(
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


def _participant_members(f: FieldConfig, *, side: str, predicate) -> MembersFromList:
    # 'reactants'/'products' (the stoichiometry-bearing field) vs. the
    # singular 'reactant'/'product' prefix the per-member fields carry.
    prefix = 'reactant' if side == 'reactants' else 'product'
    return MembersFromList(
        entity_type=model.ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.BIGG_METABOLITE,
                value=f(side, delimiter='||', map='stoich_id', preserve_indices=True),
            ),
            CV(
                term=Namespace.HMDB,
                value=f(f'{prefix}_hmdb', delimiter='||', map='split_member_values', preserve_indices=True),
            ),
            CV(
                term=Namespace.CHEBI,
                value=f(f'{prefix}_chebi', delimiter='||', map='split_member_values', preserve_indices=True),
            ),
            CV(
                term=Namespace.KEGG,
                value=f(f'{prefix}_kegg_compound', delimiter='||', map='split_member_values', preserve_indices=True),
            ),
            CV(
                term=Namespace.METANETX,
                value=f(f'{prefix}_metanetx', delimiter='||', map='split_member_values', preserve_indices=True),
            ),
            CV(
                term=Namespace.NAME,
                value=f(f'{prefix}_name', delimiter='||', preserve_indices=True),
            ),
        ),
        annotations=AnnotationsBuilder(
            conversion_direction_cv(),
            CV(
                term=slots.stoichiometry,
                value=f(side, delimiter='||', map='stoich_val', preserve_indices=True),
            ),
            CV(term=COMPARTMENT, value=f(side, delimiter='||', map='stoich_go', preserve_indices=True)),
            CV(term=CELLULAR_LOCATION, value=f(side, delimiter='||', map='stoich_comp', preserve_indices=True)),
        ),
        entity_annotations=AnnotationsBuilder(
            CV(
                term='chemrof:charge',
                value=f(f'{prefix}_charge', delimiter='||', preserve_indices=True),
            ),
            CV(
                term=slots.has_chemical_formula,
                value=f(f'{prefix}_formula', delimiter='||', preserve_indices=True),
            ),
        ),
        predicate=predicate,
    )


def reactions_schema(*, taxon: str, gpr_prefix: str) -> EntityBuilder:
    """The ``reactions`` dataset's entity shape.

    Args:
        taxon: The ``NCBITaxon:...`` value every catalyst gene in a GPR
            alternative is annotated with.
        gpr_prefix: The model-scoped identifier type for a GPR AND-group
            entity (e.g. ``'recon3d.gpr_and'``) — must be unique per model.
    """
    f = field_config()

    def gene_rule_group(genes: list[str]) -> Entity:
        return _gene_rule_group(
            genes, identifier_type=gpr_prefix, taxon=taxon, gene_namespace=Namespace.ENTREZ,
        )

    class _ReactionMembership(MembershipBuilder):
        """Attach alternative GPR clauses, retaining required genes within each."""

        def build(self, row, cache):
            members = super().build(row, cache)
            for genes in row.get('gene_rule_clauses', []):
                if not genes:
                    continue
                group = gene_rule_group(genes)
                member = group.membership[0].member if len(genes) == 1 else group
                members.append(Membership(member=member, predicate=slots.associated_with))
            return members

    return EntityBuilder(
        entity_type=model.MolecularActivity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.BIGG_REACTION, value=f('bigg_reaction_id')),
            CV(term=Namespace.METANETX_REACTION, value=f('metanetx_reaction')),
            CV(term=Namespace.NAME, value=f('name')),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.has_topic,
                value=f('ec', transform=lambda v: 'EC:' + str(v).removeprefix('EC:')),
            ),
            reaction_type_cv(f('sbo')),
            # The model's own subsystem (e.g. "Glycolysis/gluconeogenesis").
            source_reaction_type_cv(f('subsystem')),
            reaction_compartments_cv(model_reaction_compartments),
        ),
        membership=_ReactionMembership(
            _participant_members(f, side='reactants', predicate=slots.has_input),
            _participant_members(f, side='products', predicate=slots.has_output),
        ),
    )


def enzyme_complexes_schema(*, taxon: str, gpr_prefix: str):
    """The ``enzyme_complexes`` dataset's mapper: a plain row-to-entity function.

    Args: see :func:`reactions_schema`.
    """

    def mapper(row: dict) -> Entity | None:
        genes = [gene for gene in str(row.get('complex_subunits') or '').split('||') if gene]
        if not genes:
            return None
        return _gene_rule_group(
            genes, identifier_type=gpr_prefix, taxon=taxon, gene_namespace=Namespace.ENTREZ,
        )

    return mapper


def genes_schema(*, taxon: str) -> EntityBuilder:
    """The ``genes`` dataset's entity shape.

    Args:
        taxon: The ``NCBITaxon:...`` value every gene entity is annotated
            with.
    """
    f = field_config()
    return EntityBuilder(
        entity_type=model.Gene,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.ENTREZ, value=f('entrez_id')),
            CV(term=Namespace.GENESYMBOL, value=f('name')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=taxon),
        ),
    )


def bigg_resource(
    *, download: Download, taxon: str, gpr_prefix: str, orphan_filter: bool = True,
) -> dict[str, Dataset]:
    """The four datasets a BiGG model builds, over its own ``download``.

    Args:
        download: The model's ``Download`` (its own BiGG JSON URL/filename);
            shared by all four datasets, as every BiGG model module shares
            one download across its datasets today.
        taxon: The ``NCBITaxon:...`` value for this model's genes and GPR
            catalysts. No default — every model differs, and a wrong default
            would silently mislabel an organism's genes.
        gpr_prefix: This model's GPR AND-group identifier prefix (e.g.
            ``'imm1415.gpr_and'``). No default — two models sharing a prefix
            would collide on an identical gene-rule digest.
        orphan_filter: Whether the ``metabolites`` dataset excludes a
            metabolite that also appears as a reaction participant.
            Defaults to ``True`` (Recon3D's behavior: the standalone dataset
            covers only metabolites no reaction mentions). A model with no
            `reactions` dataset of its own, or one that wants its standalone
            dataset to list every metabolite regardless of reaction
            membership, passes ``False``.

    Returns:
        A dict with keys ``metabolites``, ``reactions``, ``enzyme_complexes``
        and ``genes``, each a ``Dataset`` — not a ``Resource``, so a model
        that wants to replace one of the four (iMM1415 replaces
        ``metabolites`` with its own hand-written dataset) can do so with a
        plain dict merge before constructing its ``Resource``.
    """
    return {
        'metabolites': Dataset(
            download=download,
            mapper=metabolites_schema(),
            raw_parser=lambda opener, **kwargs: _raw(
                opener, data_type='metabolites', orphan_filter=orphan_filter, **kwargs,
            ),
            allow_empty=True,
        ),
        'reactions': Dataset(
            download=download,
            mapper=reactions_schema(taxon=taxon, gpr_prefix=gpr_prefix),
            raw_parser=lambda opener, **kwargs: _raw(opener, data_type='reactions', **kwargs),
        ),
        'enzyme_complexes': Dataset(
            download=download,
            mapper=enzyme_complexes_schema(taxon=taxon, gpr_prefix=gpr_prefix),
            raw_parser=lambda opener, **kwargs: _raw(opener, data_type='enzyme_complexes', **kwargs),
        ),
        'genes': Dataset(
            download=download,
            mapper=genes_schema(taxon=taxon),
            raw_parser=lambda opener, **kwargs: _raw(opener, data_type='genes', **kwargs),
        ),
    }
