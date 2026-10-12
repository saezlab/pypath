"""
Parse the scConnect ligand-receptor interactions and emit Relation records.

scConnect predicts cell-cell communication from the interaction table of the
Guide to Pharmacology (release 2020.5), which it ships in a cleaned copy. Of
that table scConnect uses the ligands it can score from gene expression:
peptides encoded by a gene, and eight neurotransmitters. Each row is one
measurement of a ligand on a target of one species. The table's `action` gives
the ligand's effect on the target (agonist, antagonist, inhibition ...); it is
carried as in `guidetopharma.py`.
"""

from __future__ import annotations
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    CausalMechanismQualifierEnum,
    ChemicalEntity,
    DirectionQualifierEnum,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.scconnect import (
    iter_interactions,
    ligand_gene,
    pubmed_ids,
    split_values,
)

config = ResourceConfig(
    id=ResourceCv.SCCONNECT,
    name='scConnect',
    url='https://github.com/JonETJakobsson/scConnect',
    # The tables are derived from the Guide to Pharmacology, which is under
    # CC BY-SA 4.0; the MIT license of the repository, which the legacy
    # registry recorded, covers the code.
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='33974001',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Peptide and neurotransmitter ligands with their receptors and the action of the ligand, from the Guide to Pharmacology tables of scConnect.',
    description='scConnect is a Python package that infers ligand-receptor connectivity between cell populations in single-cell RNA sequencing data. Its interactions come from the Guide to Pharmacology. Peptide ligands are scored by the expression of their gene, and eight neurotransmitters by the expression of the enzymes and transporters that make, load and take up the molecule.',
)
download = Download(
    url='https://raw.githubusercontent.com/JonETJakobsson/scConnect/master/scConnect/data/GTP_tables_clean/2020-5/interactions.csv',
    filename='scconnect_gtp_interactions_2020-5.csv',
    subfolder='scconnect',
)
# `target_species` of the targets that carry gene identifiers.
SPECIES = {
    'Human': '9606',
    'Mouse': '10090',
    'Rat': '10116',
}
# The small molecules scConnect scores as ligands, in the spelling of the
# table (`scConnect/molecules.py`). The molecule itself is the ligand; the
# enzymes and transporters scConnect lists for it are not.
MOLECULE_LIGANDS = frozenset(
    {
        'L-glutamic acid',
        'GABA',
        'glycine',
        'dopamine',
        '(-)-noradrenaline',
        '(-)-adrenaline',
        '5-hydroxytryptamine',
        'acetylcholine',
    }
)
# The interaction table is the Guide to Pharmacology's.
DATA_SOURCE = 'Guide to Pharmacology'
# The `action` column, as the Guide to Pharmacology writes it, to the causal
# mechanism and the direction of the ligand's effect on the target. Every value
# of the file is listed; an unlisted one raises. The mapping is that of
# `guidetopharma.py` (`_ACTION_MECHANISM`, `_ACTION_DIRECTION`), where the same
# values come from the same database: a value with a direction or a mechanism
# other than binding is an `affects`; `Binding` and the values that state no
# effect are an `interacts_with`. Broader qualifiers ("Competitive",
# "Non-competitive") do not alone determine inhibition and carry no qualifier,
# nor do "None", the file's word for no stated action, and "Unknown".
_mechanism = CausalMechanismQualifierEnum
_direction = DirectionQualifierEnum
ACTIONS = {
    'Activation': (_mechanism.activation, _direction.increased),
    'Agonist': (_mechanism.agonism, _direction.increased),
    'Full agonist': (_mechanism.agonism, _direction.increased),
    'Partial agonist': (_mechanism.partial_agonism, _direction.increased),
    'Biased agonist': (_mechanism.biased_agonism, None),
    'Antagonist': (_mechanism.antagonism, _direction.decreased),
    'Inhibition': (_mechanism.inhibition, _direction.decreased),
    'Feedback inhibition': (_mechanism.feedback_inhibition, _direction.decreased),
    'Irreversible inhibition': (
        _mechanism.irreversible_inhibition,
        _direction.decreased,
    ),
    'Inverse agonist': (_mechanism.inverse_agonism, _direction.decreased),
    'Voltage-dependent inhibition': (None, _direction.decreased),
    'Positive': (_mechanism.positive_modulation, None),
    'Negative': (_mechanism.negative_modulation, None),
    'Potentiation': (_mechanism.potentiation, None),
    'Pore blocker': (_mechanism.molecular_channel_blockage, None),
    'Binding': (_mechanism.binding, None),
    'Mixed': (_mechanism.modulation, None),
    'Biphasic': (_mechanism.modulation, None),
    'Competitive': (None, None),
    'Non-competitive': (None, None),
    'Neutral': (None, None),
    'Slows inactivation': (None, None),
    'Irreversible agonist': (None, None),
    'Unknown': (None, None),
    'None': (None, None),
}


def exclusion(row: dict[str, Any]) -> str | None:
    """Why a row yields no relation, or None if it yields one."""
    ligand_genes = split_values(row.get('ligand_gene_symbol'))
    if not ligand_genes and row.get('ligand') not in MOLECULE_LIGANDS:
        return 'ligand is not a peptide with a gene or a scConnect molecule'
    target_genes = split_values(row.get('target_gene_symbol'))
    if not target_genes:
        return 'target has no gene'
    if len(target_genes) > 1:
        return 'target is named by several genes'
    if row.get('target_species') not in SPECIES:
        return 'target species has no taxon'
    if not ligand_genes:
        return None
    species = split_values(row.get('ligand_species'))
    if row.get('target_species') not in species:
        return 'ligand is from another species than the target'
    if len(species) == 1 and len(ligand_genes) > 1:
        return 'ligand is named by several genes'
    if ligand_gene(row) is None:
        return 'ligand gene of the target species is not identified'
    return None


def _identifiers(*items: tuple[Any, Any]) -> list[Identifier]:
    return [
        Identifier(type=namespace, value=text)
        for namespace, value in items
        if (text := str(value or '').strip())
    ]


# The table has a `ligand` column. Its other side is the `target`, which the
# Guide to Pharmacology does not call a receptor (targets include enzymes and
# transporters), so the target gets no role.
def _ligand_role() -> Annotation:
    return Annotation(term=PARTICIPANT_ROLE, value='ligand')


def _ligand(row: dict[str, Any], taxon: str) -> Entity:
    gene = ligand_gene(row)
    if gene is None:
        # A neurotransmitter: a chemical, which has no species.
        return Entity(
            type=ChemicalEntity,
            identifiers=_identifiers((Namespace.NAME, row.get('ligand'))),
            annotations=[_ligand_role()],
        )
    # The gene encodes the peptide's precursor, which may yield several
    # ligands, so the peptide's name is not a name of the protein.
    return Entity(
        type=Protein,
        identifiers=_identifiers((Namespace.GENESYMBOL, gene)),
        annotations=[
            Annotation(term=slots.in_taxon, value=taxon),
            _ligand_role(),
        ],
    )


def map_scconnect_interaction(row: dict[str, Any]) -> Relation | None:
    # A ligand or a target named by several genes may be a complex, a set of
    # alternative subunits or one protein with several genes. The table does
    # not say which, so no complex is asserted and such rows are left out.
    # The action is read first, so that a value not in the table raises for
    # every row, not only for the rows that are emitted.
    action = row.get('action')
    if action not in ACTIONS:
        raise ValueError(f'Unknown scConnect action {action!r}')
    mechanism, direction = ACTIONS[action]
    if exclusion(row) is not None:
        return None
    taxon = f'NCBITaxon:{SPECIES[row["target_species"]]}'
    # The ligand is the subject and the receptor the object.
    return Relation(
        subject=_ligand(row, taxon),
        predicate=(
            slots.affects
            if direction is not None
            or (mechanism is not None and action != 'Binding')
            else slots.interacts_with
        ),
        object=Entity(
            type=Protein,
            identifiers=_identifiers(
                (Namespace.GENESYMBOL, row.get('target_gene_symbol')),
                (Namespace.UNIPROT, row.get('target_uniprot')),
                (Namespace.ENSG, row.get('target_ensembl_gene_id')),
                (Namespace.NAME, row.get('target')),
            ),
            annotations=[Annotation(term=slots.in_taxon, value=taxon)],
        ),
        identifiers=[],
        annotations=[
            Annotation(term=slots.supporting_data_source, value=DATA_SOURCE),
            *(
                [Annotation(term=slots.causal_mechanism_qualifier, value=mechanism)]
                if mechanism is not None
                else []
            ),
            *(
                [Annotation(term=slots.object_direction_qualifier, value=direction)]
                if direction is not None
                else []
            ),
            *[
                Annotation(term=slots.publications, value=f'PMID:{pmid}')
                for pmid in pubmed_ids(row.get('pubmed_id'))
            ],
        ],
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_scconnect_interaction,
        raw_parser=iter_interactions,
    ),
)
