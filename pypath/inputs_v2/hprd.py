"""
Parse the HPRD flat files and emit enzyme-substrate, interaction and complex records.

HPRD, the Human Protein Reference Database, was curated from the literature by
the Pandey lab and the Institute of Bioinformatics. The last release, number 9
(2010), is distributed as one archive of flat files. Three of its tables are
read:

- ``enzyme_substrate``: the modification sites for which HPRD names an enzyme.
  Each becomes an enzyme-``affects``-substrate relation with the modification
  type, the site on the substrate isoform, the experiment types and the
  publications. Sites without an enzyme are not relations and are not emitted.
- ``interactions``: the binary protein-protein interactions, each an
  ``interacts_with`` relation with the experiment types and publications.
- ``complexes``: the protein complexes, each with its members.

``BINARY_PROTEIN_NONPROTEIN_INTERACTIONS.txt`` is not read: its partners are
free text ("DNA", "Calcuim", "17 al", "0") and are not identified.

Every protein is human: the database is the Human Protein Reference Database.
The archive has no species column, so the taxon is the database's, not a row's.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.molecular_forms import molecular_form_from_identifiers
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2._molecular_forms import combine_forms
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.hprd import (
    COMPLEX_MEMBER,
    INTERACTION_MEMBER,
    PTM_MEMBER,
    iter_complexes,
    iter_enzyme_sites,
    iter_interactions,
)
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import Membership
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.HPRD,
    name='HPRD',
    full='Human Protein Reference Database',
    url='http://www.hprd.org/',
    # HPRD's registration form for the download said: "I hereby declare that I
    # will not use or share the data contained in HPRD for any commercial
    # purposes and agree to not distribute the data to a third party". Its FAQ
    # said: "Everything in HPRD is free as long as it is not used for
    # commercial purposes. Commercial entities will have to pay a fee under a
    # licensing arrangement" (Internet Archive captures of hprd.org/download
    # and of the FAQ, 2023). hprd.org no longer answers (503).
    license=LicenseCV.NONCOMMERCIAL,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='18988627',
    topics=(Topic.SIGNALING,),
    content=(
        Content.ENZYME_PTM_RELATIONSHIPS,
        Content.PHYSICAL_INTERACTIONS,
        Content.PROTEIN_COMPLEXES,
    ),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human enzyme-substrate relations with modification sites, protein-'
        'protein interactions and protein complexes, with experiment types '
        'and publications, release 9 (2010).'
    ),
    description=(
        'The Human Protein Reference Database (HPRD) is a curated resource '
        'of human proteins, built by reading the published literature. For '
        'each protein it records the post-translational modifications with '
        'their sites and the enzymes responsible, the binary protein-protein '
        'interactions and the protein complexes, each with the type of '
        'experiment (in vivo, in vitro, yeast two-hybrid) and the PubMed '
        'references. The last release, number 9, is from 2010. The website '
        'is no longer served; OmniPath keeps a copy of the flat files.'
    ),
)

# hprd.org no longer answers (503). This is OmniPath's copy of the archive the
# site served as http://www.hprd.org/RELEASE9/HPRD_FLAT_FILES_041310.tar.gz,
# the file the legacy pypath module read.
_ARCHIVE_URL = 'https://rescued.omnipathdb.org/HPRD_FLAT_FILES_041310.tar.gz'


def _download(member: str) -> Download:
    return Download(
        url=_ARCHIVE_URL,
        filename='HPRD_FLAT_FILES_041310.tar.gz',
        subfolder='hprd',
        large=True,
        default_mode='rb',
        encoding=None,
        ext='tar.gz',
        needed=[member],
    )


# Identifier systems of HPRD without a Namespace member: its protein records
# ("00001"), the isoforms of those ("00001_1") and its complexes ("COM_1").
_HPRD = 'hprd'
_HPRD_ISOFORM = 'hprd_isoform'
_HPRD_COMPLEX = 'hprd_complex'

_TAXON = 'NCBITaxon:9606'

# The kind of accession of the "refseq" columns, as the parser names it. The
# nucleotide records (``genbank_nucleotide``) are not identifiers of a protein.
_ACCESSION_NAMESPACES = {
    'refseq_protein': Namespace.REFSEQ_PROTEIN,
    'genbank': Namespace.GENBANK,
    'uniprot': Namespace.UNIPROT,
}

# HPRD's experiment types, which it assigns to the kind of experiment (its own
# words: coimmunoprecipitation is "in vivo", a GST pull-down "in vitro"), and
# the Evidence & Conclusion Ontology term of each. The values are those seen
# in the 2010 archive; any other value stops the mapper.
_IN_VIVO = 'ECO:0000178'  # in vivo assay evidence
_IN_VITRO = 'ECO:0000181'  # in vitro assay evidence
_YEAST_2_HYBRID = 'ECO:0000068'  # yeast 2-hybrid evidence
EXPERIMENT_TYPES: dict[str, str | None] = {
    'in vivo': _IN_VIVO,
    'in vitro': _IN_VITRO,
    'yeast 2-hybrid': _YEAST_2_HYBRID,
    # Spellings of those three that occur in a few rows.
    'invivo': _IN_VIVO,
    'invitro': _IN_VITRO,
    'Y2H': _YEAST_2_HYBRID,
    # Codes of a few interaction rows (23 "VV", 6 "VT", ...). The archive does
    # not explain them; "in vivo", "in vitro" would be a guess, so they give
    # no evidence term and the row keeps its other types.
    'VV': None,
    'VT': None,
    'Vt': None,
}


def _evidence_terms(row: dict[str, Any]) -> list[str]:
    """ECO terms of a row's experiment types; an unknown type is an error."""
    terms = []
    for experiment_type in row['experiment_types']:
        if experiment_type not in EXPERIMENT_TYPES:
            raise ValueError(f'Unknown HPRD experiment type: {experiment_type!r}')
        term = EXPERIMENT_TYPES[experiment_type]
        if term is not None and term not in terms:
            terms.append(term)
    return terms


_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
_INCREASED = model.DirectionQualifierEnum.increased
_DECREASED = model.DirectionQualifierEnum.decreased

# HPRD's modification types, each with the aspect of the substrate the
# enzyme changes, the direction of the change and the PSI-MOD name of the
# modified residue where there is a faithful one. A writer increases the
# modification and an eraser (the "de-" types) decreases it. The direction is
# that of the modification; HPRD says nothing about the substrate's activity.
_MODIFICATION_TYPES = {
    'Phosphorylation': (_ASPECT.phosphorylation, _INCREASED, 'phosphorylated residue'),
    'Dephosphorylation': (_ASPECT.phosphorylation, _DECREASED, 'phosphorylated residue'),
    'Acetylation': (_ASPECT.acetylation, _INCREASED, 'acetylated residue'),
    'Deacetylation': (_ASPECT.acetylation, _DECREASED, 'acetylated residue'),
    'Methylation': (_ASPECT.methylation, _INCREASED, 'methylated residue'),
    'Demethylation': (_ASPECT.methylation, _DECREASED, 'methylated residue'),
    # PSI-MOD names ubiquitination and sumoylation for lysine only.
    'Ubiquitination': (_ASPECT.ubiquitination, _INCREASED, None),
    'Sumoylation': (_ASPECT.sumoylation, _INCREASED, None),
    'Glycosylation': (_ASPECT.glycosylation, _INCREASED, None),
    'Deglycosylation': (_ASPECT.glycosylation, _DECREASED, None),
    'Proteolytic Cleavage': (_ASPECT.cleavage, _INCREASED, None),
    'Palmitoylation': (_ASPECT.palmitoylation, _INCREASED, None),
    'Sulfation': (_ASPECT.sulfation, _INCREASED, None),
    'Glycation': (_ASPECT.glycation, _INCREASED, None),
    'Prenylation': (_ASPECT.prenylation, _INCREASED, None),
    'Carboxylation': (_ASPECT.carboxylation, _INCREASED, 'carboxylated residue'),
    'Myristoylation': (_ASPECT.myristoylation, _INCREASED, 'myristoylated residue'),
    'Farnesylation': (_ASPECT.farnesylation, _INCREASED, 'farnesylated residue'),
    'Glutathionylation': (_ASPECT.glutathionylation, _INCREASED, None),
    'Hydroxylation': (_ASPECT.hydroxylation, _INCREASED, None),
    'Alkylation': (_ASPECT.alkylation, _INCREASED, None),
    'ADP-Ribosylation': (getattr(_ASPECT, 'ADP-ribosylation'), _INCREASED, None),
    'Glucuronosyl transfer': (_ASPECT.glucuronidation, _INCREASED, None),
    'Deacylation': (_ASPECT.acylation, _DECREASED, None),
}
# Modification types that no aspect value of the pinned Biolink model names;
# the nearest values (nitrosation, lipidation, ubiquitination, amination) are
# broader or another modification. A row of one of them gives no relation.
MODIFICATION_TYPES_WITHOUT_ASPECT = frozenset(
    {
        'Amidation',
        'Transglutamination',
        'Neddylation',
        'Nitration',
        'S-Nitrosylation',
        'Disulfide Bridge',
        'Glycosyl phosphatidyl inositol GPI',
    }
)
# PSI-MOD names of ubiquitinated and sumoylated lysine.
_LYSINE_TERMS = {
    'Ubiquitination': 'ubiquitinylated lysine',
    'Sumoylation': 'sumoylated lysine',
}


def _modification_type(row: dict[str, Any]) -> str:
    value = row['modification_type']
    if value not in _MODIFICATION_TYPES and value not in MODIFICATION_TYPES_WITHOUT_ASPECT:
        raise ValueError(f'Unknown HPRD modification type: {value!r}')
    return value


f = FieldConfig(extract={'pmid': PMID_EXTRACT}, transform={'pmid': lambda v: f'PMID:{v}'})

_publications = CV(
    term=slots.publications,
    value=f(lambda row: row['references'], extract='pmid', transform='pmid'),
)
_evidence_types = CV(term=slots.has_evidence_of_type, value=_evidence_terms)


def _protein(select) -> Any:
    """Builder for a protein from ``select(row)``, a dict that has ``hprd_id``,
    ``gene_symbol``, ``accession``, ``accession_kind`` and, for a substrate,
    ``isoform_id`` and ``molecular_form``; None if ``select`` gives None."""

    def accession(kind: str):
        return lambda p: p['accession'] if p['accession_kind'] == kind else None

    builder = EntityBuilder(
        entity_type=model.Protein,
        molecular_form=lambda p: p.get('molecular_form'),
        identifiers=IdentifiersBuilder(
            *(
                CV(term=namespace, value=accession(kind))
                for kind, namespace in _ACCESSION_NAMESPACES.items()
            ),
            CV(term=_HPRD, value=lambda p: p['hprd_id']),
            CV(term=_HPRD_ISOFORM, value=lambda p: p.get('isoform_id')),
            CV(term=Namespace.GENESYMBOL, value=lambda p: p['gene_symbol']),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=lambda p: _TAXON),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        participant = select(row)
        return builder.build(participant) if participant else None

    return build


def _interactor(number: int):
    def select(row: dict[str, Any]) -> dict[str, Any]:
        return {
            'hprd_id': row[f'interactor_{number}_hprd_id'],
            'gene_symbol': row[f'interactor_{number}_gene_symbol'],
            'accession': row[f'interactor_{number}_accession'],
            'accession_kind': row[f'interactor_{number}_accession_kind'],
        }

    return select


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The modified residue, a feature of the substrate isoform in this record.

    The position is the source's, in the sequence of ``substrate_accession``;
    it is not checked against that sequence here. A RefSeq accession names
    the isoform's sequence; a GenBank or UniProt accession is kept as the
    coordinate reference only.
    """
    kind = row['substrate_accession_kind']
    # A nucleotide record in the protein column is no protein sequence.
    reference = (
        {'ns': str(_ACCESSION_NAMESPACES[kind]), 'id': row['substrate_accession']}
        if kind in _ACCESSION_NAMESPACES
        else None
    )
    modification_type = row['modification_type']
    _, _, term = _MODIFICATION_TYPES[modification_type]
    if term is None and row['residue'] == 'K':
        term = _LYSINE_TERMS.get(modification_type)
    # A site written "0" is no position of a sequence: the residue is kept
    # without one, and the cell as the description.
    position = row['position'] or None
    description = modification_type if term is None else None
    if position is None:
        description = f'{description or modification_type} (site {row["site"]})'
    site = {
        'modifications': [
            {
                'term': term,
                'residue': row['residue'],
                'position': position,
                'end_position': position,
                'coordinate_reference': {
                    'identifier': reference,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                # Without a PSI-MOD name the entry keeps HPRD's type.
                'description': description,
            }
        ],
    }
    return combine_forms(
        molecular_form_from_identifiers([reference] if reference else []), site
    )


def _substrate(row: dict[str, Any]) -> dict[str, Any]:
    return {
        'hprd_id': row['substrate_hprd_id'],
        'isoform_id': row['substrate_isoform_id'],
        'gene_symbol': row['substrate_gene_symbol'],
        'accession': row['substrate_accession'],
        'accession_kind': row['substrate_accession_kind'],
        'molecular_form': _substrate_form(row),
    }


def _enzyme(row: dict[str, Any]) -> dict[str, Any]:
    # HPRD names the enzyme by gene symbol and HPRD identifier only.
    return {
        'hprd_id': row['enzyme_hprd_id'],
        'gene_symbol': row['enzyme_name'],
        'accession': None,
        'accession_kind': None,
    }


_enzyme_substrate_schema = RelationBuilder(
    subject=_protein(_enzyme),
    predicate=slots.affects,
    object=_protein(_substrate),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: _MODIFICATION_TYPES[row['modification_type']][1],
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=lambda row: _MODIFICATION_TYPES[row['modification_type']][0],
        ),
        _evidence_types,
        _publications,
    ),
)


def map_enzyme_substrate(row: dict[str, Any]) -> Any:
    """Map a row that names an enzyme; None if its modification type has no
    Biolink aspect.

    An unknown modification type or experiment type raises ``ValueError``, as
    does a second position (a pair of residues, as in a disulfide bridge),
    which no row that names an enzyme has.
    """
    if _modification_type(row) in MODIFICATION_TYPES_WITHOUT_ASPECT:
        return None
    _evidence_terms(row)
    if row['second_position'] is not None:
        raise ValueError(f'HPRD site with two positions: {row["site"]!r}')
    return _enzyme_substrate_schema(row)


_interactions_schema = RelationBuilder(
    subject=_protein(_interactor(1)),
    # HPRD records physical interactions and does not say which are direct.
    predicate=slots.interacts_with,
    object=_protein(_interactor(2)),
    annotations=AnnotationsBuilder(_evidence_types, _publications),
)


def map_interaction(row: dict[str, Any]) -> Any:
    """Map one row of the binary interactions; an unknown experiment type raises."""
    _evidence_terms(row)
    return _interactions_schema(row)


_complex_schema = EntityBuilder(
    entity_type=model.MacromolecularComplex,
    identifiers=IdentifiersBuilder(
        CV(term=_HPRD_COMPLEX, value=lambda row: row['complex_id']),
    ),
    annotations=AnnotationsBuilder(_evidence_types, _publications),
)
_complex_member = _protein(lambda member: member)


def map_complex(row: dict[str, Any]) -> Any:
    """Map a complex of at least two identified members.

    A complex with a member that has no HPRD identifier is not mapped: the
    source states a member that cannot be identified, and the others alone
    would misstate the composition. HPRD gives no name and no stoichiometry.
    """
    _evidence_terms(row)
    if row['unidentified_rows'] or len(row['members']) < 2:
        return None
    head = _complex_schema(row)
    return head._replace(
        membership=[
            Membership(member=_complex_member(member)) for member in row['members']
        ]
    )


resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=_download(PTM_MEMBER),
        mapper=map_enzyme_substrate,
        raw_parser=iter_enzyme_sites,
    ),
    interactions=Dataset(
        download=_download(INTERACTION_MEMBER),
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
    complexes=Dataset(
        download=_download(COMPLEX_MEMBER),
        mapper=map_complex,
        raw_parser=iter_complexes,
    ),
)
