"""
Membranes and classification of the protein structures of OPM.

OPM, the Orientations of Proteins in Membranes database, places the
experimental structures of membrane proteins and peptides in the lipid bilayer
(by the PPM method) and classifies them. Each of its primary structures is one
PDB entry with the membrane it belongs to, the family in OPM's own
classification (type > class > superfamily > family) and the calculated
orientation in the membrane: the hydrophobic thickness, the tilt angle, the
transfer energy and the radius of curvature.

The module imports:

* ``terms``: the membranes and the four levels of the classification, as OPM
  numbers and names them (the tables have no ontology identifiers), each level
  in a namespace of its own, with the classification hierarchy as ``subclass_of``;
* ``membrane_localizations``: a statement that the structure is in a membrane, with the
  orientation of the structure in it;
* ``family_annotations``: a statement that the structure belongs to a family.

The structure is the subject, because that is what OPM classifies. A structure
is identified by what OPM gives and nothing is translated: its UniProt code
(an entry name such as ``OMPA_ECOLI``) when the entry has one protein, the PDB
code and the UniProt codes of the proteins as members when it has several, the
PDB code alone when OPM states no UniProt code. The legacy input looked the
protein names up in UniProt and kept the first match; the species of a
structure is a name in OPM and is not asserted.
"""

from __future__ import annotations

from functools import partial
import re
from typing import Any

from biolink_model.datamodel.model import (
    MacromolecularComplex,
    OntologyClass,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.opm import iter_structures, iter_terms
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    Membership,
    OntologyRelation,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.OPM,
    name='OPM',
    full='Orientations of Proteins in Membranes',
    url='https://opm.phar.umich.edu/',
    # License as recorded in the legacy registry (resources.json); the source
    # states none today (the home page and the script of its single page
    # application, checked 2026-10-12, carry no terms of use).
    license=LicenseCV.APACHE_2_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='21890895',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Membranes and families of the membrane protein structures of the PDB, '
        'with their calculated orientation in the membrane.'
    ),
    description=(
        'OPM (Orientations of Proteins in Membranes) is a database of the '
        'spatial positions of transmembrane, monotopic and peripheral '
        'proteins and peptides of known three-dimensional structure in the '
        'lipid bilayer, calculated with the PPM method. Its structures are '
        'classified by type, class, superfamily and family, and assigned to '
        'the biological membrane they belong to.'
    ),
)

# ================================== DOWNLOAD ==================================

# The legacy URLs (lomize-group-opm.herokuapp.com) are gone; the site's page
# application reads its tables from this backend.
_API = 'https://opm-back.cc.lehigh.edu/opm-backend'


def _table(name: str, endpoint: str) -> Download:
    return Download(
        url=f'{_API}/{endpoint}?fileFormat=csv',
        filename=f'opm_{name}.csv',
        subfolder='opm',
        large=True,
    )


download_membranes = _table('membranes', 'membranes')
download_types = _table('types', 'types')
download_classes = _table('classes', 'classtypes')
download_superfamilies = _table('superfamilies', 'superfamilies')
download_families = _table('families', 'families')
download_structures = _table('primary_structures', 'primary_structures')

# =================================== SCHEMA ===================================

# OPM numbers each table on its own, so each level has its own namespace.
MEMBRANE = 'opm_membrane'
TYPE = 'opm_type'
CLASS = 'opm_class'
SUPERFAMILY = 'opm_superfamily'
FAMILY = 'opm_family'
TERM_NAMESPACES = {
    'membrane': MEMBRANE,
    'type': TYPE,
    'class': CLASS,
    'superfamily': SUPERFAMILY,
    'family': FAMILY,
}

# Attributes of the vocabularies and of the statements; Biolink has no slot for
# them. The sides of a membrane are named as OPM writes them (``cytoplasmic
# side``, ``periplasm``), the first the inside of the topology OPM draws.
TOPOLOGY_IN = 'opm:topology_in'
TOPOLOGY_OUT = 'opm:topology_out'
PDB_ID = 'opm:pdb_id'
# Orientation of the structure in the membrane, as the site's tables label
# them: thickness and its error in angstroms, the tilt angle and its error in
# degrees, the transfer energy in kcal/mol, the radius of curvature in
# angstroms. The error columns carry the unit of their value.
HYDROPHOBIC_THICKNESS = 'opm:hydrophobic_thickness'
HYDROPHOBIC_THICKNESS_ERROR = 'opm:hydrophobic_thickness_error'
TILT_ANGLE = 'opm:tilt_angle'
TILT_ANGLE_ERROR = 'opm:tilt_angle_error'
TRANSFER_ENERGY = 'opm:transfer_energy'
CURVATURE_RADIUS = 'opm:curvature_radius'
_ORIENTATION = (
    (HYDROPHOBIC_THICKNESS, 'thickness', 'Å'),
    (HYDROPHOBIC_THICKNESS_ERROR, 'thickness_error', 'Å'),
    (TILT_ANGLE, 'tilt', '°'),
    (TILT_ANGLE_ERROR, 'tilt_error', '°'),
    (TRANSFER_ENERGY, 'transfer_energy', 'kcal/mol'),
    (CURVATURE_RADIUS, 'curvature_radius', 'Å'),
)

# The "membranes" of OPM that are not a membrane of a cell or a virion: the
# table also lists the secreted proteins, whose position is calculated in a
# generic bilayer, and the proteins of no known membrane (ids 20 and 23, named
# "Secreted" and "Undefined"). A structure is not located in them; it is
# associated with them.
_NOT_A_MEMBRANE = {'20', '23'}

# The cross-references of the classification, with the prefix of the CURIE.
_XREF_PREFIXES = {'pfam': 'Pfam', 'interpro': 'InterPro', 'tcdb': 'TCDB'}

# UniProt accessions among the codes (``Q9RP17``); the others are entry names.
_ACCESSION = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)
_ENTRY_NAME = re.compile(r'[A-Z0-9]+_[A-Z0-9]+')


def _term_annotations(row: dict[str, Any]) -> list[Annotation]:
    annotations = []
    if row.get('comment'):
        annotations.append(Annotation(slots.description, row['comment']))
    for database, accession in sorted((row.get('xrefs') or {}).items()):
        annotations.append(
            Annotation(slots.xref, f'{_XREF_PREFIXES[database]}:{accession}')
        )
    if row['kind'] == 'membrane':
        for term, key in ((TOPOLOGY_IN, 'inside'), (TOPOLOGY_OUT, 'outside')):
            if row.get(key):
                annotations.append(Annotation(term, row[key]))
    return annotations


def map_term(row: dict[str, Any]) -> Entity:
    """A membrane or a level of the classification, below its parent."""
    kind = row['kind']
    if kind not in TERM_NAMESPACES:
        raise ValueError(f'Unknown OPM term kind: {kind!r}')
    identifiers = [
        Identifier(TERM_NAMESPACES[kind], row['id']),
        Identifier(Namespace.NAME, row['name']),
    ]
    # The short name of a membrane ("Gram-neg. inner") is the label of the
    # tables of structures.
    if row.get('short_name'):
        identifiers.append(Identifier(Namespace.SYNONYM, row['short_name']))
    relations = None
    if row.get('parent_id'):
        relations = [
            OntologyRelation(
                predicate=slots.subclass_of,
                object=EntityRef(
                    OntologyClass,
                    TERM_NAMESPACES[row['parent_kind']],
                    row['parent_id'],
                ),
            )
        ]
    return Entity(
        type=OntologyClass,
        identifiers=identifiers,
        annotations=_term_annotations(row) or None,
        ontology_relations=relations,
    )


def _protein(code: str) -> Entity:
    """The protein of a UniProt code of OPM: an accession or an entry name."""
    if _ACCESSION.fullmatch(code):
        namespace = Namespace.UNIPROT
    elif _ENTRY_NAME.fullmatch(code):
        namespace = Namespace.UNIPROT_ENTRY
    else:
        raise ValueError(f'Unknown OPM UniProt code: {code!r}')
    return Entity(type=Protein, identifiers=[Identifier(namespace, code)])


def structure_entity(row: dict[str, Any]) -> Entity:
    """The protein, or the entry with several proteins, that a structure is of.

    One UniProt code: that protein. Several: the PDB entry, as a complex of the
    proteins, the way the PDB module makes it; OPM does not say which of them
    is in the membrane (the entries of antibody fragments and other partners
    list them too) and the membrane is not asserted of each. None: the PDB
    code, with the name of the structure.
    """
    entries = row['entries']
    if len(entries) == 1:
        return _protein(entries[0])
    identifiers = [Identifier(Namespace.PDB, row['pdb_id'])]
    if row.get('name'):
        identifiers.append(Identifier(Namespace.NAME, row['name']))
    if not entries:
        return Entity(type=Protein, identifiers=identifiers)
    return Entity(
        type=MacromolecularComplex,
        identifiers=identifiers,
        membership=[
            Membership(member=_protein(code), predicate=slots.has_member)
            for code in entries
        ],
    )


def _orientation(row: dict[str, Any]) -> list[Annotation]:
    annotations = []
    for term, key, unit in _ORIENTATION:
        if row.get(key) is None:
            continue
        value = measurement(row[key], unit=unit, source_field=key)
        if value is None:
            raise ValueError(f'OPM structure {row["id"]}: {key} is {row[key]!r}')
        annotations.append(Annotation(term, value))
    return annotations


def map_membrane(row: dict[str, Any]) -> Relation:
    """The structure is in a membrane, in the orientation OPM calculates."""
    predicate = (
        slots.associated_with
        if row['membrane_id'] in _NOT_A_MEMBRANE
        else slots.located_in
    )
    return Relation(
        subject=structure_entity(row),
        predicate=predicate,
        object=EntityRef(OntologyClass, MEMBRANE, row['membrane_id']),
        annotations=[Annotation(PDB_ID, row['pdb_id']), *_orientation(row)],
    )


def map_family(row: dict[str, Any]) -> Relation:
    """The structure belongs to a family; its superfamily, class and type are above it."""
    return Relation(
        subject=structure_entity(row),
        predicate=slots.associated_with,
        object=EntityRef(OntologyClass, FAMILY, row['family_id']),
        annotations=[Annotation(PDB_ID, row['pdb_id'])],
    )


# ================================= RESOURCE ===================================


resource = Resource(
    config,
    terms=Dataset(
        download=download_membranes,
        mapper=map_term,
        raw_parser=partial(
            iter_terms,
            types=download_types,
            classes=download_classes,
            superfamilies=download_superfamilies,
            families=download_families,
        ),
    ),
    membrane_localizations=Dataset(
        download=download_structures,
        mapper=map_membrane,
        raw_parser=iter_structures,
    ),
    family_annotations=Dataset(
        download=download_structures,
        mapper=map_family,
        raw_parser=iter_structures,
    ),
)

# ================================== REFERENCE =================================

# primary_structures
#   id  ordering  family_name_cache  species_name_cache  membrane_name_cache
#   name  description  comments  pdbid  resolution  topology_subunit
#   topology_show_in  thickness  thicknesserror  subunit_segments  tilt
#   tilterror  gibbs  tau  verification  membrane_id  species_id  family_id
#   superfamily_id  classtype_id  type_id  secondary_representations_count
#   structure_subunits_count  citations_count  created_at  updated_at
#   uniprotcode  interpro
#   e.g. 1  6024.0  OmpA family  Escherichia coli  Gram-neg. outer
#        Outer membrane protein A (OMPA), disordered loops  ...  ="1qjp"  1.65
#        A  true  25.4  1.5  8  11  1  -29.5  ""  ...  3  9  1  1  2  ...  OMPA_ECOLI
