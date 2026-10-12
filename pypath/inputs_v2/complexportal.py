"""
Map the Complex Portal ComplexTAB files of human complexes to complexes.

The Complex Portal (EMBL-EBI) describes macromolecular complexes of a few
species. Every complex has a ``CPX`` accession, a recommended name, aliases, the
molecules it consists of with their stoichiometry, an evidence code, GO
annotations and cross references. The molecules can be proteins (UniProt
accessions, with isoforms and chains), small molecules (ChEBI), RNAs
(RNAcentral) and other complexes; a defined set of molecules that cannot be
told apart (paralogous genes) is one member. All of them are members of the
complex as they are written: a sub-complex stays a member, nothing is
flattened to proteins, and no pairs are made. A stoichiometry of 0 means that
it is not known and is not emitted.

Two files of human complexes are mapped, one dataset each:

- ``human_complexes``: the manually curated complexes (``9606.tsv``);
- ``human_predicted_complexes``: the complexes predicted by machine learning
  from hu.MAP 3.0 and MuSIC (``9606_predicted.tsv``).

The other species of the Complex Portal have a file of the same layout
(``{taxon}.tsv``, for example 10090 for mouse); a dataset for one is a
``_complextab`` call with its file name.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import (
    ChemicalEntity,
    MacromolecularComplex,
    NucleicAcidEntity,
    OntologyClass,
    Protein,
    ProteinFamily,
    RNAProduct,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.complexportal import iter_complextab
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.COMPLEXPORTAL,
    name='Complex Portal',
    url='https://www.ebi.ac.uk/complexportal/',
    # The Complex Portal site states CC0 for its data (checked when the module
    # was written, October 2026); the legacy registry recorded CC BY 4.0.
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='39558156',
    content=(Content.PROTEIN_COMPLEXES,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Human macromolecular complexes with their members, stoichiometry, '
        'evidence code and GO annotations.'
    ),
    description=(
        'The Complex Portal is a manually curated, encyclopaedic resource of '
        'stable macromolecular complexes of known function from key model '
        'organisms. Complexes can contain proteins, nucleic acids and small '
        'molecules, and other complexes. For the human complexome it also '
        'holds high confidence complexes predicted by machine learning from '
        'hu.MAP 3.0 and MuSIC.'
    ),
)

_FTP = 'https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/'
_UNIPROT = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+|-PRO_\d+)?'
)
_COMPLEX = re.compile(r'CPX-\d+')
_CHEBI = re.compile(r'CHEBI:(\d+)')
_RNACENTRAL = re.compile(r'URS[0-9A-F]{10}_(\d+)')
_INSDC = re.compile(r'[A-Z]\d{5}')
_EC = re.compile(r'\d+(?:\.(?:n?\d+|-))*')

# The curating database named in the "Source" column (observed values).
_SOURCES = {
    name: name
    for name in (
        'IntAct', 'UniProt', 'ceitec', 'bhf-ucl', 'NTNU', 'PDBe', 'matrixdb',
        'HuMap', 'MuSIC',
    )
}

# Cross references with the qualifier "identity" (the same object in another
# database): source database -> accession pattern, CURIE prefix.
_IDENTICAL = {
    'reactome': (re.compile(r'R-[A-Z]{3}-\d+(?:\.\d+)?'), 'REACTOME'),
    'signor': (re.compile(r'SIGNOR-C\d+'), 'SIGNOR'),
    'humap': (re.compile(r'huMAP\d_\d+(?:\.\d+)?'), 'HUMAP'),
    'music': (re.compile(r'[A-Z]?\d+'), 'MUSIC'),
    'qsproteome': (re.compile(r'QS\d+'), 'QSPROTEOME'),
    'wwpdb': (re.compile(r'[0-9][A-Za-z0-9]{3}'), 'PDB'),
    'emdb': (re.compile(r'EMD-\d+'), 'EMDB'),
    'chembl target': (re.compile(r'CHEMBL\d+'), 'CHEMBL.TARGET'),
    'matrixdb': (re.compile(r'MULT_.+'), 'MATRIXDB'),
    'protein ontology': (re.compile(r'PR[:_](\d+)'), 'PR'),
    'intact': (re.compile(r'EBI-\d+'), 'INTACT'),
}


def _member_entity(accession: str) -> Entity:
    """The member that a molecule accession names; its kind is read from the accession."""
    if _COMPLEX.fullmatch(accession):
        return Entity(
            type=MacromolecularComplex,
            identifiers=[Identifier(Namespace.COMPLEXPORTAL, accession)],
        )
    if match := _CHEBI.fullmatch(accession):
        return Entity(
            type=ChemicalEntity,
            identifiers=[Identifier(Namespace.CHEBI, match[1])],
        )
    if match := _RNACENTRAL.fullmatch(accession):
        # An RNAcentral ID ends with the taxon of the sequence.
        return Entity(
            type=RNAProduct,
            identifiers=[Identifier(Namespace.RNACENTRAL, accession)],
            annotations=[Annotation(slots.in_taxon, f'NCBITaxon:{match[1]}')],
        )
    if _UNIPROT.fullmatch(accession):
        return Entity(
            type=Protein,
            identifiers=[Identifier(Namespace.UNIPROT, accession)],
        )
    if _INSDC.fullmatch(accession):
        # A nucleotide sequence record (ENA/GenBank); the file does not say
        # more than that it is a nucleic acid, which the accession shows.
        return Entity(
            type=NucleicAcidEntity,
            identifiers=[Identifier(Namespace.GENBANK, accession)],
        )
    raise ValueError(f'Unknown kind of Complex Portal molecule: {accession!r}')


def _molecule(member: dict[str, Any]) -> Entity:
    """A molecule, or a defined set of molecules as a family of its members."""
    if member['set'] is None:
        return _member_entity(member['accession'])
    members = [_member_entity(accession) for accession in member['set']]
    if any(entity.type is not Protein for entity in members):
        raise ValueError(f'Complex Portal set of other than proteins: {member["set"]!r}')
    return Entity(
        type=ProteinFamily,
        # The file writes a set as the list of its members; that is its name.
        identifiers=[Identifier(Namespace.NAME, '[' + ','.join(member['set']) + ']')],
        membership=[Membership(member=entity) for entity in members],
    )


def _membership(member: dict[str, Any]) -> Membership:
    return Membership(
        member=_molecule(member),
        # 0 is "no information about stoichiometry" in the Complex Portal.
        annotations=(
            [Annotation(slots.stoichiometry, str(member['stoichiometry']))]
            if member['stoichiometry']
            else None
        ),
    )


def _xref(xref: dict[str, str]) -> Annotation | None:
    """An identical object in another database."""
    if xref['qualifier'] != 'identity' or xref['database'] not in _IDENTICAL:
        return None
    pattern, prefix = _IDENTICAL[xref['database']]
    match = pattern.fullmatch(xref['accession'])
    if not match:
        return None
    accession = xref['accession']
    if prefix == 'PDB':
        accession = accession.upper()
    elif prefix == 'PR':
        accession = match[1]
    return Annotation(slots.xref, f'{prefix}:{accession}')


def _publication(xref: dict[str, str]) -> Annotation | None:
    """Additional literature of the complex."""
    if (
        xref['database'] == 'pubmed'
        and xref['qualifier'] == 'see-also'
        and (curie := pmid_curie(xref['accession']))
    ):
        return Annotation(slots.publications, curie)
    return None


def _enzyme_class(xref: dict[str, str]) -> Annotation | None:
    """The enzymatic activity of the complex, an EC number."""
    if (
        xref['database'] == 'intenz'
        and xref['qualifier'] == 'identity'
        and _EC.fullmatch(xref['accession'])
    ):
        return Annotation(slots.has_topic, f'EC:{xref["accession"]}')
    return None


def _identifiers(row: dict[str, Any]) -> list[Identifier]:
    identifiers = [Identifier(Namespace.COMPLEXPORTAL, row['accession'])]
    for xref in row['xrefs']:
        if xref['qualifier'] != 'secondary-ac':
            continue
        # The accession of a complex that was merged into this one, and the
        # IntAct accession of the complex.
        if xref['database'] == 'complex portal':
            identifiers.append(Identifier(Namespace.COMPLEXPORTAL, xref['accession']))
        elif xref['database'] == 'intact':
            identifiers.append(Identifier(Namespace.INTACT, xref['accession']))
    if row['name']:
        identifiers.append(Identifier(Namespace.NAME, row['name']))
    identifiers.extend(
        Identifier(Namespace.SYNONYM, synonym)
        for synonym in dict.fromkeys(row['synonyms'])
        if synonym != row['name']
    )
    return identifiers


def map_complex(row: dict[str, Any]) -> Entity:
    """A complex with its molecules, evidence code, GO annotations and cross references."""
    source = row['source']['name']
    if source not in _SOURCES:
        raise ValueError(f'Unknown Complex Portal source: {source!r}')
    xrefs = row['xrefs']
    annotations = [
        Annotation(slots.in_taxon, f'NCBITaxon:{row["taxon"]}'),
        Annotation(slots.has_evidence_of_type, row['evidence_code']),
        Annotation(slots.supporting_data_source, _SOURCES[source]),
        *([Annotation(slots.description, row['description'])] if row['description'] else []),
        *filter(None, map(_publication, xrefs)),
        *filter(None, map(_enzyme_class, xrefs)),
        *filter(None, map(_xref, xrefs)),
    ]
    return Entity(
        type=MacromolecularComplex,
        identifiers=_identifiers(row),
        annotations=annotations,
        membership=[_membership(member) for member in row['members']],
        associations=[
            Association(
                object=EntityRef(OntologyClass, Namespace.GO, term),
                predicate=slots.associated_with,
            )
            for term in row['go_terms']
        ],
    )


def _complextab(filename: str) -> Dataset:
    """A dataset of one ComplexTAB file; the species is the one of each row."""
    return Dataset(
        download=Download(
            url=_FTP + filename,
            filename=f'complexportal_{filename}',
            subfolder='complexportal',
            ext='tsv',
        ),
        mapper=map_complex,
        raw_parser=iter_complextab,
    )


resource = Resource(
    config,
    human_complexes=_complextab('9606.tsv'),
    human_predicted_complexes=_complextab('9606_predicted.tsv'),
)
