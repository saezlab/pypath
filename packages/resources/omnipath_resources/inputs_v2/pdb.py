"""
Parse the PDB structures that contain several UniProt proteins and emit Entity records.

The Protein Data Bank (wwPDB) archives the experimentally determined
structures of biological macromolecules. This module makes a complex of the
UniProt proteins that the chains of one PDB entry map to, as the SIFTS table
``pdb_chain_uniprot`` (PDBe and UniProt) gives them. The complex is the
PDB entry: its accession is the PDB identifier and its members are the
distinct UniProt accessions of its chains.

What the source states and what the grouping infers. SIFTS states which
UniProt entry each chain of a structure maps to. That the proteins of one
entry form a complex is inferred from the entry: they were crystallised or
imaged together, which is not always a physiological complex. An entry can
hold a fusion construct, an antibody or nanobody used to crystallise the
target, a peptide, a bound cofactor protein, or proteins that are in the
crystal by accident of the preparation. Nothing in the table tells which. The
stoichiometry is not asserted either: the chains of the entry are those of
the deposited model (often the asymmetric unit, not the biological
assembly), and a protein with several chains is a member once.

Entries with a single UniProt protein (homo-oligomers and monomers) and
entries whose chains all carry the same proteins (one chain made of segments
of several proteins, or copies of the same fusion) are not complexes of
several proteins and are not emitted. The structure method and
resolution are not in the table, and are not represented.

The legacy input read the same table and made a complex of the chains of
each entry with more than one chain, lower-casing the PDB identifier; it kept
homo-oligomers (the same protein once per chain) and merged complexes with
the same components.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import MacromolecularComplex, Protein
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.pdb import iter_pdb_entries
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Entity, Identifier, Membership

config = ResourceConfig(
    id=ResourceCv.PDB,
    name='PDB',
    full='Protein Data Bank',
    url='https://www.wwpdb.org/',
    # The files of the PDB archive are CC0 1.0 (wwPDB usage policies,
    # https://www.wwpdb.org/about/usage-policies). The SIFTS mapping states no
    # license of its own; EMBL-EBI's terms of use make its data freely available.
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='14634627',
    content=(Content.PROTEIN_COMPLEXES,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL,),
    summary=(
        'Sets of UniProt proteins found together in one experimental PDB structure '
        '(derived groupings, not necessarily physiological complexes).'
    ),
    description=(
        'The Protein Data Bank is the archive of experimentally determined '
        'three-dimensional structures of proteins, nucleic acids and their complexes, '
        'run by the wwPDB (RCSB PDB, PDBe, PDBj, BMRB). This resource groups the UniProt '
        'proteins that the chains of one PDB entry map to, with the SIFTS mapping of '
        'PDBe and UniProt. The groups say which proteins were crystallised or imaged '
        'together; they include fusion partners, antibodies and peptides, and are not '
        'all physiological complexes.'
    ),
)


def map_complex(row: dict[str, Any]) -> Entity | None:
    """The complex of a PDB entry with its UniProt proteins; None below two proteins."""
    accessions = row.get('accessions') or []
    if len(accessions) < 2:
        return None
    return Entity(
        type=MacromolecularComplex,
        identifiers=[Identifier(type=Namespace.PDB, value=row['pdb_id'])],
        membership=[
            Membership(
                member=Entity(
                    type=Protein,
                    identifiers=[Identifier(type=Namespace.UNIPROT, value=accession)],
                )
            )
            for accession in accessions
        ],
    )


# The current release of the table; the release is in its first line.
resource = Resource(
    config,
    complexes=Dataset(
        download=Download(
            url=(
                'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/'
                'pdb_chain_uniprot.tsv.gz'
            ),
            filename='pdb_chain_uniprot.tsv.gz',
            subfolder='pdb',
            large=True,
            ext='gz',
        ),
        mapper=map_complex,
        raw_parser=iter_pdb_entries,
    ),
)
