"""
Parse the PhosphoSitePlus kinase-substrate table and emit Relation records.

Source file (``Kinase_Substrate_Dataset``, tab separated). Three lines precede
the column names: the release date, the license statement and an empty line.

    050522
    Data extracted from PhosphoSitePlus(R), created by Cell Signaling ...

    GENE  KINASE  KIN_ACC_ID  KIN_ORGANISM  SUBSTRATE  SUB_GENE_ID  SUB_ACC_ID ...
    Pak2  PAK2    Q64303      rat           MEK1       170851       Q01986     ...

Each row is one site-level observation: a kinase, a substrate and the residue
of the substrate the kinase phosphorylates, with the surrounding sequence and
whether the reaction was seen in vivo, in vitro or both. Kinase and substrate
may be of different species; the row is kept as written, with no translation
to another organism.

The publications and methods come from a second file of the same release, the
BioPAX export ``Kinase_substrates.owl`` (see ``parsers/phosphosite.py``), which
lists the same sites. A row takes the PubMed IDs and PSI-MI method terms of the
BioPAX catalysis with the same kinase accession, substrate accession, residue
and position (21,782 of its 21,837 catalyses match a row; 21,782 of the 22,073
rows find one). The BioPAX file attaches its evidence to the phosphorylation
of the substrate site, so every kinase of a site carries the same publications
and methods; the file does not say which of them concern the kinase. It does
not pair a publication with its method either: the relation holds the union of
each. Publication records without a PubMed ID are not carried. Catalyses
without a matching row are dropped, rows without a catalysis have no
publications and no PSI-MI method.

Both files are the copies of the May 2022 release on the OmniPath rescued-data
server (both files were put there within a minute of each other, 16 May 2022). The table carries
the CC BY-NC-SA 3.0 statement PhosphoSitePlus released it under; the BioPAX
file states no license, and the module takes it under the same terms as the
table of the same release. The terms of the PhosphoSitePlus site (July 2026)
no longer include downloads, so later releases need a written agreement with
Cell Signaling Technology.
"""

from __future__ import annotations

import functools
import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.molecular_forms import molecular_form_from_identifiers
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._molecular_forms import combine_forms
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phosphosite import iter_kinase_substrate_evidence
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

_UNIPROT = (
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?'
)
# The accession columns hold a UniProt accession, with the isoform number where
# the record is an isoform, or, for a protein without one, a RefSeq or GenBank
# protein accession.
_ACCESSION_REGEXES = (
    (Namespace.UNIPROT, re.compile(rf'^({_UNIPROT})$')),
    (Namespace.REFSEQ_PROTEIN, re.compile(r'^([NXY]P_\d+(?:\.\d+)?)$')),
    (Namespace.GENBANK, re.compile(r'^([A-Z]{3}\d{5}(?:\.\d+)?)$')),
)
# PhosphoSitePlus marks a record on a UniProt sequence variant
# ("Q12809_VAR_014388") or on a sequence UniProt has withdrawn ("P63058_OBS").
# The accession before the mark is the UniProt entry; the positions are not
# those of its current sequence.
_MARKED_UNIPROT_REGEX = re.compile(rf'^({_UNIPROT})_(VAR_\d+|OBS)$')

# Organism labels that name one species, with its NCBI Taxonomy ID.
ORGANISMS = {
    'human': '9606',
    'mouse': '10090',
    'rat': '10116',
    'rabbit': '9986',
    'cow': '9913',
    'pig': '9823',
    'dog': '9615',
    'sheep': '9940',
    'chicken': '9031',
    'turkey': '9103',
    'guinea pig': '10141',
    'ferret': '9669',
}
# Labels that name a group of species, or none: the participant gets no taxon.
# "hamster" covers golden, Chinese and black-bellied hamster proteins.
ORGANISMS_WITHOUT_TAXON = frozenset(
    {
        'hamster',
        'frog',
        'quail',
        'marmoset',
        'torpedo',
        'starfish',
        'bacteria',
        'herpesvirus',
        'no information',
    }
)
# Every organism label of the file: the reader rejects any other.
ORGANISM_LABELS = frozenset(ORGANISMS) | ORGANISMS_WITHOUT_TAXON
# Evidence & Conclusion Ontology terms of the two reaction columns.
_IN_VIVO = 'ECO:0000178'  # in vivo assay evidence
_IN_VITRO = 'ECO:0000181'  # in vitro assay evidence
# PSI-MI methods of the BioPAX evidence, with their labels in that file: the
# reader rejects any other code, and a label that differs.
METHODS = {
    'MI:0421': 'identification by antibody',
    'MI:0427': 'identification by mass spectrometry',
    'MI:0113': 'western blot',
    'MI:0074': 'mutation analysis',
}


def _accession(value: object) -> tuple[Namespace, str] | None:
    """The namespace and accession of an accession cell; None if it has neither."""
    text = str(value or '').strip()
    for namespace, regex in _ACCESSION_REGEXES:
        if regex.match(text):
            return namespace, text
    marked = _MARKED_UNIPROT_REGEX.match(text)
    return (Namespace.UNIPROT, marked[1]) if marked else None


def _accession_value(column: str, namespace: Namespace):
    def value(row: dict[str, Any]) -> str | None:
        accession = _accession(row.get(column))
        return accession[1] if accession and accession[0] == namespace else None

    return value


def _taxon(column: str):
    def value(row: dict[str, Any]) -> str | None:
        # The reader has rejected labels outside the two sets above.
        taxon = ORGANISMS.get(row.get(column))
        return f'NCBITaxon:{taxon}' if taxon else None

    return value


def _sequence_form(column: str):
    """The isoform or RefSeq sequence the accession cell names, if it names one."""

    def form(row: dict[str, Any]) -> dict | None:
        text = str(row.get(column) or '').strip()
        accession = _accession(text)
        # A marked accession is not the entry's sequence or one of its isoforms.
        if not accession or accession[1] != text:
            return None
        return molecular_form_from_identifiers(
            [{'ns': str(accession[0]), 'id': accession[1]}]
        )

    return form


_substrate_sequence_form = _sequence_form('SUB_ACC_ID')


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The phosphorylated residue, a feature of the substrate in this record.

    The position is the source's, in the sequence of the substrate's accession;
    it is not checked against that sequence here. A record on a variant or a
    withdrawn sequence, or without a recognised accession, has no reference
    sequence. The description is the site and its sequence window as written.
    """
    if row.get('position') is None:
        return None
    text = str(row.get('SUB_ACC_ID') or '').strip()
    accession = _accession(text)
    reference = (
        {'ns': str(accession[0]), 'id': accession[1]}
        if accession and accession[1] == text
        else None
    )
    marked = _MARKED_UNIPROT_REGEX.match(text)
    site = {
        'modifications': [
            {
                'term': 'phosphorylated residue',
                'residue': row['residue'],
                'position': row['position'],
                'end_position': row['position'],
                'coordinate_reference': {
                    'identifier': reference,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                'description': f'{row["SUB_MOD_RSD"]} ({row["SITE_+/-7_AA"]})',
            }
        ],
        # The UniProt variant the record is on, by its feature identifier.
        'variants': [{'description': marked[2]}]
        if marked and marked[2].startswith('VAR_')
        else None,
    }
    return combine_forms(_substrate_sequence_form(row), site)


def _evidence_types(row: dict[str, Any]) -> list[str]:
    """The in vivo and in vitro terms of the table, then the PSI-MI methods."""
    return [
        term
        for key, term in (('in_vivo', _IN_VIVO), ('in_vitro', _IN_VITRO))
        if row.get(key)
    ] + list(row.get('evidence_methods') or ())


def _publications(row: dict[str, Any]) -> list[str]:
    # The reader has kept the well-formed PubMed IDs.
    return [f'PMID:{pmid}' for pmid in row.get('publications') or ()]


config = ResourceConfig(
    id=ResourceCv.PHOSPHOSITE,
    name='PhosphoSite',
    url='https://www.phosphosite.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_NC_SA_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22135298',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Kinase-substrate phosphorylation sites with their sequence context, '
        'in vivo or in vitro evidence, and the publications and PSI-MI methods '
        'of the BioPAX file, May 2022 release.'
    ),
    description=(
        'PhosphoSitePlus is a resource of Cell Signaling Technology for '
        'experimentally observed post-translational modifications, mainly '
        'of human, mouse and rat proteins, curated from the literature and '
        'from mass spectrometry data. Its kinase-substrate table lists the '
        'kinases reported to phosphorylate each site, the site with its '
        'surrounding sequence, and whether the reaction was observed in '
        'vivo or in vitro. A BioPAX export of the same release gives the '
        'publications and detection methods recorded for each site; they '
        'belong to the phosphorylation of the site, not to one kinase.'
    ),
)

f = FieldConfig()


def _protein(
    accession: str,
    gene_symbol: str,
    name: str,
    organism: str,
    gene_id: str | None = None,
    molecular_form=None,
):
    builder = EntityBuilder(
        entity_type=model.Protein,
        molecular_form=molecular_form,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=_accession_value(accession, Namespace.UNIPROT),
            ),
            CV(
                term=Namespace.REFSEQ_PROTEIN,
                value=_accession_value(accession, Namespace.REFSEQ_PROTEIN),
            ),
            CV(
                term=Namespace.GENBANK,
                value=_accession_value(accession, Namespace.GENBANK),
            ),
            CV(term=Namespace.GENESYMBOL, value=f(gene_symbol)),
            *([CV(term=Namespace.ENTREZ, value=f(gene_id))] if gene_id else []),
            # The PhosphoSitePlus name of the protein, e.g. "Tau iso8".
            CV(term=Namespace.NAME, value=f(name)),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=_taxon(organism)),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        # An accession or a gene symbol identifies the participant; the
        # PhosphoSitePlus name alone does not.
        identified = _accession(row.get(accession)) or (row.get(gene_symbol) or '').strip()
        return builder.build(row) if identified else None

    return build


# A row whose kinase has neither an accession nor a gene symbol yields no
# relation.
enzyme_substrate_schema = RelationBuilder(
    subject=_protein(
        'KIN_ACC_ID',
        'GENE',
        'KINASE',
        'KIN_ORGANISM',
        molecular_form=_sequence_form('KIN_ACC_ID'),
    ),
    predicate=slots.affects,
    # The phosphorylated residue is in the substrate's molecular form.
    object=_protein(
        'SUB_ACC_ID',
        'SUB_GENE',
        'SUBSTRATE',
        'SUB_ORGANISM',
        gene_id='SUB_GENE_ID',
        molecular_form=_substrate_form,
    ),
    annotations=AnnotationsBuilder(
        # The kinase increases the substrate's phosphorylation; the table says
        # nothing about the substrate's activity.
        CV(
            term=slots.object_direction_qualifier,
            value=model.DirectionQualifierEnum.increased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=model.GeneOrGeneProductOrChemicalEntityAspectEnum.phosphorylation,
        ),
        CV(term=slots.has_evidence_of_type, value=_evidence_types),
        CV(term=slots.publications, value=_publications),
    ),
)

# The BioPAX export of the same release, read by the parser of the table.
evidence_download = Download(
    url=(
        'https://rescued.omnipathdb.org/phosphosite/'
        'Kinase_substrates.owl.gz'
    ),
    filename='Kinase_substrates.owl.gz',
    subfolder='phosphosite',
    large=True,
    ext='gz',
)

resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url=(
                'https://rescued.omnipathdb.org/phosphosite/'
                'Kinase_Substrate_Dataset.gz'
            ),
            filename='Kinase_Substrate_Dataset.gz',
            subfolder='phosphosite',
            large=True,
            # One byte of the license statement is not UTF-8.
            encoding='iso-8859-1',
            ext='gz',
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=functools.partial(
            iter_kinase_substrate_evidence,
            organisms=ORGANISM_LABELS,
            methods=METHODS,
            evidence_download=evidence_download,
        ),
    ),
)
