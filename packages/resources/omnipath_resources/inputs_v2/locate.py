"""
Parse LOCATE and emit the subcellular localizations of human and mouse proteins.

LOCATE (Sprenger et al. 2008) is a database of the subcellular localization of
mammalian proteins, with the membrane topology class of each protein. A protein
has localizations of three kinds, which the file keeps apart and the statements
carry as the attribute ``locate:evidence_type``:

- ``literature``: manually curated from a publication, with its PubMed
  identifier, the organism the study used and the curator's notes;
- ``external annotation``: taken from another database (HPRD, UniProt, LIFEdb or
  Ensembl), with the database, the protein's record in it and, for a few, a GO
  evidence code;
- ``prediction``: a prediction by CELLO, MultiLoc, Proteome Analyst, WoLFPSORT or
  pTarget, with the method's label for the location and its score.

The experimental images of LOCATE (``experimental_data``) are empty in this
release. Each statement is ``located_in`` a GO term, which the file gives. The
protein's class (``cytoplasmic``, ``secretome``, ``typeI``, ``typeII`` or
``mtmp``) is an ``associated_with`` statement to a term of LOCATE.

The legacy reader returned the first two kinds for human proteins whose
accession it could translate to a Swiss-Prot entry, and the predictions only on
request, and kept the location names, not the GO terms. Here the protein is
identified by the accessions of the file, and the GO term is the term.
The HGNC id (human) and the MGI id (mouse) of the record are secondary
identifiers of the protein, for the many records without a UniProt accession.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.locate import (
    EXTERNAL,
    LITERATURE,
    PREDICTION,
    PROTEIN_CLASS,
    iter_localizations,
    iter_protein_classes,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation

config = ResourceConfig(
    id=ResourceCv.LOCATE,
    name='LOCATE',
    full='LOCATE subcellular localization database',
    url='http://locate.imb.uq.edu.au/',
    # The site is gone and the registry records the journal's license for the paper
    # (Nucleic Acids Res 2008, open access).
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.NAR_OPEN_ACCESS,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='17986452',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Subcellular localizations of human and mouse proteins as Gene Ontology '
        'terms, from the literature, other databases and predictions, with the '
        'protein topology class.'
    ),
    description=(
        'LOCATE is a curated database of the subcellular localization of '
        'mammalian proteins. Its entries combine localizations curated from the '
        'literature, annotations from other databases, predictions of the '
        'localization from the sequence, the membrane organization of each '
        'protein (cytoplasmic, secreted, type I or II, multi-pass membrane) '
        'and, in principle, images of experiments.'
    ),
)

TAXA = {'Human': 'NCBITaxon:9606', 'Mouse': 'NCBITaxon:10090'}
CLASS_NAMESPACE = 'locate_protein_class'
CLASSES = frozenset({'cytoplasmic', 'mtmp', 'secretome', 'typeI', 'typeII'})
# Databases of the external annotations and methods of the predictions, as the
# file writes them. A name that is not here raises.
EXTERNAL_DATABASES = frozenset(
    {
        'Ensembl-Peptide Human', 'Ensembl-Peptide Mouse', 'HPRD', 'LIFEdb', 'MGI ID',
        'UniProt/SPTrEMBL',
    }
)
METHODS = frozenset({'CELLO', 'MultiLoc', 'Proteome_Analyst', 'WoLFPSORT', 'pTarget'})

# The GO evidence codes of the external annotations (the root codes of the GAF
# format) and the ECO terms the Evidence & Conclusion Ontology assigns them by
# default: http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt (the version
# of January 2026; the file lists these 26 codes). A code that is not here raises.
ECO_OF_GO_EVIDENCE_CODE = {
    'EXP': 'ECO:0000269', 'HDA': 'ECO:0007005', 'HEP': 'ECO:0007007', 'HGI': 'ECO:0007003',
    'HMP': 'ECO:0007001', 'HTP': 'ECO:0006056', 'IBA': 'ECO:0000318', 'IBD': 'ECO:0000319',
    'IC': 'ECO:0000305', 'IDA': 'ECO:0000314', 'IEA': 'ECO:0000501', 'IEP': 'ECO:0000270',
    'IGC': 'ECO:0000317', 'IGI': 'ECO:0000316', 'IKR': 'ECO:0000320', 'IMP': 'ECO:0000315',
    'IPI': 'ECO:0000353', 'IRD': 'ECO:0000321', 'ISA': 'ECO:0000247', 'ISM': 'ECO:0000255',
    'ISO': 'ECO:0000266', 'ISS': 'ECO:0000250', 'NAS': 'ECO:0000303', 'ND': 'ECO:0000307',
    'RCA': 'ECO:0000245', 'TAS': 'ECO:0000304',
}

# Source attributes of a statement. The prediction score has the scale of its
# method (in the human file: CELLO 0.72 to 5, MultiLoc and Proteome_Analyst 0 to 1,
# pTarget 6.3 to 100, WoLFPSORT 1 to 32), so it is not a confidence score; it is a
# Measurement whose source field names the method.
_EVIDENCE_TYPE = 'locate:evidence_type'
_STUDY_ORGANISM = 'locate:study_organism'
_RECORD = 'locate:source_record'
_EVIDENCE_CODE = 'locate:evidence_code'
_METHOD = 'locate:prediction_method'
_PREDICTED_LOCATION = 'locate:predicted_location'
_SCORE = 'locate:prediction_score'

# The identifier namespaces of the keys of the row's `ids`, in the order the
# protein lists them.
_NAMESPACES = (
    ('uniprot', Namespace.UNIPROT),
    ('uniprot_entry', Namespace.UNIPROT_ENTRY),
    ('entrez', Namespace.ENTREZ),
    ('ensp', Namespace.ENSP),
    ('refseq_protein', Namespace.REFSEQ_PROTEIN),
    ('genbank_protein', Namespace.GENBANK),
)
# The HGNC and MGI ids are gene-level identifiers on a protein record, as the
# Entrez Gene id above is: LOCATE gives them for the records that have no UniProt
# accession as well (40% of the human and 61% of the mouse proteins have none).
# They come after all the other identifiers, LOCATE's own included, so the
# leading identifier of a protein does not change.
_GENE_NAMESPACES = (
    ('hgnc', Namespace.HGNC),
    ('mgi', Namespace.MGI),
)


def _protein(row: dict[str, Any]) -> Entity:
    if row['organism'] not in TAXA:
        raise ValueError(f'Unknown LOCATE organism: {row["organism"]!r}')
    identifiers = [
        Identifier(type=namespace, value=value)
        for key, namespace in _NAMESPACES
        for value in row['ids'][key]
    ]
    # LOCATE's own identifier, for the protein records that have no other.
    identifiers.append(Identifier(type='locate', value=row['ids']['locate']))
    identifiers.extend(
        Identifier(type=namespace, value=value)
        for key, namespace in _GENE_NAMESPACES
        for value in row['ids'][key]
    )
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXA[row['organism']])],
    )


def _checked(value: str | None, allowed: frozenset[str], what: str) -> str:
    if value not in allowed:
        raise ValueError(f'Unknown LOCATE {what}: {value!r}')
    return value


def _prediction_score(row: dict[str, Any]):
    """The score of a prediction, with the method whose scale it is on."""
    score = measurement(row['score'], source_field=f'LOCATE score of {row["method"]}')
    if score is None:
        raise ValueError(f'Not a prediction score: {row["score"]!r}')
    return score


def map_localization(row: dict[str, Any]) -> Relation:
    """A protein located in a GO term, or associated with its LOCATE class."""
    if row['kind'] == PROTEIN_CLASS:
        return Relation(
            subject=_protein(row),
            predicate=slots.associated_with,
            object=Entity(
                type=OntologyClass,
                identifiers=[
                    Identifier(type=CLASS_NAMESPACE, value=_checked(row['class'], CLASSES, 'class'))
                ],
            ),
        )
    annotations = [Annotation(term=_EVIDENCE_TYPE, value=row['kind'])]
    if row['kind'] == LITERATURE:
        if row['pmid']:
            curie = pmid_curie(row['pmid'])
            if curie is None:
                raise ValueError(f'Not a PubMed identifier: {row["pmid"]!r}')
            annotations.append(Annotation(term=slots.publications, value=curie))
        if row['study_organism']:
            annotations.append(Annotation(term=_STUDY_ORGANISM, value=row['study_organism']))
        if row['notes']:
            annotations.append(Annotation(term=slots.description, value=row['notes']))
    elif row['kind'] == EXTERNAL:
        database = _checked(row['database'], EXTERNAL_DATABASES, 'database')
        annotations.append(Annotation(term=slots.supporting_data_source, value=database))
        if row['record']:
            annotations.append(Annotation(term=_RECORD, value=f'{database}:{row["record"]}'))
        if row['evidence_code']:
            annotations.append(Annotation(term=_EVIDENCE_CODE, value=row['evidence_code']))
            annotations.append(
                Annotation(
                    term=slots.has_evidence_of_type,
                    value=ECO_OF_GO_EVIDENCE_CODE[
                        _checked(row['evidence_code'], frozenset(ECO_OF_GO_EVIDENCE_CODE), 'evidence code')
                    ],
                )
            )
    elif row['kind'] == PREDICTION:
        annotations.extend(
            [
                Annotation(term=_METHOD, value=_checked(row['method'], METHODS, 'method')),
                Annotation(term=_PREDICTED_LOCATION, value=row['location']),
                Annotation(term=_SCORE, value=_prediction_score(row)),
            ]
        )
    else:
        raise ValueError(f'Unknown kind of LOCATE statement: {row["kind"]!r}')
    return Relation(
        subject=_protein(row),
        predicate=slots.located_in,
        object=Entity(
            type=OntologyClass,
            identifiers=[Identifier(type=Namespace.GO, value=row['go'])],
        ),
        annotations=annotations,
    )


def map_class(row: dict[str, Any]) -> Entity:
    """A topology class of LOCATE."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=CLASS_NAMESPACE, value=_checked(row['name'], CLASSES, 'class')),
            Identifier(type=Namespace.NAME, value=row['name']),
        ],
    )


def _download(species: str) -> Download:
    """The XML file of a species. The human file is the OmniPath rescue copy of
    the file of the site, which is gone; the mouse file is the Internet Archive's
    copy of the site's file (captured 2019-03-06), which the rescue server lacks.
    """
    name = f'LOCATE_{species}_v6_20081121.xml'
    if species == 'human':
        url = f'https://rescued.omnipathdb.org/{name}.zip'
    else:
        url = (
            'https://web.archive.org/web/20190306074407id_/'
            f'http://locate.imb.uq.edu.au/info_files/{name}.zip'
        )
    return Download(
        url=url,
        filename=f'{name}.zip',
        subfolder='locate',
        large=True,
        default_mode='rb',
        encoding=None,
        ext='zip',
        needed=[name],
    )


human_download = _download('human')
mouse_download = _download('mouse')

resource = Resource(
    config,
    terms=Dataset(
        download=human_download,
        mapper=map_class,
        raw_parser=iter_protein_classes,
    ),
    human=Dataset(
        download=human_download,
        mapper=map_localization,
        raw_parser=iter_localizations,
    ),
    mouse=Dataset(
        download=mouse_download,
        mapper=map_localization,
        raw_parser=iter_localizations,
    ),
)
