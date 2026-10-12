"""
Parse the TRIP Database export and emit TRP channel interaction records.

TRIP (mammalian TRansient receptor potential channel-Interacting Protein
database) offered its curated table as one CSV file. Each line is one
observation from one publication: a TRP channel, an interacting protein, the
method, the details of one evidence category (screening, in vitro validation,
in vivo validation, characterization or functional consequence) and a PubMed
ID. The file names proteins the way TRIP does ("TRPP1", "PKD1/Polycystin-1");
the site's protein list gives the UniProt accession of each name in human,
mouse and rat.

trpchannel.org is gone, so both files are read from the Internet Archive: the
last export (6 August 2015) and the last capture of the protein list (2021).
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import (
    DirectionQualifierEnum,
    GeneOrGeneProductOrChemicalEntityAspectEnum,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.trip import iter_trip_observations
from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)

_ASPECT = GeneOrGeneProductOrChemicalEntityAspectEnum
_DIRECTION = DirectionQualifierEnum

# Method labels -> PSI-MI interaction detection method, where PSI-MI has the
# method under this name or as a synonym. A label mapped to None is known but
# has no such term (functional and modification assays, mostly).
METHOD_TERMS = {
    'Co-immunoprecipitation': 'MI:0019',  # coimmunoprecipitation
    'Yeast two-hybrid': 'MI:0018',  # two hybrid
    'Split-Ubiquitin Yeast two-hybrid': 'MI:0112',  # ubiquitin reconstruction
    'Mammalian two-hybrid': 'MI:0728',  # gal4 vp16 complementation
    'Fusion protein-pull down assay': 'MI:0096',  # pull down
    'Affinity purification-mass spectrometry': 'MI:0004',  # affinity chromatography technology
    'Fluorescence resonance energy transfer': 'MI:0055',
    'Biomolecular fluorescence complementation': 'MI:0809',
    'Co-immunofluorescence staining': 'MI:0403',  # colocalization
    'X-ray diffraction-structural analysis': 'MI:0114',  # x-ray crystallography
    'Atomic force microscopy(AFM) analysis': 'MI:0872',
    'Fast performance liquid chromatography': None,
    'Native-gel electrophoresis': None,
    'Single-molecule photobleaching': None,
    'Fluorescence probe labeling': None,
    'Cell surface biotinylation': None,
    'Patch clamp': None,
    'Calcium measurement': None,
    'Ion influx measurement': None,
    'In vitro PTM assay': None,
    'In vivo PTM assay': None,
    # Non-experimental screening.
    'Inference': None,
    'Speculation': None,
}

# Functional consequence of the interaction on the TRP channel -> (aspect,
# direction). A direction is given only where the label states one: the level
# changes of "Subcellular trafficking" are relative to a compartment that the
# qualifiers cannot name, and "Glycosylation" does not say gain or loss.
ACTIVITY_EFFECTS = {
    'Activation': (_ASPECT.activity, _DIRECTION.increased),
    'Inhibition': (_ASPECT.activity, _DIRECTION.decreased),
    'Sensitization': (_ASPECT.activity, None),
    'Desensitization': (_ASPECT.activity, None),
    # Two channel subunits form a channel with new properties.
    'New channel creation': (None, None),
}
TRAFFICKING_EFFECTS = {
    'Increase in plasma membrane level': (_ASPECT.localization, None),
    'Decrease in plasma membrane level': (_ASPECT.localization, None),
    'Increase in lysosomal membrane level': (_ASPECT.localization, None),
    'Retain in the endoplasmic reticulum': (_ASPECT.localization, None),
    'Internalization from membrane by ligand': (_ASPECT.localization, None),
}
PTM_EFFECTS = {
    'Phosphorylation': (_ASPECT.phosphorylation, _DIRECTION.increased),
    'Phosphorylayion': (_ASPECT.phosphorylation, _DIRECTION.increased),  # sic
    'Dephosphorylation': (_ASPECT.phosphorylation, _DIRECTION.decreased),
    'Ubiquitination': (_ASPECT.ubiquitination, _DIRECTION.increased),
    'De-ubiquitination': (_ASPECT.ubiquitination, _DIRECTION.decreased),
    'Glycosylation': (_ASPECT.glycosylation, None),
}

# Modification type -> PSI-MOD term of the modified residue, where PSI-MOD
# names the residue modification under this name (as in ``hprd.py``). The
# others have no such term and the entry keeps the source's type as its
# description. A dephosphorylation acts on a phosphorylated residue.
PTM_RESIDUE_TERMS = {
    'Phosphorylation': 'phosphorylated residue',
    'Phosphorylayion': 'phosphorylated residue',
    'Dephosphorylation': 'phosphorylated residue',
    'Ubiquitination': None,
    'De-ubiquitination': None,
    'Glycosylation': None,
}
# Three-letter amino acid names of the residue in the modification column.
AMINO_ACIDS = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Gln': 'Q',
    'Glu': 'E',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
}

EXCLUDED_NON_EXPERIMENTAL = 'non-experimental screening (prediction)'
EXCLUDED_MIRRORED = 'second listing of a TRP-TRP observation'

config = ResourceConfig(
    id=ResourceCv.TRIP,
    name='TRIP',
    full='TRIP Database',
    url='http://www.trpchannel.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_NC_SA_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='23071747',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS, Content.CAUSAL_REGULATION),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Literature-curated interactions of mammalian TRP channels with other '
        'proteins, one record per method and publication.'
    ),
    description=(
        'The TRIP Database is a manually curated collection of '
        'protein-protein interactions of the 28 mammalian transient receptor '
        'potential (TRP) channels. For each interaction reported in a '
        'publication it records how the interaction was found (screening), '
        'validated in vitro and in vivo and characterized, and its functional '
        'consequence on the channel. It was maintained at Seoul National '
        'University and Myongji University; the website is no longer online.'
    ),
)


def trip_exclusion(row: dict[str, Any]) -> str | None:
    """Why an observation is not mapped, or None if it is."""
    flag = row.get('screening_non_experimental')
    if flag == 'Prediction':
        # Interactions the authors inferred or speculated on, without an
        # experiment (methods "Inference" and "Speculation").
        return EXCLUDED_NON_EXPERIMENTAL
    if flag:
        raise ValueError(f'Unknown TRIP non-experimental screening: {flag!r}')
    if row.get('mirrored'):
        return EXCLUDED_MIRRORED
    return None


def functional_effect(row: dict[str, Any]) -> tuple[Any, Any] | None:
    """(aspect, direction) of a functional consequence line, else None."""
    effects = []
    if row.get('activity'):
        effects.append((ACTIVITY_EFFECTS, row['activity'], 'activity'))
    if row.get('trafficking'):
        effects.append((TRAFFICKING_EFFECTS, row['trafficking'], 'trafficking'))
    if row.get('ptm'):
        effects.append((PTM_EFFECTS, row.get('ptm_type'), 'modification'))
    if not effects:
        return None
    if len(effects) > 1:
        raise ValueError(
            f'TRIP line {row.get("line")}: several functional consequences'
        )
    crosswalk, value, column = effects[0]
    if value not in crosswalk:
        raise ValueError(f'Unknown TRIP {column}: {row.get(column, value)!r}')
    return crosswalk[value]


def channel_form(row: dict[str, Any]) -> dict | None:
    """The modified residue of the TRP channel, if the line gives one.

    The modification column of a functional consequence line names the
    modification of the channel by the interactor, and for most lines the
    residue (``Phosphorylation (Tyr-273)``). The position is the source's; the
    line names no sequence it refers to, so the coordinate has no reference.
    """
    if not row.get('ptm') or not row.get('ptm_residue'):
        return None
    amino_acid = row['ptm_amino_acid']
    if amino_acid not in AMINO_ACIDS:
        raise ValueError(f'Unknown TRIP amino acid: {row["ptm_residue"]!r}')
    position = row['ptm_position']
    if not position:
        # A position of a sequence starts at 1.
        return None
    term = PTM_RESIDUE_TERMS[row['ptm_type']]
    return {
        'modifications': [
            {
                'term': term,
                'residue': AMINO_ACIDS[amino_acid],
                'position': position,
                'end_position': position,
                'coordinate_reference': {
                    'identifier': None,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                'description': row['ptm_type'] if term is None else None,
            }
        ],
    }


def _participant(row: dict[str, Any], role: str) -> Entity | None:
    name = row.get(f'{role}_name')
    if not name:
        return None
    identifiers = []
    if row.get(f'{role}_uniprot'):
        identifiers.append(
            Identifier(type=Namespace.UNIPROT, value=row[f'{role}_uniprot'])
        )
    identifiers.append(Identifier(type=Namespace.NAME, value=name))
    taxon = row.get(f'{role}_taxon')
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=(
            [Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon}')]
            if taxon
            else []
        ),
    )


def map_trip_observation(row: dict[str, Any]) -> Relation | None:
    """One observation as a relation from the interactor to the TRP channel.

    Screening, validation and characterization lines are evidence of a
    physical interaction. A functional consequence line states an effect of
    the interactor on the channel; the residue of a modification is a feature
    of the channel.
    """
    if trip_exclusion(row):
        return None
    interactor = _participant(row, 'interactor')
    channel = _participant(row, 'channel')
    if interactor is None or channel is None:
        return None
    method = row.get('method')
    if method not in METHOD_TERMS:
        raise ValueError(f'Unknown TRIP method: {method!r}')
    annotations = []
    effect = functional_effect(row)
    # The effect is checked first: its modification type is known.
    channel = channel._replace(molecular_form=channel_form(row))
    if effect:
        aspect, direction = effect
        if direction is not None:
            annotations.append(
                Annotation(term=slots.object_direction_qualifier, value=direction)
            )
        if aspect is not None:
            annotations.append(
                Annotation(term=slots.object_aspect_qualifier, value=aspect)
            )
    if METHOD_TERMS[method]:
        annotations.append(
            Annotation(term=slots.has_evidence_of_type, value=METHOD_TERMS[method])
        )
    if curie := pmid_curie(row.get('reference')):
        annotations.append(Annotation(term=slots.publications, value=curie))
    return Relation(
        subject=interactor,
        predicate=slots.affects if effect else slots.interacts_with,
        object=channel,
        annotations=annotations,
    )


# The site's list of proteins with their UniProt accessions, last capture.
_proteins_download = Download(
    url='https://web.archive.org/web/20210418195036id_/http://trpchannel.org/proteins',
    filename='trip_proteins_20210418.html',
    subfolder='trip',
    large=False,
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            # The export of 6 August 2015, the last one the site offered.
            url='https://web.archive.org/web/20171217101001id_/http://www.trpchannel.org:80/20150806.csv',
            filename='trip_20150806.csv',
            subfolder='trip',
            large=False,
            encoding='utf-8-sig',
        ),
        mapper=map_trip_observation,
        raw_parser=functools.partial(
            iter_trip_observations, proteins_download=_proteins_download
        ),
    ),
)
