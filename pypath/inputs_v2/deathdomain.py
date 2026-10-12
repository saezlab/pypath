"""
Parse the DeathDomain database pages and emit physical interaction records.

The DeathDomain database curated the interactions among proteins of the death
domain superfamily. It had one page per subfamily (CARD, DD, DED, PYD), each
with four tables: "In vitro interaction", "In vivo interaction",
"Characterization" and "Functional Role". A table line is one observation from
one publication: two proteins, a method, the details of that table and a
PubMed ID. Proteins are written with the names the database uses ("Apaf1",
"Caspase9", "Bcl-10"); the tables have no accessions.

Every observation between two different proteins is listed twice, once from
each protein; the second listing is not mapped. "Functional Role" lines give
the cellular outcome of an interaction (apoptosis, NF-kB activation), not an
observation of the interaction, and are not mapped either.

deathdomain.org now hosts another site, so the pages are read from the
Internet Archive (captures of April 2016; the database was last updated in
November 2012).
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.deathdomain import iter_deathdomain_observations
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

# Subfamily pages and the Internet Archive capture of each.
FAMILY_CAPTURES = {
    'CARD': '20160429003344',
    'DD': '20160428183051',
    'DED': '20160429144820',
    'PYD': '20160429003400',
}

# Species labels -> NCBI Taxonomy.
SPECIES_TAXA = {
    'Human': '9606',
    'Mouse': '10090',
    'Rat': '10116',
}
# Values of the species columns that name no single species: missing values,
# groups of species and the cell lysates an endogenous protein came from.
_NO_SPECIES = frozenset(
    {
        '',
        'Not specified',
        'Not used',
        'Murine',
        'Rodent',
        'Drosophila',
        'Herpesvirus2',
        'H460cell lyaste',
        'HEK293 lysate',
        'HeLa lysate',
        'HuH-7 lysate',
        'Jurkat lysate',
        'THP-1 lysate',
        'U937 lysate',
    }
)

# Method labels of the three interaction tables -> PSI-MI interaction
# detection method, where PSI-MI has the method under this name or as a
# synonym. A label mapped to None is known but has no such term.
_ASSEMBLY = 'In vitro purification protein assembly'
METHOD_TERMS = {
    'Co-immunoprecipitation': 'MI:0019',  # coimmunoprecipitation
    f'{_ASSEMBLY}(immunoprecipitation)': 'MI:0019',
    'Yeast two-hybrid': 'MI:0018',  # two hybrid
    'Two-hybrid': 'MI:0018',
    'GST fusion protein pull-down': 'MI:0096',  # pull down
    'His fusion protein pull-down': 'MI:0096',
    'His Fusion protein pull-down': 'MI:0096',
    f'{_ASSEMBLY}(Pull down)': 'MI:0096',
    'Affinity purification-mass spectrometry': 'MI:0004',  # affinity chromatography technology
    'Fluorescence resonance energy transfer': 'MI:0055',
    'Bimolecular fluorescence complementation': 'MI:0809',
    'Co-Immunofluorescence staining': 'MI:0403',  # colocalization
    'Size-exclusion chromatography': 'MI:0071',  # molecular sieving
    'Gel filtration': 'MI:0071',
    f'{_ASSEMBLY}(Size-exclusion chromatography)': 'MI:0071',
    'Surface plasmon resonance': 'MI:0107',
    f'{_ASSEMBLY}(Surface plasmon resonance)': 'MI:0107',
    'Isothermal titration calorimety': 'MI:0065',  # sic
    f'{_ASSEMBLY}(Isothermal titration calorimety)': 'MI:0065',
    f'{_ASSEMBLY}(chemical crosslinking)': 'MI:0030',  # cross-linking study
    f'{_ASSEMBLY}(Native gel-electrophoresis)': None,
    f'{_ASSEMBLY}(2D gel electrophoresis)': None,
    'Structure': None,
    'Mass spectrometry, MALS': None,
    'Fluorescence probe labeling': None,
    '': None,
}

EXCLUDED_FUNCTIONAL_ROLE = 'functional role of an interaction'
EXCLUDED_MIRRORED = 'second listing of an observation'

config = ResourceConfig(
    id=ResourceCv.DEATHDOMAIN,
    name='DeathDomain',
    full='Death Domain database',
    url='http://www.deathdomain.org/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.NAR_OPEN_ACCESS,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22135292',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Literature-curated physical interactions among death domain '
        'superfamily proteins, one record per method and publication.'
    ),
    description=(
        'The Death Domain database is a manually curated collection of '
        'protein-protein interactions of the death domain superfamily: '
        'proteins with a death domain (DD), death effector domain (DED), '
        'caspase recruitment domain (CARD) or pyrin domain (PYD), which '
        'assemble the signaling complexes of apoptosis, necrosis and '
        'inflammation. It was maintained at Yeungnam University, Seoul '
        'National University and Myongji University; the original website '
        'is no longer online.'
    ),
)


def _family_download(family: str) -> Download:
    return Download(
        url=(
            f'https://web.archive.org/web/{FAMILY_CAPTURES[family]}id_/'
            f'http://deathdomain.org/proteins/show?family={family}'
        ),
        filename=f'deathdomain_{family}_{FAMILY_CAPTURES[family]}.html',
        subfolder='deathdomain',
        large=False,
    )


_family_downloads = {family: _family_download(family) for family in FAMILY_CAPTURES}


def deathdomain_exclusion(row: dict[str, Any]) -> str | None:
    """Why an observation is not mapped, or None if it is."""
    if row.get('table') == 'functional_role':
        return EXCLUDED_FUNCTIONAL_ROLE
    if row.get('mirrored'):
        return EXCLUDED_MIRRORED
    return None


def _taxon(row: dict[str, Any], side: str) -> str | None:
    label = row.get(f'{side}_species')
    if label is None or label in _NO_SPECIES:
        return None
    if label not in SPECIES_TAXA:
        raise ValueError(f'Unknown DeathDomain species: {label!r}')
    return SPECIES_TAXA[label]


def _protein(row: dict[str, Any], side: str) -> Entity | None:
    name = (row.get(f'protein_{side}') or '').strip()
    if not name:
        return None
    taxon = _taxon(row, side)
    return Entity(
        type=Protein,
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=(
            [Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon}')]
            if taxon
            else []
        ),
    )


def map_deathdomain_observation(row: dict[str, Any]) -> Relation | None:
    """One observation of a physical interaction between two proteins."""
    if deathdomain_exclusion(row):
        return None
    first = _protein(row, 'a')
    second = _protein(row, 'b')
    if first is None or second is None:
        return None
    method = row.get('method')
    if method not in METHOD_TERMS:
        raise ValueError(f'Unknown DeathDomain method: {method!r}')
    annotations = []
    if METHOD_TERMS[method]:
        annotations.append(
            Annotation(term=slots.has_evidence_of_type, value=METHOD_TERMS[method])
        )
    reference = (row.get('reference') or '').strip()
    if curie := pmid_curie(reference):
        annotations.append(Annotation(term=slots.publications, value=curie))
    return Relation(
        subject=first,
        predicate=slots.interacts_with,
        object=second,
        annotations=annotations,
    )


resource = Resource(
    config,
    interactions=Dataset(
        # The four subfamily pages make one table; the parser opens them.
        download=None,
        mapper=map_deathdomain_observation,
        raw_parser=functools.partial(
            iter_deathdomain_observations, family_downloads=_family_downloads
        ),
    ),
)
