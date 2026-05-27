"""Rhea reaction database — inputs_v2 module.

Provides a :class:`~pypath.inputs_v2.base.Resource` with four datasets parsed
from the Rhea web API and FTP export:

Datasets:
    reactions: REACTION entities for all master reactions (metabolic and
        transport).  ChEBI participants are included as SMALL_MOLECULE
        sub-members with REACTANT/PRODUCT role and SUBCELLULAR_LOCATION
        annotations derived from the equation.  Cross-references to EC,
        PubMed, GO, EcoCyc, MetaCyc, KEGG, and Reactome are stored as
        annotations.
    metabolic_reactions: As ``reactions`` but restricted to metabolic
        (non-transport) reactions only.
    transport_reactions: TRANSPORT entities for reactions that move a
        molecule across a compartment boundary.  Identified by the presence
        of ``[compartment]`` suffixes in the equation string.  Compartment
        labels are stored as SUBCELLULAR_LOCATION annotations on each
        participant.
    catalysis: CATALYSIS entities linking a UniProt PROTEIN as CONTROLLER
        to a Rhea master REACTION as CONTROLLED.  Reaction direction is
        derived from the DIRECTION column of the rhea2uniprot FTP export
        (LR → LEFT-TO-RIGHT, RL → RIGHT-TO-LEFT, BI → REVERSIBLE).
"""

from __future__ import annotations

import re

from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.rhea import _raw


config = ResourceConfig(
    id = ResourceCv.RHEA,
    name = 'Rhea',
    url = 'https://www.rhea-db.org/',
    license = LicenseCV.CC_BY_4_0,
    update_category = UpdateCategoryCV.REGULAR,
    pubmed = '33479744',
    primary_category = 'reactions',
    description = (
        'Rhea is an expert-curated knowledgebase of chemical and transport '
        'reactions of biological interest.'
    ),
)


_direction_map = {
    'LR': 'LEFT-TO-RIGHT',
    'RL': 'RIGHT-TO-LEFT',
    'BI': 'REVERSIBLE',
    'UN': None,
}

_role_map = {
    'reactant': BiologicalRoleCv.REACTANT,
    'product': BiologicalRoleCv.PRODUCT,
}

_CHEBI_RE = re.compile(r'(?:CHEBI:)?(\d+)')
_RHEA_ID_RE = re.compile(r'(?:RHEA:)?(\d+)')


# ── reactions ─────────────────────────────────────────────────────────────────

f = FieldConfig(
    delimiter = ';',
    map = {
        'role': lambda value: _role_map.get(value),
    },
)


reactions_download = Download(
    url = (
        'https://www.rhea-db.org/rhea/?query=&columns=rhea-id,equation,chebi,'
        'chebi-id,ec,uniprot,go,pubmed,reaction-xref(EcoCyc),reaction-xref(MetaCyc),'
        'reaction-xref(KEGG),reaction-xref(Reactome),reaction-xref(M-CSA)'
        '&format=tsv&limit=1000000'
    ),
    filename = 'rhea_reactions.tsv',
    subfolder = 'rhea',
    ext = 'tsv',
    default_mode = 'r',
)


reactions_schema = EntityBuilder(
    entity_type = EntityTypeCv.REACTION,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.RHEA_ID, value = f('rhea_id')),
        CV(term = IdentifierNamespaceCv.NAME, value = f('equation')),
    ),
    annotations = AnnotationsBuilder(
        CV(term = MoleculeAnnotationsCv.EC_NUMBER, value = f('ec')),
        CV(term = IdentifierNamespaceCv.PUBMED, value = f('pubmed')),
        CV(term = IdentifierNamespaceCv.CV_TERM_ACCESSION, value = f('go')),
        CV(term = IdentifierNamespaceCv.ECOCYC, value = f('ecocyc')),
        CV(term = IdentifierNamespaceCv.METACYC, value = f('metacyc')),
        CV(term = IdentifierNamespaceCv.KEGG, value = f('kegg')),
        CV(term = IdentifierNamespaceCv.REACTOME_STABLE_ID, value = f('reactome')),
    ),
    membership = MembershipBuilder(
        MembersFromList(
            entity_type = EntityTypeCv.SMALL_MOLECULE,
            identifiers = IdentifiersBuilder(
                CV(
                    term = IdentifierNamespaceCv.CHEBI,
                    value = f('participant_chebi', delimiter = '||', extract = _CHEBI_RE),
                ),
                CV(
                    term = IdentifierNamespaceCv.NAME,
                    value = f('participant_display_name', delimiter = '||'),
                ),
            ),
            annotations = AnnotationsBuilder(
                CV(term = f('participant_role', delimiter = '||', map = 'role')),
            ),
        ),
    ),
)


# ── transport_reactions ───────────────────────────────────────────────────────

transport_reactions_schema = EntityBuilder(
    entity_type = EntityTypeCv.TRANSPORT,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.RHEA_ID, value = f('rhea_id')),
        CV(term = IdentifierNamespaceCv.NAME, value = f('equation')),
    ),
    annotations = AnnotationsBuilder(
        CV(term = MoleculeAnnotationsCv.EC_NUMBER, value = f('ec')),
        CV(term = IdentifierNamespaceCv.PUBMED, value = f('pubmed')),
        CV(term = IdentifierNamespaceCv.CV_TERM_ACCESSION, value = f('go')),
        CV(term = IdentifierNamespaceCv.ECOCYC, value = f('ecocyc')),
        CV(term = IdentifierNamespaceCv.METACYC, value = f('metacyc')),
        CV(term = IdentifierNamespaceCv.KEGG, value = f('kegg')),
        CV(term = IdentifierNamespaceCv.REACTOME_STABLE_ID, value = f('reactome')),
    ),
    membership = MembershipBuilder(
        MembersFromList(
            entity_type = EntityTypeCv.SMALL_MOLECULE,
            identifiers = IdentifiersBuilder(
                CV(
                    term = IdentifierNamespaceCv.CHEBI,
                    value = f('participant_chebi', delimiter = '||', extract = _CHEBI_RE),
                ),
                CV(
                    term = IdentifierNamespaceCv.NAME,
                    value = f('participant_display_name', delimiter = '||'),
                ),
            ),
            annotations = AnnotationsBuilder(
                CV(term = f('participant_role', delimiter = '||', map = 'role')),
                CV(
                    term = ParticipantMetadataCv.MEMBRANE_SIDE,
                    value = f('participant_compartment', delimiter = '||'),
                ),
            ),
        ),
    ),
)


# ── catalysis ─────────────────────────────────────────────────────────────────

g = FieldConfig(
    map = {
        'direction': lambda value: _direction_map.get(value),
    },
)


catalysis_download = Download(
    url = 'https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv',
    filename = 'rhea2uniprot.tsv',
    subfolder = 'rhea',
    ext = 'tsv',
    default_mode = 'r',
)


catalysis_schema = EntityBuilder(
    entity_type = EntityTypeCv.CATALYSIS,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.RHEA_ID, value = g('MASTER_ID', extract = _RHEA_ID_RE)),
    ),
    annotations = AnnotationsBuilder(
        CV(
            term = InteractionMetadataCv.CONVERSION_DIRECTION,
            value = g('DIRECTION', map = 'direction'),
        ),
    ),
    membership = MembershipBuilder(
        Member(
            entity = EntityBuilder(
                entity_type = EntityTypeCv.PROTEIN,
                identifiers = IdentifiersBuilder(
                    CV(term = IdentifierNamespaceCv.UNIPROT, value = g('ID')),
                ),
            ),
            annotations = AnnotationsBuilder(CV(term = BiologicalRoleCv.CONTROLLER)),
        ),
        Member(
            entity = EntityBuilder(
                entity_type = EntityTypeCv.REACTION,
                identifiers = IdentifiersBuilder(
                    CV(
                        term = IdentifierNamespaceCv.RHEA_ID,
                        value = g('MASTER_ID', extract = _RHEA_ID_RE),
                    ),
                ),
            ),
            annotations = AnnotationsBuilder(CV(term = BiologicalRoleCv.CONTROLLED)),
        ),
    ),
)


# ── resource ─────────────────────────────────────────────────────────────────

resource = Resource(
    config,
    reactions = Dataset(
        download = reactions_download,
        mapper = reactions_schema,
        raw_parser = _raw,
    ),
    metabolic_reactions = Dataset(
        download = reactions_download,
        mapper = reactions_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'metabolic_reactions', **kwargs),
    ),
    transport_reactions = Dataset(
        download = reactions_download,
        mapper = transport_reactions_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'transport_reactions', **kwargs),
    ),
    catalysis = Dataset(
        download = catalysis_download,
        mapper = catalysis_schema,
        raw_parser = iter_tsv,
    ),
)
