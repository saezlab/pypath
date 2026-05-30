"""
Parse TCDB (Transporter Classification Database) data and emit Entity records.

Two datasets are provided:

- transporters: proteins annotated with their TC number and transporter family
  name, derived from acc2tcid.py (UniProt → TCID) joined with families.py
  (family prefix → family name).

- substrates: transporter-substrate interactions linking UniProt accessions to
  CHEBI-identified small molecules, derived from getSubstrates.py (TCID →
  substrates) joined with acc2tcid.py (TCID → UniProt).
"""

from __future__ import annotations

import functools

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers import tcdb as _parsers
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    MoleculeAnnotationsCv,
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
)


config = ResourceConfig(
    id=ResourceCv.TCDB,
    name='TCDB',
    url='https://www.tcdb.org/',
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33170213',
    primary_category='annotations',
    description=(
        'The Transporter Classification Database (TCDB) provides a '
        'comprehensive classification system for membrane transport proteins '
        'based on both function and phylogeny. It assigns TC numbers to '
        'transporters and curates their substrates using CHEBI identifiers.'
    ),
)

f = FieldConfig(
    extract={
        'chebi': r'CHEBI:(\d+)',
    },
)

_transporters_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('uniprot')),
        CV(term=IdentifierNamespaceCv.TCDB, value=f('tcid')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.PROTEIN_FAMILY, value=f('family_name')),
    ),
)

_transport_schema = EntityBuilder(
    entity_type=EntityTypeCv.TRANSPORT,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.TCDB, value=f('tcid')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('transporter_uniprot')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.CHEMICAL,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.CHEBI,
                        value=f('substrate_id', extract='chebi'),
                    ),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('substrate_name')),
                ),
            ),
        ),
    ),
)

resource = Resource(
    config,
    transporters=Dataset(
        download=(_acc2tc := Download(
            url='http://www.tcdb.org/cgi-bin/projectv/public/acc2tcid.py',
            filename='tcdb_acc2tc.tsv',
            subfolder='tcdb',
            ext='tsv',
        )),
        mapper=_transporters_schema,
        raw_parser=functools.partial(
            _parsers.transporters,
            families_download=Download(
                url='http://www.tcdb.org/cgi-bin/projectv/public/families.py',
                filename='tcdb_families.html',
                subfolder='tcdb',
                ext='html',
            ),
        ),
    ),
    transports=Dataset(
        download=Download(
            url='https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py',
            filename='tcdb_substrates.tsv',
            subfolder='tcdb',
            ext='tsv',
        ),
        mapper=_transport_schema,
        raw_parser=functools.partial(
            _parsers.transport_substrates,
            acc2tc_download=_acc2tc,
        ),
    ),
)
