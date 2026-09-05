"""
Parse MRClinksDB metabolite-receptor interaction data and emit Entity records.

MRClinksDB provides curated metabolite–ligand-receptor (L-R) interactions
with compound classification (ClassyFire taxonomy), receptor UniProt IDs,
and supporting literature references.
"""

from __future__ import annotations
import re
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    MacromolecularComplex,
    ChemicalEntity,
    slots,
)
from omnipath_core.naming import Namespace
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    MembersFromList,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mrclinksdb import iter_mrclinksdb_interactions

config = ResourceConfig(
    id=ResourceCv.MRCLINKSDB,
    name='MRClinksDB',
    url='https://www.cellknowledge.com.cn/mrclinkdb/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    primary_category='interactions',
    pubmed='38978014',
    description='MRClinksDB is a curated database of metabolite–receptor interactions, covering ligand-receptor pairs with compound classification based on the ClassyFire taxonomy, receptor UniProt identifiers, and PubMed literature support.',
)
f = FieldConfig(
    transform={
        'pubchem_cid': lambda v: [
            token.strip()[4:]
            for token in v.split(';')
            if re.fullmatch('CID:\\d+', token.strip())
        ],
        'pubchem_sid': lambda v: [
            token.strip()[4:]
            for token in v.split(';')
            if re.fullmatch(r'SID:\d+', token.strip())
        ],
        'primary_name': lambda v: v.split('(')[0].strip(),
        'alt_names': lambda v: re.findall('\\(([^)]+)\\)', v),
    }
)


def _if_single(field: str):
    """Return a callable selector that yields *field* only for single-protein receptors.

    Gates on the absence of ``_`` in ``receptor_uniprot_id``. When the row
    represents a heteromeric complex the selector returns ``None``, causing
    the enclosing ``EntityBuilder`` to produce no identifiers and be silently
    skipped.

    The returned callable can be passed directly to ``f()`` alongside any
    transform, e.g. ``f(_if_single('protein_name'), transform='primary_name')``.
    """
    return (
        lambda row, _f=field: row.get(_f)
        if '_' not in str(row.get('receptor_uniprot_id') or '')
        else None
    )


def _if_complex(field: str):
    """Return a callable selector that yields *field* only for heteromeric complex receptors.

    Gates on the presence of ``_`` in ``receptor_uniprot_id``. When the row
    represents a single-protein receptor the selector returns ``None``, causing
    the enclosing ``EntityBuilder`` to produce no identifiers and be silently
    skipped.

    The returned callable can be passed directly to ``f()`` alongside any
    transform, e.g. ``f(_if_complex('protein_name'), transform='primary_name')``.
    """
    return (
        lambda row, _f=field: row.get(_f)
        if '_' in str(row.get('receptor_uniprot_id') or '')
        else None
    )


def _iter_with_taxon(opener, *, taxon_id: str, **kwargs):
    for row in iter_mrclinksdb_interactions(opener, **kwargs):
        row['taxon_id'] = taxon_id
        yield row


metabolite_builder = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.HMDB, value=f('hmdb_id')),
        CV(
            term=Namespace.PUBCHEM,
            value=f('pubchem_cid_sid', transform='pubchem_cid'),
        ),
        CV(
            term=Namespace.PUBCHEM_SUBSTANCE,
            value=f('pubchem_cid_sid', transform='pubchem_sid'),
        ),
        CV(term=Namespace.SMILES, value=f('canonical_smiles')),
        CV(term=Namespace.NAME, value=f('metabolite_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('molecular_formula'))
    ),
)
single_receptor_builder = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f(_if_single('receptor_uniprot_id'))),
        CV(term=Namespace.ENTREZ, value=f(_if_single('receptor_gene_id'))),
        CV(term=Namespace.GENESYMBOL, value=f(_if_single('receptor_symbol'))),
        CV(
            term=Namespace.NAME,
            value=f(_if_single('protein_name'), transform='primary_name'),
        ),
        CV(
            term=Namespace.SYNONYM,
            value=f(_if_single('protein_name'), transform='alt_names'),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f('taxon_id', transform=lambda v: f'NCBITaxon:{v}'),
        )
    ),
)
complex_receptor_builder = EntityBuilder(
    entity_type=MacromolecularComplex,
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.NAME,
            value=f(_if_complex('protein_name'), transform='primary_name'),
        ),
        CV(
            term=Namespace.SYNONYM,
            value=f(_if_complex('protein_name'), transform='alt_names'),
        ),
    ),
    annotations=AnnotationsBuilder(),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=Protein,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.UNIPROT,
                    value=f('receptor_uniprot_id', delimiter='_'),
                ),
                CV(
                    term=Namespace.ENTREZ,
                    value=f('receptor_gene_id', delimiter='_'),
                ),
                CV(
                    term=Namespace.GENESYMBOL,
                    value=f('receptor_symbol', delimiter='_'),
                ),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term=slots.in_taxon,
                    value=f('taxon_id', transform=lambda v: f'NCBITaxon:{v}'),
                )
            ),
        )
    ),
)


def _mrclinksdb_receptor_builder(row: dict[str, object]) -> object | None:
    if '_' in str(row.get('receptor_uniprot_id') or ''):
        return complex_receptor_builder.build(row)
    return single_receptor_builder.build(row)


interactions_schema = RelationBuilder(
    subject=metabolite_builder,
    predicate=slots.interacts_with,
    object=_mrclinksdb_receptor_builder,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.MRCLINKSDB, value=f('mrid'))
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.supporting_data_source,
            value=f('other_db', delimiter=';'),
        ),
        CV(
            term=slots.publications,
            value=f(
                'pmid',
                delimiter=';',
                transform=lambda v: f'PMID:{v.removeprefix("PMID:")}',
            ),
        ),
    ),
)


def _transporter_schema(species: str) -> RelationBuilder:
    fields = {
        'human': ('uniprot_id', 'human_geneid', 'gene_name'),
        'mouse': ('mouse_uniprot', 'mouse_geneid', 'mouse_gene_symbol'),
    }[species]
    return RelationBuilder(
        subject=EntityBuilder(
            entity_type=Protein,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.UNIPROT, value=f(fields[0])),
                CV(term=Namespace.ENTREZ, value=f(fields[1])),
                CV(term=Namespace.GENESYMBOL, value=f(fields[2])),
                CV(
                    term=Namespace.NAME,
                    value=f('enzyme_name', transform='primary_name'),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.in_taxon,
                    value=f('taxon_id', transform=lambda v: f'NCBITaxon:{v}'),
                )
            ),
        ),
        predicate=slots.associated_with,
        object=EntityBuilder(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.HMDB, value=f('hmdb_id')),
                CV(term=Namespace.NAME, value=f('metabolite_name')),
            ),
            annotations=AnnotationsBuilder(),
        ),
    )


human_transporters_schema = _transporter_schema('human')
mouse_transporters_schema = _transporter_schema('mouse')
resource = Resource(
    config,
    human_interactions=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Homo%20sapiens%20metabolite%20L-R%20interaction.txt',
            filename='mrclinksdb_human_interactions.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=interactions_schema,
        raw_parser=lambda opener, **kwargs: _iter_with_taxon(
            opener, taxon_id='9606', **kwargs
        ),
    ),
    mouse_interactions=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Mus%20musculus%20metabolite%20L-R%20interaction.txt',
            filename='mrclinksdb_mouse_interactions.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=interactions_schema,
        raw_parser=lambda opener, **kwargs: _iter_with_taxon(
            opener, taxon_id='10090', **kwargs
        ),
    ),
    human_transporters=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Homo%20sapiens%20transporter%20protein.txt',
            filename='mrclinksdb_human_transporters.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=human_transporters_schema,
        raw_parser=lambda opener, **kwargs: _iter_with_taxon(
            opener, taxon_id='9606', **kwargs
        ),
    ),
    mouse_transporters=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Mus%20musculus%20transporter%20protein.txt',
            filename='mrclinksdb_mouse_transporters.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=mouse_transporters_schema,
        raw_parser=lambda opener, **kwargs: _iter_with_taxon(
            opener, taxon_id='10090', **kwargs
        ),
    ),
)
