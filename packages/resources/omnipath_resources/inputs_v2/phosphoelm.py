"""
Parse the Phospho.ELM phosphorylation sites with their kinases and emit Relation records.

Source file: the dump of the vertebrate data (``phosphoELM_vertebrate_latest.dump.tgz``,
April 2015), a tab separated table read by ``parsers/phosphoelm.py``. Each line
is one observation: a phosphorylated residue of a protein, the kinase named for
it (on most lines none), the publication, whether the study was low-throughput
(LTP) or high-throughput (HTP), and the species of the protein.

    acc     sequence  position  code  pmids     kinases  source  species       entry_date
    P00519  MGQQP...  412       Y     10964922  Abl      LTP     Homo sapiens  2005-03-14 ...

Only the lines that name a kinase are enzyme-substrate relations; the other
phosphorylation sites have no enzyme and are not read.

The kinases are Phospho.ELM's own names, not gene symbols. A name ending in
``_group`` ("PKA_group", "CK2_group") stands for a family of kinases and is a
``ProteinFamily``; the table does not say which member of it phosphorylates the
site, and the family is not replaced by its members. The names are not translated
to other identifiers. Where the kinase list on the Phospho.ELM site
(http://phospho.elm.eu.org/kinases.html) gives a name a UniProt accession, that
accession is an identifier of the kinase: it is Phospho.ELM's own assignment.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.molecular_forms import molecular_form_from_identifiers
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._molecular_forms import combine_forms
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phosphoelm import iter_phosphoelm
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

# The species of the dump, by the Latin name it uses, with the NCBI Taxonomy ID.
# The reader rejects any other name.
ORGANISMS = {
    'Homo sapiens': '9606',
    'Mus musculus': '10090',
    'Rattus norvegicus': '10116',
    'Bos taurus': '9913',
    'Gallus gallus': '9031',
    'Oryctolagus cuniculus': '9986',
    'Xenopus laevis': '8355',
    'Mesocricetus auratus': '10036',
    'Sus scrofa': '9823',
    'Ovis aries': '9940',
    # NCBI Taxonomy 9615 is "Canis lupus familiaris"; "Canis familiaris" is a
    # synonym of it.
    'Canis familiaris': '9615',
    'Didelphis virginiana': '9267',
    'Cricetulus griseus': '10029',
}
# The source column of the dump: how the study that found the site was done.
# The value is an attribute of the evidence of the relation.
_THROUGHPUT = 'phosphoelm:throughput'
THROUGHPUT_VALUES = {
    'LTP': 'Low Throughput',
    'HTP': 'High Throughput',
}
# The accession column holds a UniProt accession, with the isoform number where
# the site is on an isoform, or an Ensembl protein.
_ACCESSION_NAMESPACES = {
    'uniprot': Namespace.UNIPROT,
    'ensembl_protein': Namespace.ENSP,
}


def _taxon(row: dict[str, Any]) -> str:
    # The reader has rejected names outside ORGANISMS.
    return f'NCBITaxon:{ORGANISMS[row["species"]]}'


def _throughput(row: dict[str, Any]) -> str:
    value = row.get('source')
    if value not in THROUGHPUT_VALUES:
        raise ValueError(f'Unknown Phospho.ELM throughput: {value!r}')
    return THROUGHPUT_VALUES[value]


def _substrate_accession(namespace: Namespace):
    def value(row: dict[str, Any]) -> str | None:
        same = _ACCESSION_NAMESPACES[row['accession_database']] == namespace
        return row['accession'] if same else None

    return value


def _substrate_reference(row: dict[str, Any]) -> dict[str, str]:
    return {
        'ns': str(_ACCESSION_NAMESPACES[row['accession_database']]),
        'id': row['accession'],
    }


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The phosphorylated residue, a feature of the substrate in this record.

    The position is the source's, in the sequence of the accession; the reader
    has checked that the residue is at that position of the sequence in the
    dump. A UniProt accession with an isoform number, and an Ensembl protein,
    also name the sequence the substrate is, which the form keeps.
    """
    reference = _substrate_reference(row)
    site = {
        'modifications': [
            {
                'term': 'phosphorylated residue',
                'residue': row['code'],
                'position': row['position'],
                'end_position': row['position'],
                'coordinate_reference': {
                    'identifier': reference,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
            }
        ],
    }
    return combine_forms(molecular_form_from_identifiers([reference]), site)


config = ResourceConfig(
    id=ResourceCv.PHOSPHOELM,
    name='Phospho.ELM',
    url='http://phospho.elm.eu.org/',
    # The data are under the Phospho.ELM academic license agreement (EMBLEM;
    # academic research only):
    # http://phospho.elm.eu.org/dumps/Phospho.Elm_AcademicLicense.pdf
    # The download page (http://phospho.elm.eu.org/dataset.html) says: "The
    # data in phospho.ELM should not be used or shared for any commercial
    # purposes and should not be distributed to a third party without prior
    # consent."
    license=LicenseCV.NONCOMMERCIAL,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='21062810',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Kinase-substrate phosphorylation sites with their publications and '
        'low- or high-throughput origin, April 2015 dump.'
    ),
    description=(
        'Phospho.ELM is a database of experimentally verified serine, '
        'threonine and tyrosine phosphorylation sites of eukaryotic '
        'proteins, curated from the literature and from high-throughput '
        'studies, and part of the Eukaryotic Linear Motif (ELM) resource. '
        'Each site is given with its position in the protein sequence, the '
        'publication, whether the evidence is from a low-throughput or a '
        'high-throughput study, and, for a part of the sites, the kinase '
        'that phosphorylates it.'
    ),
)

f = FieldConfig()

# The kinase by its Phospho.ELM name; a group name is a protein family. The
# UniProt accession is the one the kinase list of Phospho.ELM gives the name.
# The dump does not give the species of a kinase.
_enzyme = EntityBuilder(
    entity_type=lambda row: (
        model.ProteinFamily if row.get('enzyme_is_group') else model.Protein
    ),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('enzyme_uniprot')),
        CV(term=Namespace.NAME, value=f('enzyme')),
    ),
)

_substrate = EntityBuilder(
    entity_type=model.Protein,
    # The phosphorylated residue is in the substrate's molecular form.
    molecular_form=_substrate_form,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=_substrate_accession(Namespace.UNIPROT)),
        CV(term=Namespace.ENSP, value=_substrate_accession(Namespace.ENSP)),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_taxon)),
)

enzyme_substrate_schema = RelationBuilder(
    subject=_enzyme,
    predicate=slots.affects,
    object=_substrate,
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
        # A line has one publication, or none ("N.N.").
        CV(
            term=slots.publications,
            value=f('pmid', transform=lambda pmid: f'PMID:{pmid}'),
        ),
        CV(term=_THROUGHPUT, value=_throughput),
    ),
)

resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url=(
                'http://phospho.elm.eu.org/dumps/'
                'phosphoELM_vertebrate_latest.dump.tgz'
            ),
            filename='phosphoELM_vertebrate_latest.dump.tgz',
            subfolder='phosphoelm',
            large=True,
            ext='tar.gz',
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=functools.partial(
            iter_phosphoelm,
            organisms=frozenset(ORGANISMS),
            # The kinase names with the accessions Phospho.ELM assigns to them.
            kinases_download=Download(
                url='http://phospho.elm.eu.org/kinases.html',
                filename='kinases.html',
                subfolder='phosphoelm',
                large=False,
            ),
        ),
    ),
)
