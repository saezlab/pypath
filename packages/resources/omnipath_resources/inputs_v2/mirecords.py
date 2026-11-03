"""
Parse miRecords microRNA - target gene interactions and emit Relation records.

miRecords is a literature-curated database of experimentally validated
interactions between animal microRNAs and their target genes. Its table lists
one literature account per line: the publication, the miRNA and the species it
is from, the target gene with its RefSeq accession and species, and the
experiments (reporter assays, methods, whether the target's mRNA level changed).

Each account becomes a relation, the miRNA ``affects`` the target gene. The
table does not state whether the miRNA raised or lowered the target: its columns
name how the miRNA level was changed, the methods and whether the mRNA level
"changed" or "unchanged", never the direction, and the descriptions are the
papers' prose. So the record has no direction and no aspect.

The miRNA is named as the table writes it, with no accession. In 221 rows the
name is between square brackets (``[hsa-miR-34a]``, ``[miR-129]``); the
documentation of the database does not explain them. They are on names that the
curators took from the paper as it gives them, often a family or precursor name
without a -3p/-5p arm or an organism prefix, and not a mature name of miRBase
(two of the rows say "The miRNA can not be found in miRBase"). The brackets
are not part of the name and are removed; the name is not otherwise changed or
translated between miRBase versions, and -3p and -5p forms stay separate.

A target whose gene name is missing or was turned into a date by Excel is
identified by its accession alone; with neither, the row gives no relation.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import Gene, MicroRNA, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mirecords import iter_interactions
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

config = ResourceConfig(
    id=ResourceCv.MIRECORDS,
    name='miRecords',
    # The site is offline; the file comes from the OmniPath rescue server.
    url='http://c1.accurascience.com/miRecords/',
    # The disclaimer of the site (archived at
    # https://web.archive.org/web/20161123115811/http://c1.accurascience.com/miRecords/copyright.php)
    # says: "you may only download the content for your own personal
    # non-commercial use", "The authorized use of the information on this
    # website is limited to informational and educational purposes only, and NOT
    # for operational or commercial purposes", and "You are not permitted to
    # copy, broadcast, download, store (in any medium), transmit, [...] adapt or
    # change in any way the content of these web pages for any other purpose
    # whatsoever without the prior written permission of Biolead.org". The
    # project owner decided to use the rescued copy of the download.
    license=LicenseCV.NONCOMMERCIAL,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='18996891',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Experimentally validated animal microRNA - target gene interactions '
        'curated from the literature (release of April 2013).'
    ),
    description=(
        'miRecords is a resource for animal microRNA - target interactions '
        'developed at Biolead.org. Its Validated Targets component is a '
        'literature-curated database of microRNA - target pairs with the '
        'evidence of each paper: whether it concerns the target gene, a '
        'region of the transcript or a single target site, how the microRNA '
        'level was changed, the reporter assays and other methods, and '
        'whether the target mRNA level changed. This module imports the last '
        'release of the Validated Targets table (27 April 2013), from the '
        'copy kept on the OmniPath rescue server because the site is offline. '
        'The Predicted Targets component is not part of the download.'
    ),
)

# Species names of the two species columns, after the repairs of the parser ->
# NCBI Taxonomy. 'human' is the source's own spelling in three rows.
TAXA = {
    'Bombyx mori': '7091',
    'Bos taurus': '9913',
    'Caenorhabditis elegans': '6239',
    'Ciona intestinalis': '7719',
    'Cricetulus griseus': '10029',
    'Danio rerio': '7955',
    'Drosophila melanogaster': '7227',
    'Gallus gallus': '9031',
    'Homo sapiens': '9606',
    'human': '9606',
    'Human immunodeficiency virus 1': '11676',
    'Kaposi sarcoma-associated herpesvirus': '37296',  # human gammaherpesvirus 8
    'Mus musculus': '10090',
    'Ovis aries': '9940',
    'Rattus norvegicus': '10116',
    'Trichomonas vaginalis': '5722',
    'Vesicular stomatitis Indiana virus': '11277',
    'Xenopus laevis': '8355',
    'Xenopus tropicalis': '8364',
}

# Cells of the species columns that do not name a species: they get no taxon.
# 'Xenopus' is a genus, 'us musculus' and 'Mus musculus miR-24' are damaged
# names, and '[Virus]' and the RefSeq accession are not names of organisms.
NO_SPECIES = frozenset({
    '[Virus]',
    'NM_213232.1',
    'Xenopus',
    'us musculus',
    'Mus musculus miR-24',
})


def _taxon(column: str):
    def curie(row: dict[str, Any]) -> str | None:
        name = row.get(column)
        if name in NO_SPECIES or name is None:
            return None
        if name not in TAXA:
            raise ValueError(f'Unknown miRecords species: {name!r}')
        return f'NCBITaxon:{TAXA[name]}'

    return curie


# The column is named for RefSeq accessions of the transcript (NM_, XM_, NR_, XR_)
# and holds some of the genome, contig and genomic region (NC_, NW_, NG_). The
# version is kept where the source has it. A few cells hold a GenBank accession
# instead (AF000673, BC083072.1); NM213137 is a RefSeq accession written without
# its underscore, which is not repaired.
_REFSEQ = r'^((?:NM|NR|XM|XR|NC|NG|NW)_\d+(?:\.\d+)?)$'
_GENBANK = r'^((?![NX][MR])[A-Z]{1,2}\d{5,6}(?:\.\d+)?)$'

f = FieldConfig(
    extract={
        'refseq': _REFSEQ,
        'genbank': _GENBANK,
        # Two rows have 0 where they have no publication.
        'pmid': PMID_EXTRACT,
    },
    transform={
        'pmid_curie': lambda v: f'PMID:{v}',
    },
)

# The miRNA is named by its mature name only, as the table writes it. The
# species is the species column's, whatever the prefix of the name.
mirna = EntityBuilder(
    entity_type=MicroRNA,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('mirna')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_taxon('mirna_species')),
    ),
)

# The gene name first, then the accession of the target's transcript.
target_gene = EntityBuilder(
    entity_type=Gene,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.GENESYMBOL, value=f('target_gene')),
        CV(term=Namespace.REFSEQ, value=f('target_refseq', extract='refseq')),
        CV(term=Namespace.GENBANK, value=f('target_refseq', extract='genbank')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_taxon('target_species')),
    ),
)

interactions_schema = RelationBuilder(
    subject=mirna,
    predicate=slots.affects,
    object=target_gene,
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=f('pmid', extract='pmid', transform='pmid_curie')),
    ),
)


def has_target_identifier(row: dict[str, Any]) -> bool:
    """Whether the target has a gene name or an accession of a known kind."""
    accession = row.get('target_refseq') or ''
    return bool(row.get('target_gene')) or any(
        re.match(pattern, accession) for pattern in (_REFSEQ, _GENBANK)
    )


def map_interaction(row: dict[str, Any]):
    """The relation of one account; None where the target has no identifier."""
    if not has_target_identifier(row):
        return None
    return interactions_schema(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/miRecords_v4.xls',
            filename='miRecords_v4.xls',
            subfolder='mirecords',
            large=True,
            default_mode='rb',
            encoding=None,
        ),
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
)
