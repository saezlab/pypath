"""
Parse BioGRID interactions and emit Relation records.

BioGRID publishes its curated interactions of all organisms in its own
tab-separated format, BioGRID Tab 3.0, with a header line:
https://wiki.thebiogrid.org/doku.php/biogrid_tab_version_3.0
One row is one interaction between two genes reported by one publication
with one experimental system. The experimental system type says whether the
evidence is physical or genetic. These are different assertions, so they are
separate datasets, which read the same file.

Each interaction also carries the experimental context the file states: the
ontology terms (cell line, tissue, disease, chemical, phenotype, ...) and their
qualifiers, the type of the phenotype terms (``biogrid:phenotype_type``) and
the modification a biochemical activity observed (``biogrid:modification``).
They are carried in two ways, by what a term says:

* A phenotype term of a genetic interaction is an association of both genes to
  the term as an entity (``associated_with`` an ontology class identified in
  the vocabulary of the term, as UniProt links a protein to its GO terms), so
  that the term can be searched and its genes listed. BioGRID states the
  phenotype for the pair of genes; ``associated_with`` asserts nothing causal.
* Every other term is an annotation of the interaction whose value is the
  accession (``biogrid:cell_line = BTO:0000007``), under one attribute per kind
  of term (the few phenotype terms of physical interactions too:
  ``biogrid:phenotype``, with no link), and every qualifier of a term an
  annotation ``biogrid:<kind>_qualifier``
  (``biogrid:cell_line_qualifier = DOID:3702``, also for the qualifiers of the
  phenotype terms that are links).

They describe the evidence, not a different interaction, so they are not
Biolink qualifiers and do not change the identity of the relation. The kind of
a term is BioGRID's category, repaired where the category is wrong (it labels
38,398 BTO and 1,014 CL terms of physical rows ``phenotype``): see ``_KINDS``.
The names of the terms are not carried: the ontologies name them.

The release is pinned: this module reads release 5.0.262 (the file is dated
29 September 2026) of BioGRID's release archive until the URL is changed, so
that the data do not change under the module when BioGRID publishes a new
release. The file has every organism; the organism of each participant is
the one the row names (a row can have participants of two organisms).
"""

from __future__ import annotations

import functools
from typing import Any, NamedTuple

from biolink_model.datamodel.model import Gene, OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2 import _mitab
from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.biogrid import (
    PHENOTYPES,
    QUALIFIERS_BY_KIND,
    SCORE,
    TERMS,
    TERMS_BY_KIND,
    THROUGHPUTS,
    UNIPROT_ACCESSION,
    iter_interactions,
)


config = ResourceConfig(
    id=ResourceCv.BIOGRID,
    name='BioGRID',
    full='Biological General Repository for Interaction Datasets',
    url='https://thebiogrid.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.MIT,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33070389',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Physical and genetic interactions of genes and proteins of all '
        'organisms curated by BioGRID from the literature (release '
        '5.0.262), with the experimental system, throughput, score and '
        'publication of each, the experimental context (cell line, tissue, '
        'disease and other ontology terms, with their qualifiers) and the '
        'modification observed, and the phenotypes of genetic interactions '
        'as associations of the genes to the phenotype terms.'
    ),
    description=(
        'BioGRID is a repository of protein, genetic and chemical '
        'interactions curated from the primary literature for the major '
        'model organisms and humans. Each interaction is recorded with the '
        'publication that reports it and the experimental system that '
        'supports it, from low-throughput studies and from high-throughput '
        'datasets. An interaction in BioGRID is evidence of an association '
        'between two genes or their products; the experimental system tells '
        'how direct that evidence is. This module reads the BioGRID Tab 3.0 '
        'file of all organisms of release 5.0.262, with the organism of '
        'each participant. It also keeps the ontology terms BioGRID '
        'records for an interaction. The phenotype terms of a genetic '
        'interaction (APO, HP, GO, PATO, MP and WBPhenotype terms) link both '
        'genes to the term, an ontology class, so that a phenotype can be '
        'searched and its genes shown. The other terms (the cell line, cell '
        'type or tissue of the experiment, the disease, the chemical, and '
        'the like) are attributes of the interaction whose value is the '
        'accession of the term, each with the accessions of its qualifiers, '
        'next to the type of a phenotype and the modification a biochemical '
        'activity observed. It does not read the free-text qualification '
        'column, the tag column or the names of the terms; BioGRID labels the '
        'category of a term inconsistently, so a few categories are repaired '
        '(see the module).'
    ),
)

# The release the module reads. The link is the one the BioGRID download
# page shows for the file (the page of the
# "downloads.thebiogrid.org/File/..." address of the release archive is an
# HTML page that opens it); the archive has one text member.
_RELEASE = '5.0.262'
_FILENAME = f'BIOGRID-ALL-{_RELEASE}.tab3.zip'
_URL = (
    'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/'
    f'BIOGRID-{_RELEASE}/{_FILENAME}'
)
# BioGRID interaction records have their own numbers, which are not the
# numbers of its gene records (``Namespace.BIOGRID``).
_BIOGRID_INTERACTION = 'biogrid_interaction'
# Attributes of the BioGRID ontologies (github.com/BioGRID/BioGRID-Ontologies):
# the experimental system and the throughput of an interaction. Their values
# are the names BioGRID gives.
_EXPERIMENTAL_SYSTEM = 'BIOGRID_SYSTEM:0000001'
_THROUGHPUT = 'BIOGRID_TP:0000001'
# A study of both throughputs has both values (joined by ``|`` in the file).
_THROUGHPUT_VALUES = {name: name for name in THROUGHPUTS}
# The modification of a biochemical activity row, as the file writes it (the
# 22 values of release 5.0.262, all on Biochemical Activity rows of the
# physical dataset). ``No Modification`` is a statement of the source and is
# carried as written. A value that is not here raises.
_MODIFICATION = 'biogrid:modification'
_MODIFICATIONS = frozenset({
    'Acetylation', 'Deacetylation', 'Demethylation', 'Dephosphorylation',
    'Desumoylation', 'Deubiquitination', 'Deneddylation', 'FAT10ylation',
    'Glycosylation', 'ISGylation', 'Methylation', 'Nedd(Rub1)ylation',
    'Neddylation', 'No Modification', 'Phosphorylation', 'Prenylation',
    'Proteolytic Processing', 'Ribosylation', 'Sumoylation',
    'Ubiquitination', 'Ufmylation', 'de-ISGylation',
})
# The type BioGRID gives a term of a genetic interaction, as written (``-``, no
# type, is not carried). It is the type of a phenotype term (also, in the
# file, of a few terms of other categories), so it is one attribute of the
# relation; a type that is not here raises.
_PHENOTYPE_TYPE = 'biogrid:phenotype_type'
_TERM_TYPES = frozenset({'partial rescue', 'undetermined', 'wild type'})
# The score of a row, with the scale of the publication that reports it: the
# scales differ between publications and the file names no unit. It is
# therefore a source attribute, not ``has_confidence_score``.
_SCORE = 'biogrid:score'


class _System(NamedTuple):
    """What an experimental system means for the mapped record."""

    # PSI-MI interaction detection method and interaction type.
    detection_method: str
    interaction_type: str
    # Class of the two participants.
    participant: type = Protein


# Experimental system -> PSI-MI terms. BioGRID's own PSI-MI TAB export gives
# each system one detection method and one interaction type; this table
# repeats them. These are the 31 systems of the all-organism file of release
# 5.0.262; a system that is not here raises.
#
# BioGRID identifies every participant by a gene record. The physical
# systems are defined between proteins, so their participants are the
# proteins of those genes. Two systems are between a protein and an RNA, and
# the file does not say which participant is the RNA: there the participants
# stay the genes.
_PHYSICAL_SYSTEMS: dict[str, _System] = {
    'Affinity Capture-Luminescence': _System('MI:0004', 'MI:0915'),
    'Affinity Capture-MS': _System('MI:0004', 'MI:0915'),
    'Affinity Capture-RNA': _System('MI:0004', 'MI:0915', Gene),
    'Affinity Capture-Western': _System('MI:0004', 'MI:0915'),
    'Biochemical Activity': _System('MI:0415', 'MI:0407'),
    'Co-crystal Structure': _System('MI:0114', 'MI:0914'),
    'Co-fractionation': _System('MI:0401', 'MI:0403'),
    'Co-localization': _System('MI:0428', 'MI:0403'),
    'Co-purification': _System('MI:0401', 'MI:0914'),
    'Cross-Linking-MS (XL-MS)': _System('MI:0030', 'MI:0915'),
    'FRET': _System('MI:0055', 'MI:0407'),
    'Far Western': _System('MI:0047', 'MI:0407'),
    'PCA': _System('MI:0090', 'MI:0407'),
    'Protein-RNA': _System('MI:0686', 'MI:0407', Gene),
    'Protein-peptide': _System('MI:0686', 'MI:0407'),
    'Proximity Label-MS': _System('MI:1313', 'MI:0915'),
    'Reconstituted Complex': _System('MI:0096', 'MI:0407'),
    'Surface Display': _System('MI:0686', 'MI:0914'),
    'Thermal Shift Assay': _System('MI:0686', 'MI:0914'),
    'Two-hybrid': _System('MI:0018', 'MI:0407'),
}
# Genetic systems: the detection method is genetic interference and each
# system has its own "sensu BioGRID" interaction type. The participants are
# the genes.
_GENETIC_SYSTEMS: dict[str, _System] = {
    'Dosage Growth Defect': _System('MI:0254', 'MI:2378', Gene),
    'Dosage Lethality': _System('MI:0254', 'MI:2377', Gene),
    'Dosage Rescue': _System('MI:0254', 'MI:2376', Gene),
    'Negative Genetic': _System('MI:0254', 'MI:2373', Gene),
    'Phenotypic Enhancement': _System('MI:0254', 'MI:2368', Gene),
    'Phenotypic Suppression': _System('MI:0254', 'MI:2374', Gene),
    'Positive Genetic': _System('MI:0254', 'MI:2371', Gene),
    'Synthetic Growth Defect': _System('MI:0254', 'MI:2369', Gene),
    'Synthetic Haploinsufficiency': _System('MI:0254', 'MI:2372', Gene),
    'Synthetic Lethality': _System('MI:0254', 'MI:2370', Gene),
    'Synthetic Rescue': _System('MI:0254', 'MI:2375', Gene),
}
_SYSTEMS = {'physical': _PHYSICAL_SYSTEMS, 'genetic': _GENETIC_SYSTEMS}
# BioGRID gives co-fractionation and co-localization the PSI-MI type
# colocalization and classifies both systems as physical interaction evidence.
_COLOCALIZATION = 'MI:0403'


# Ontology terms. The kind of a term is BioGRID's category (cell line, tissue,
# disease, ...), but the category column is unreliable: it labels 38,398 BTO and
# 1,014 CL terms of physical rows (and 188 BTO, 75 ZFA and 1 WBbt terms of
# genetic rows) "phenotype", although they are tissues, cell lines and cell
# types. The table has every combination of category and vocabulary (prefix of
# the accession) of release 5.0.262; one that is not here raises, in the parser.
# The repairs are explicit entries:
#   * an anatomy vocabulary (BTO, UBERON, EMAP, ZFA, WBbt) under "phenotype" is
#     a tissue, and a cell vocabulary (CL) a cell type;
#   * a BTO term under "phenotype" is a cell line when the file labels the same
#     accession "cell line" in any row (``_CELL_LINE_REPAIRS``; HEK-293 cell is
#     both), and otherwise a tissue. The file's own labels are the only source:
#     a BTO cell line that is never labelled "cell line" stays a tissue.
# The accession BIOGRID:0000011 ("BRENDA tissue", a placeholder of BioGRID's
# ontology, not a term of a published ontology) is carried as written under the
# kind of its category.
_KINDS: dict[tuple[str, str], str] = {
    ('biological process', 'GO'): 'biological_process',
    ('cell line', 'BTO'): 'cell_line',
    ('cell type', 'BIOGRID'): 'cell_type',
    ('cell type', 'BTO'): 'cell_type',
    ('cell type', 'CL'): 'cell_type',
    ('cellular component', 'GO'): 'cellular_component',
    ('chemical', 'CHEBI'): 'chemical',
    ('development', 'UBERON'): 'development',
    ('disease', 'DOID'): 'disease',
    ('environment', 'PATO'): 'environment',
    ('tissue', 'BTO'): 'tissue',
    ('tissue', 'EMAP'): 'tissue',
    ('tissue', 'UBERON'): 'tissue',
    # Labelled phenotype and a phenotype.
    ('phenotype', 'APO'): 'phenotype',
    ('phenotype', 'BIOGRID'): 'phenotype',
    ('phenotype', 'GO'): 'phenotype',
    ('phenotype', 'HP'): 'phenotype',
    ('phenotype', 'MP'): 'phenotype',
    ('phenotype', 'PATO'): 'phenotype',
    ('phenotype', 'WBPhenotype'): 'phenotype',
    # Labelled phenotype and an anatomy or cell term.
    ('phenotype', 'BTO'): 'tissue',
    ('phenotype', 'EMAP'): 'tissue',
    ('phenotype', 'UBERON'): 'tissue',
    ('phenotype', 'WBbt'): 'tissue',
    ('phenotype', 'ZFA'): 'tissue',
    ('phenotype', 'CL'): 'cell_type',
}
# The combinations whose terms are cell lines where the file calls the accession
# a cell line elsewhere.
_CELL_LINE_REPAIRS = frozenset({('phenotype', 'BTO')})
# The kind whose terms of a genetic interaction are associations of the two
# genes, and the namespace of each vocabulary cited (by the prefix of the
# accessions). The identifier is the accession as written (``APO:0000112``), as
# the ontology modules write the identifiers of their terms (apo, hpo and go
# keep the prefix); PATO, MP and WBPhenotype have no ``Namespace`` member and no
# ontology module, so their namespace is the lowercase prefix. A phenotype
# vocabulary that is not here raises, in the parser.
_PHENOTYPE = 'phenotype'
_PHENOTYPE_VOCABULARIES: dict[str, Namespace | str] = {
    'APO': Namespace.APO,
    'GO': Namespace.GO,
    'HP': Namespace.HPO,
    'MP': 'mp',
    'PATO': 'pato',
    'WBPhenotype': 'wbphenotype',
}
# Attributes of the interaction: the accessions of the terms of each kind and
# of their qualifiers, in the order of first appearance in the table.
_TERM_KINDS = tuple(dict.fromkeys(_KINDS.values()))


def _attribute(kind: str, qualifier: bool = False) -> str:
    return f'biogrid:{kind}_qualifier' if qualifier else f'biogrid:{kind}'


def _accessions(key: str, kind: str):
    """The accessions of the terms (or qualifiers) of one kind in the row, as
    the parser grouped them."""

    def _kind(row: dict[str, Any]) -> list[str]:
        return row.get(key, {}).get(kind, [])

    return _kind


def _phenotypes(prefix: str, suffix: str):
    """The phenotype terms of a vocabulary that a participant is linked to.
    A gene that interacts with itself is both participants and is linked
    once."""

    def _terms(row: dict[str, Any]) -> list[str]:
        if (
            suffix == 'B'
            and row.get('BioGRID ID Interactor A')
            == row.get('BioGRID ID Interactor B')
        ):
            return []
        return row.get(PHENOTYPES, {}).get(prefix, [])

    return _terms


def _term_types(row: dict[str, Any]) -> list[str]:
    """The distinct types BioGRID gives the terms of the row, as written; a
    type that the set lacks raises."""
    types = []
    for term in row[TERMS]:
        type_ = term['type']
        if type_ is None:
            continue
        if type_ not in _TERM_TYPES:
            raise ValueError(
                f'Unknown BioGRID ontology term type {type_!r} '
                f'in interaction {row.get("#BioGRID Interaction ID")}'
            )
        types.append(type_)
    return list(dict.fromkeys(types))


def _modification(row: dict[str, Any]) -> str | None:
    """The modification of the row as written; a value the table lacks raises."""
    cell = row.get('Modification')
    if cell == '-':
        return None
    if cell not in _MODIFICATIONS:
        raise ValueError(
            f'Unknown BioGRID modification {cell!r} '
            f'in interaction {row.get("#BioGRID Interaction ID")}'
        )
    return cell


def _system(row: dict[str, Any]) -> _System:
    """The system of the row; one that the table lacks raises."""
    type_ = row.get('Experimental System Type')
    name = row.get('Experimental System')
    try:
        return _SYSTEMS[type_][name]
    except KeyError:
        raise ValueError(
            f'Unknown BioGRID {type_} experimental system {name!r} '
            f'in interaction {row.get("#BioGRID Interaction ID")}'
        ) from None


def _score(row: dict[str, Any]):
    """The score of the row as a number without unit; a cell that is neither
    ``-`` nor a number raises."""
    cell = row.get('Score')
    if cell == '-':
        return None
    if not isinstance(cell, str) or not SCORE.fullmatch(cell):
        raise ValueError(
            f'Unknown BioGRID score {cell!r} '
            f'in interaction {row.get("#BioGRID Interaction ID")}'
        )
    return measurement(cell, source_field='Score')


def _publication(value: str) -> str | None:
    """``PUBMED:9006895`` or ``DOI:10.1101/...`` as a publication CURIE."""
    source, _, accession = value.partition(':')
    if source == 'PUBMED':
        return _mitab.pmid_curie(accession)
    if source == 'DOI' and accession.startswith('10.'):
        return f'doi:{accession}'
    return None


def _refseq_pairs(suffix: str):
    column = f'REFSEQ Accessions Interactor {suffix}'

    def _pairs(row: dict[str, Any]) -> list[tuple[Any, str]]:
        accessions = str(row.get(column) or '').split('|')
        return [
            (namespace, accession)
            for accession in dict.fromkeys(map(str.strip, accessions))
            if (namespace := _mitab.refseq_namespace(accession))
        ]

    return _pairs


f = FieldConfig(
    extract={
        'number': r'^(\d+)$',
        'uniprot': f'^({UNIPROT_ACCESSION})$',
    },
    map={
        'detection_method': {
            name: system.detection_method
            for systems in _SYSTEMS.values()
            for name, system in systems.items()
        },
        'interaction_type': {
            name: system.interaction_type
            for systems in _SYSTEMS.values()
            for name, system in systems.items()
        },
        'throughput': _THROUGHPUT_VALUES,
    },
    transform={
        'taxon': _mitab.taxon_curie,
        'publication': _publication,
    },
    delimiter='|',
)


def _participant(suffix: str, phenotypes: bool = False) -> EntityBuilder:
    """A gene record of BioGRID with the cross-references the row lists; with
    ``phenotypes``, associated with the phenotype terms of the row."""
    return EntityBuilder(
        entity_type=lambda row: _system(row).participant,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.ENTREZ,
                value=f(f'Entrez Gene Interactor {suffix}', extract='number'),
            ),
            CV(
                term=Namespace.BIOGRID,
                value=f(f'BioGRID ID Interactor {suffix}', extract='number'),
            ),
            CV(
                term=Namespace.UNIPROT,
                value=f(
                    f'SWISS-PROT Accessions Interactor {suffix}',
                    extract='uniprot',
                ),
            ),
            CV(
                term=Namespace.UNIPROT_TREMBL,
                value=f(
                    f'TREMBL Accessions Interactor {suffix}',
                    extract='uniprot',
                ),
            ),
            CV.from_pairs(_refseq_pairs(suffix)),
            CV(
                term=Namespace.GENESYMBOL,
                value=f(f'Official Symbol Interactor {suffix}'),
            ),
            CV(
                term=Namespace.SYNONYM,
                value=f(f'Systematic Name Interactor {suffix}'),
            ),
            CV(
                term=Namespace.SYNONYM,
                value=f(f'Synonyms Interactor {suffix}'),
            ),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=f(f'Organism ID Interactor {suffix}', transform='taxon'),
            ),
        ),
        associations=AssociationsBuilder(
            *(
                AssociationBuilder(
                    predicate=slots.associated_with,
                    object_entity_type=OntologyClass,
                    object_identifier_type=namespace,
                    object_identifier=_phenotypes(prefix, suffix),
                )
                for prefix, namespace in _PHENOTYPE_VOCABULARIES.items()
            )
        )
        if phenotypes
        else None,
    )


def _physical_predicate(row: dict[str, Any]):
    """
    Predicate for the PSI-MI interaction type BioGRID gives the system.

    Direct interaction is a direct physical interaction; association and
    physical association say only that the participants interact, and so
    does colocalization as BioGRID uses it.
    """
    interaction_type = _system(row).interaction_type
    if interaction_type == _COLOCALIZATION:
        return slots.interacts_with
    return _mitab.physical_predicate([interaction_type])


def _schema(predicate, phenotypes: bool = False) -> RelationBuilder:
    return RelationBuilder(
        subject=_participant('A', phenotypes),
        predicate=predicate,
        object=_participant('B', phenotypes),
        identifiers=IdentifiersBuilder(
            CV(
                term=_BIOGRID_INTERACTION,
                value=f('#BioGRID Interaction ID', extract='number'),
            ),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.has_evidence_of_type,
                value=f('Experimental System', map='detection_method'),
            ),
            CV(
                term=slots.original_predicate,
                value=f('Experimental System', map='interaction_type'),
            ),
            # Several systems share a PSI-MI method, so the system is kept.
            CV(term=_EXPERIMENTAL_SYSTEM, value=f('Experimental System')),
            # An interaction of both kinds of study has both values.
            CV(term=_THROUGHPUT, value=f('Throughput', map='throughput')),
            CV(
                term=slots.publications,
                value=f('Publication Source', transform='publication'),
            ),
            CV(term=_SCORE, value=_score),
            # The experimental context: the modification observed, the
            # accessions of the ontology terms and of their qualifiers by
            # kind and the type of the phenotype terms. The phenotype terms
            # of a genetic interaction are associations of the participants.
            CV(term=_MODIFICATION, value=_modification),
            *(
                CV(term=_attribute(kind), value=_accessions(TERMS_BY_KIND, kind))
                for kind in _TERM_KINDS
            ),
            *(
                CV(
                    term=_attribute(kind, qualifier=True),
                    value=_accessions(QUALIFIERS_BY_KIND, kind),
                )
                for kind in _TERM_KINDS
            ),
            CV(term=_PHENOTYPE_TYPE, value=_term_types),
        ),
    )


physical_schema = _schema(_physical_predicate)
genetic_schema = _schema(slots.genetically_interacts_with, phenotypes=True)


def _interaction(schema: RelationBuilder):
    def _map(row: dict[str, Any]):
        # The builders drop what they cannot read, so a system, a score, a
        # modification or a term of an unknown form is checked here, where it
        # raises.
        _system(row)
        _score(row)
        _modification(row)
        _term_types(row)
        return schema(row)

    return _map


map_physical_interaction = _interaction(physical_schema)
map_genetic_interaction = _interaction(genetic_schema)


def _dataset(system_type: str, mapper) -> Dataset:
    return Dataset(
        download=Download(
            url=_URL,
            filename=_FILENAME,
            subfolder='biogrid',
            large=True,
            ext='zip',
        ),
        mapper=mapper,
        raw_parser=functools.partial(
            iter_interactions,
            system_type=system_type,
            release=_RELEASE,
            systems={kind: set(names) for kind, names in _SYSTEMS.items()},
            term_types=_TERM_TYPES,
            modifications=_MODIFICATIONS,
            term_kinds=_KINDS,
            cell_line_repairs=_CELL_LINE_REPAIRS,
            phenotype_kind=_PHENOTYPE,
            phenotype_vocabularies=frozenset(_PHENOTYPE_VOCABULARIES),
        ),
    )


resource = Resource(
    config,
    physical=_dataset('physical', map_physical_interaction),
    genetic=_dataset('genetic', map_genetic_interaction),
)
