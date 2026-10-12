"""
Parse ELM motif-domain interactions and emit Relation records.

ELM, the Eukaryotic Linear Motif resource, curates short linear motifs: motif
classes and their experimentally validated instances in proteins. Its
interaction table lists, for a motif instance, the protein whose domain
binds it: the motif class and the residues of the instance, the domain
(a Pfam, InterPro or SMART accession) and where it lies in the partner, a
measured affinity where one was reported, and the publications.

This module emits the interaction between the two proteins, with the motif
instance and the domain as binding regions of the participants. It reads two
exports of the same records: the interaction table, and the PSI-MI TAB
export for the record numbers and the detection methods that the table
lacks, and checks the pairing of the two against the PSI-MI XML export of the
records, which names the motif and the domain of every record number. In that
XML every interaction has the type "direct interaction" (MI:0407) and every
motif and domain the feature type "direct binding region" (MI:1125), which is
where the predicate and the region type used here come from.

Motif classes and instances without a binding partner, and the other ELM
tables (switches, pathways, diseases), are not read.
"""

from __future__ import annotations

from collections.abc import Callable
import functools
import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2 import _mitab
from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.elm import (
    iter_interactions,
    parse_position,
    split_methods,
)


config = ResourceConfig(
    id=ResourceCv.ELM,
    name='ELM',
    full='Eukaryotic Linear Motif resource',
    url='http://elm.eu.org/',
    # ELM Software License Agreement (EMBLEM): academic research only; anyone else
    # needs a commercial license. http://elm.eu.org/media/Elm_academic_license.pdf
    license=LicenseCV.NONCOMMERCIAL,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='37962385',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Interactions between linear motif instances and the protein '
        'domains that bind them, with detection methods and publications.'
    ),
    description=(
        'The Eukaryotic Linear Motif (ELM) resource is a manually curated '
        'database of short linear motifs: compact, degenerate sequence '
        'modules, mostly in disordered regions of proteins, that mediate '
        'docking, targeting, degradation, cleavage and post-translational '
        'modification. ELM defines motif classes by a sequence pattern and '
        'annotates experimentally validated instances of them in proteins. '
        'For many instances it also records the interaction partner, the '
        'globular domain that recognises the motif, and the measured binding '
        'affinity.'
    ),
)

# The ELM server builds an export when it is requested and sends nothing until
# it is ready, which takes about two minutes.
_SLOW_SERVER = {'timeout': 600}
download_interactions = Download(
    url='http://elm.eu.org/interactions/as_tsv',
    filename='elm_interactions.tsv',
    subfolder='elm',
    large=True,
    download_kwargs=_SLOW_SERVER,
)
# The same records, in the same order, as PSI-MI TAB 2.5.
download_evidence = Download(
    url='http://elm.eu.org/interactions/as_mitab',
    filename='elm_interactions.mitab',
    subfolder='elm',
    large=True,
    download_kwargs=_SLOW_SERVER,
)

# The same records, in PSI-MI XML: by record number, with the motif and the
# domain of each. It checks that the pairing of the two exports by position
# gave every table row its own record.
download_features = Download(
    url='http://elm.eu.org/interactions/as_psimi',
    filename='elm_interactions.psimi.xml',
    subfolder='elm',
    large=True,
    download_kwargs=_SLOW_SERVER,
    encoding=None,
    default_mode='rb',
)

# Identifier systems of the domain column, by the syntax of the accession.
_DOMAIN_NAMESPACES = (
    (re.compile(r'PF\d{5}'), 'pfam'),
    (re.compile(r'IPR\d{6}'), 'interpro'),
    (re.compile(r'SM\d{5}'), 'smart'),
)
# An ELM class identifier starts with the type of the class: cleavage site,
# degron, docking site, ligand binding site, modification site, targeting site.
_CLASS_RE = re.compile(r'(?:CLV|DEG|DOC|LIG|MOD|TRG)_[\w-]+')
_ELM_CLASS = 'elm'
# ELM's numbers of its interaction records (``ELMINT000120``).
_ELM_INTERACTION = 'elm_interaction'
_INTERACTION_RE = re.compile(r'ELMINT\d+')
# A UniProt accession, with the isoform number some records give.
_ACCESSION_RE = re.compile(
    r'^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?)$'
)
# The feature type of motifs and domains in ELM's PSI-MI XML export (MI:1125).
_BINDING_REGION = 'direct binding region'


def _region(
    row: dict[str, Any],
    accession_column: str,
    start_column: str,
    stop_column: str,
    identifier: dict[str, str] | None,
) -> dict[str, Any] | None:
    start, stop = parse_position(row.get(start_column)), parse_position(row.get(stop_column))
    if start is None or stop is None or stop < start:
        start = stop = None
    if identifier is None and start is None:
        return None
    # The sequence is named only when the participant is an isoform: an entry
    # accession does not say which of its sequences the positions count.
    accession = str(row.get(accession_column) or '').strip()
    isoform = _ACCESSION_RE.match(accession) and '-' in accession
    return {
        'regions': [
            {
                'type': _BINDING_REGION,
                'identifier': identifier,
                'position': start,
                'end_position': stop,
                'coordinate_reference': {
                    'identifier': (
                        {'ns': 'uniprot', 'id': accession} if isoform else None
                    ),
                    'coordinate_system': 'protein',
                    'position_base': 1,
                }
                if start is not None
                else None,
            }
        ],
    }


def _motif_form(row: dict[str, Any]) -> dict[str, Any] | None:
    """The motif instance: its class and its residues in the protein."""
    motif_class = str(row.get('Elm') or '').strip()
    return _region(
        row,
        'interactorElm',
        'StartElm',
        'StopElm',
        {'ns': _ELM_CLASS, 'id': motif_class}
        if _CLASS_RE.fullmatch(motif_class)
        else None,
    )


def _domain_form(row: dict[str, Any]) -> dict[str, Any] | None:
    """The domain that binds the motif, and its residues where ELM gives them."""
    domain = str(row.get('Domain') or '').strip()
    namespace = next(
        (ns for pattern, ns in _DOMAIN_NAMESPACES if pattern.fullmatch(domain)),
        None,
    )
    return _region(
        row,
        'interactorDomain',
        'StartDomain',
        'StopDomain',
        {'ns': namespace, 'id': domain} if namespace else None,
    )


def _gene_name(column: str) -> Callable[[dict[str, Any]], list[str]]:
    """The gene name of the alias column (``uniprotkb:MCD1(gene name)``)."""
    return lambda row: [
        item.accession
        for item in _mitab.parse_mitab_cell(row.get(column))
        if item.database == 'uniprotkb' and item.description == 'gene name'
    ]


def _methods(row: dict[str, Any]) -> list[str]:
    return list(
        dict.fromkeys(
            term
            for method in split_methods(row.get('Methods'))
            for term in _mitab.mi_terms(method)
        )
    )


def _interaction_record(row: dict[str, Any]) -> str | None:
    record = str(row.get('Interaction_identifier') or '').strip()
    return record if _INTERACTION_RE.fullmatch(record) else None


f = FieldConfig(
    extract={
        'accession': _ACCESSION_RE,
        # ``"9606"(Homo sapiens)``; the table reader removes the quotes.
        'taxon': r'^"?(\d+)"?\(',
        'pmid': PMID_EXTRACT,
    },
    transform={
        'taxon': lambda v: f'NCBITaxon:{v}',
        'pmid': lambda v: f'PMID:{v}',
    },
)


def _protein(
    role: str, aliases: str, form: Callable[[dict[str, Any]], Any]
) -> Callable[[dict[str, Any]], Any]:
    builder = EntityBuilder(
        entity_type=Protein,
        molecular_form=form,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=f(f'interactor{role}', extract='accession'),
            ),
            CV(term=Namespace.GENESYMBOL, value=_gene_name(aliases)),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=f(f'taxonomy{role}', extract='taxon', transform='taxon'),
            ),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        accession = str(row.get(f'interactor{role}') or '').strip()
        return builder.build(row) if _ACCESSION_RE.match(accession) else None

    return build


# The protein with the motif is the subject and the protein with the domain
# the object. ELM types every one of these records a direct interaction.
interactions_schema = RelationBuilder(
    subject=_protein('Elm', 'Aliases for A', _motif_form),
    predicate=slots.directly_physically_interacts_with,
    object=_protein('Domain', 'Aliases for B', _domain_form),
    identifiers=IdentifiersBuilder(
        CV(term=_ELM_INTERACTION, value=_interaction_record),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_evidence_of_type, value=_methods),
        CV(
            term=slots.publications,
            value=f('PMID', delimiter=',', extract='pmid', transform='pmid'),
        ),
    ),
)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=functools.partial(
            iter_interactions,
            evidence=download_evidence,
            features=download_features,
        ),
    ),
)
