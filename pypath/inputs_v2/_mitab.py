"""
Shared handling of PSI-MI TAB (MITAB) fields.

MITAB is the tab-separated interaction format of the HUPO Proteomics
Standards Initiative (versions 2.5 to 2.8, each adding columns to the
previous one):
https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md

A cell holds items separated by ``|``. An item is written
``database:accession(description)``. The accession and the description are
quoted when they contain a reserved character, as in
``psi-mi:"MI:0018"(two hybrid)``, and ``-`` stands for no value.

This module has the parts that do not depend on one source: reading a cell,
the identifier systems behind the PSI-MI database names, PubMed, taxonomy and
PSI-MI term extraction, the PSI-MI interaction types that state a direct
interaction, and annotation definitions for the standard columns. Column
names differ between sources, so every helper takes the column name. What a
source writes in its own way stays in that source's module.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Generator, Iterable, Mapping
from typing import Any, NamedTuple

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2._publications import pmid_curie
from pypath.internals.tabular_builder import CV


class MitabItem(NamedTuple):
    """One ``database:accession(description)`` item of a MITAB cell."""

    # Lowercase database name; None for an item written without one.
    database: str | None
    accession: str
    description: str | None


def _split_items(cell: str) -> list[str]:
    """Split a cell at the ``|`` that are outside double quotes."""
    # An odd number of quotes is not valid MITAB; take every ``|`` literally.
    if '"' not in cell or cell.count('"') % 2:
        return cell.split('|')
    items, start, quoted, position = [], 0, False, 0
    while position < len(cell):
        character = cell[position]
        if character == '\\' and quoted:
            position += 1
        elif character == '"':
            quoted = not quoted
        elif character == '|' and not quoted:
            items.append(cell[start:position])
            start = position + 1
        position += 1
    items.append(cell[start:])
    return items


def _closing_quote(text: str) -> int:
    """Position of the quote that closes the one at position 0, or -1."""
    position = 1
    while position < len(text):
        if text[position] == '\\':
            position += 1
        elif text[position] == '"':
            return position
        position += 1
    return -1


def _unquote(text: str) -> str:
    text = text.strip()
    if len(text) >= 2 and text[0] == '"' and _closing_quote(text) == len(text) - 1:
        text = text[1:-1]
    return text.replace('\\"', '"').strip()


def _opening_parenthesis(text: str) -> int:
    """Position of the parenthesis that the last character closes, or -1."""
    depth = 0
    for position in range(len(text) - 1, -1, -1):
        if text[position] == ')':
            depth += 1
        elif text[position] == '(':
            depth -= 1
            if depth == 0:
                return position
    return -1


def parse_mitab_item(text: object) -> MitabItem | None:
    """
    Read one item of a MITAB cell.

    The database name ends at the first colon. The description is the
    parenthesis that ends the item, so an accession may itself contain
    parentheses (``p38(MAPK)(gene name synonym)``). An item without a colon,
    such as the bare ``DIP-617N`` of DIP files, has no database.
    """
    text = str(text if text is not None else '').strip()
    if not text or text == '-':
        return None
    database, rest = None, text
    if text[0] == '"':
        closing = _closing_quote(text)
        if closing > 0 and text[closing + 1 : closing + 2] == ':':
            database, rest = text[1:closing], text[closing + 2 :]
    else:
        match = re.match(r'([^:"(]+):', text)
        if match:
            database, rest = match.group(1), text[match.end() :]
    rest = rest.strip()
    accession, description = rest, None
    if rest.startswith('"'):
        closing = _closing_quote(rest)
        tail = rest[closing + 1 :].strip() if closing > 0 else ''
        if closing > 0 and (not tail or (tail[0] == '(' and tail[-1] == ')')):
            accession = rest[: closing + 1]
            description = tail[1:-1] if tail else None
    elif rest.endswith(')'):
        opening = _opening_parenthesis(rest)
        if opening > 0:
            accession, description = rest[:opening], rest[opening + 1 : -1]
    accession = _unquote(accession)
    if not accession or accession == '-':
        return None
    return MitabItem(
        database=database.strip().lower() if database else None,
        accession=accession,
        description=(_unquote(description) or None) if description else None,
    )


def parse_mitab_cell(cell: object) -> list[MitabItem]:
    """The items of a MITAB cell, in file order; ``-`` and empty give none."""
    text = str(cell if cell is not None else '').strip()
    if not text or text == '-':
        return []
    return [
        item
        for item in map(parse_mitab_item, _split_items(text))
        if item is not None
    ]


def _unique(values: Iterable[Any]) -> list[Any]:
    return list(dict.fromkeys(value for value in values if value is not None))


def checked_rows(
    rows: Iterable[dict[str, Any]],
    columns: Iterable[str],
    source: str,
) -> Generator[dict[str, Any], None, None]:
    """
    Pass on the rows of a table that has a header line, checking its layout.

    Every row must have the named columns and no field beyond the header.
    Anything else is an error: it means the download is not the table, or
    that its layout changed, and must not end as an empty dataset.
    """
    columns = tuple(columns)
    for number, row in enumerate(rows, start=2):
        missing = [column for column in columns if row.get(column) is None]
        if missing or None in row:
            problem = f'missing {missing}' if missing else 'fields beyond the header'
            raise ValueError(
                f'{source} line {number}: expected {len(columns)} named '
                f'columns, {problem}'
            )
        yield row


# PSI-MI terms ------------------------------------------------------------

_MI_RE = re.compile(r'(?:MI:)?(\d{4})')


def mi_term(item: MitabItem) -> str | None:
    """
    The PSI-MI accession an item names, as ``MI:0018``.

    Sources write ``psi-mi:"MI:0018"(two hybrid)``, ``psi-mi:MI:0018(...)``
    or ``MI:0018(...)``. An item of another database, or a name in place of
    an accession, gives None.
    """
    if item.database not in {'psi-mi', 'mi'}:
        return None
    match = _MI_RE.fullmatch(item.accession)
    return f'MI:{match.group(1)}' if match else None


def mi_terms(cell: object) -> list[str]:
    """The distinct PSI-MI accessions of a cell, in file order."""
    return _unique(map(mi_term, parse_mitab_cell(cell)))


# Publications and taxonomy -----------------------------------------------


def pmid_curies(cell: object) -> list[str]:
    """
    The distinct PubMed identifiers of a publication cell.

    Other entries of the cell contribute nothing: IMEx identifiers, DOIs, and
    the placeholders some sources write under the ``pubmed`` database name
    (``pubmed:DIP-17838S``, ``pubmed:unassigned1304``).
    """
    return _unique(
        pmid_curie(item.accession)
        for item in parse_mitab_cell(cell)
        if item.database == 'pubmed'
    )


def taxon_curie(value: object) -> str | None:
    """
    ``NCBITaxon:9606`` for an NCBI Taxonomy identifier.

    MITAB also uses the taxid field for what is not an organism: -1 in vitro,
    -2 chemical synthesis, -3 unknown, -4 in vivo, -5 in silico, and some
    sources write 0. These give None.
    """
    text = str(value if value is not None else '').strip()
    return f'NCBITaxon:{int(text)}' if text.isdigit() and int(text) > 0 else None


def taxon_curies(cell: object) -> list[str]:
    """
    The distinct organisms of a taxid cell.

    A cell often names one organism twice, with a short and a full label
    (``taxid:9606(human)|taxid:9606(Homo sapiens)``).
    """
    return _unique(
        taxon_curie(item.accession)
        for item in parse_mitab_cell(cell)
        if item.database == 'taxid'
    )


# Identifiers -------------------------------------------------------------

# PSI-MI database names, as sources write them, and their identifier system.
# A database that is absent has no registered namespace: its items are left
# out, and a source that needs them maps them itself.
MITAB_NAMESPACES: dict[str, Namespace] = {
    'bind smid': Namespace.BIND,
    'biogrid': Namespace.BIOGRID,
    'cas registry number': Namespace.CAS,
    'chebi': Namespace.CHEBI,
    'chembl': Namespace.CHEMBL,
    'chembl compound': Namespace.CHEMBL,
    'complex portal': Namespace.COMPLEXPORTAL,
    'ddbj/embl/genbank': Namespace.GENBANK,
    'dip': Namespace.DIP,
    'ensembl': Namespace.ENSEMBL,
    'ensemblgenomes': Namespace.ENSEMBL_GENOMES,
    'entrez gene/locuslink': Namespace.ENTREZ,
    'entrezgene/locuslink': Namespace.ENTREZ,
    'flybase': Namespace.FLYBASE,
    'genbank identifier': Namespace.GENBANK,
    'genbank_nucl_gi': Namespace.GENBANK_GI,
    'genbank_protein_gi': Namespace.GENBANK_GI,
    'imex': Namespace.IMEX,
    'intact': Namespace.INTACT,
    'ipi': Namespace.IPI,
    'mint': Namespace.MINT,
    'mirbase': Namespace.MIRBASE,
    'refseq': Namespace.REFSEQ,
    'rfam': Namespace.RFAM,
    'rnacentral': Namespace.RNACENTRAL,
    'uniparc': Namespace.UNIPARC,
    'uniprot/swiss-prot': Namespace.UNIPROT,
    'uniprot/trembl': Namespace.UNIPROT_TREMBL,
    'uniprotkb': Namespace.UNIPROT,
}

_ENSEMBL_RE = re.compile(r'ENS[A-Z]*([GTP])\d+(?:\.\d+)?')
_ENSEMBL_NAMESPACES = {
    'G': Namespace.ENSG,
    'T': Namespace.ENST,
    'P': Namespace.ENSP,
}
_REFSEQ_RE = re.compile(r'([A-Z]{2})_\d+(?:\.\d+)?')
_REFSEQ_PROTEIN_PREFIXES = frozenset({'AP', 'NP', 'XP', 'YP', 'WP', 'ZP'})
_CHEBI_RE = re.compile(r'(?:CHEBI:)?(\d+)', re.IGNORECASE)


def ensembl_namespace(accession: str) -> Namespace | None:
    """Gene, transcript or protein namespace of an Ensembl stable identifier."""
    match = _ENSEMBL_RE.fullmatch(accession)
    return _ENSEMBL_NAMESPACES[match.group(1)] if match else None


def refseq_namespace(accession: str) -> Namespace | None:
    """Protein or nucleotide namespace of a RefSeq accession."""
    match = _REFSEQ_RE.fullmatch(accession)
    if not match:
        return None
    if match.group(1) in _REFSEQ_PROTEIN_PREFIXES:
        return Namespace.REFSEQ_PROTEIN
    return Namespace.REFSEQ


def identifier_pair(
    item: MitabItem,
    namespaces: Mapping[str, Any] = MITAB_NAMESPACES,
) -> tuple[Any, str] | None:
    """
    The namespace and accession of an identifier item.

    The database name selects the namespace. Where the identifier system has
    a fixed accession syntax, an accession that does not follow it gives
    None: sources also put names under these database names, such as the
    gene symbol in ``entrez gene/locuslink:MAP2K4``.
    """
    namespace = namespaces.get(item.database or '')
    if namespace is None:
        return None
    accession = item.accession
    if namespace in (Namespace.ENTREZ, Namespace.GENBANK_GI):
        return (namespace, accession) if accession.isdigit() else None
    if namespace == Namespace.ENSEMBL:
        specific = ensembl_namespace(accession)
        return (specific, accession) if specific else None
    if namespace == Namespace.REFSEQ:
        specific = refseq_namespace(accession)
        return (specific, accession) if specific else None
    if namespace == Namespace.CHEBI:
        match = _CHEBI_RE.fullmatch(accession)
        return (namespace, match.group(1)) if match else None
    return namespace, accession


def identifier_pairs(
    cell: object,
    namespaces: Mapping[str, Any] = MITAB_NAMESPACES,
) -> list[tuple[Any, str]]:
    """The distinct ``(namespace, accession)`` pairs of an identifier cell."""
    return _unique(
        identifier_pair(item, namespaces) for item in parse_mitab_cell(cell)
    )


# Interactor and interaction types ----------------------------------------

# PSI-MI interactor type -> Biolink class.
INTERACTOR_TYPES: dict[str, type] = {
    'MI:0326': model.Protein,  # protein
    'MI:0327': model.Protein,  # peptide
    'MI:0314': model.MacromolecularComplex,  # complex
    'MI:0328': model.ChemicalEntity,  # small molecule
    'MI:0320': model.RNAProduct,  # ribonucleic acid
    'MI:0250': model.Gene,  # gene
    'MI:0318': model.NucleicAcidEntity,  # nucleic acid
    'MI:0319': model.NucleicAcidEntity,  # deoxyribonucleic acid
    'MI:0681': model.NucleicAcidEntity,  # double stranded dna
    'MI:0317': model.MolecularEntity,  # interactor (unspecified)
}


def interactor_type(
    cell: object,
    types: Mapping[str, type] = INTERACTOR_TYPES,
) -> type | None:
    """The Biolink class of an interactor type cell; None for an unknown type."""
    terms = mi_terms(cell)
    return types.get(terms[0]) if terms else None


# Interaction types that leave open whether the participants touch.
ASSOCIATION_TYPES = frozenset(
    {
        'MI:0914',  # association
        'MI:0915',  # physical association
    }
)
# Direct interaction (MI:0407) and every type below it in the PSI-MI ontology
# (psi-mi.obo of 2021-04-15): covalent binding, enzymatic reactions and their
# kinds (phosphorylation, cleavage, ubiquitination...), self interaction.
DIRECT_INTERACTION_TYPES = frozenset(
    f'MI:{number}'
    for number in (
        '0407 '
        '0192 0193 0194 0195 0197 0198 0199 0200 0201 0202 0203 0204 0206 '
        '0207 0209 0210 0211 0212 0213 0214 0216 0217 0220 0408 0414 0556 '
        '0557 0558 0559 0566 0567 0568 0569 0570 0571 0572 0701 0844 0871 '
        '0881 0882 0883 0902 0910 0945 0971 0985 0986 0987 1027 1126 1127 '
        '1139 1140 1143 1146 1148 1230 1237 1250 1251 1310 1327 1355 2252 '
        '2263 2272 2273 2280 2356'
    ).split()
)


def physical_predicate(interaction_types: Iterable[str]):
    """
    The predicate that the PSI-MI interaction types of a record support.

    ``directly_physically_interacts_with`` when a type states a direct
    interaction, ``interacts_with`` for an association or physical
    association. Any other type gives None: colocalization, genetic
    interaction and predicted interaction are different assertions, and an
    accession that is not an interaction type states none. The source module
    decides what to do with those.
    """
    interaction_types = set(interaction_types)
    if interaction_types & DIRECT_INTERACTION_TYPES:
        return slots.directly_physically_interacts_with
    if interaction_types & ASSOCIATION_TYPES:
        return slots.interacts_with
    return None


def is_negative(value: object) -> bool:
    """True when the ``negative`` column (MITAB 2.6 and later) says so."""
    return str(value if value is not None else '').strip().lower() == 'true'


# Annotation definitions --------------------------------------------------

# PSI-MI participant attributes without a Biolink slot; values are MI terms.
BIOLOGICAL_ROLE = 'psi_mi:biological_role'
EXPERIMENTAL_ROLE = 'psi_mi:experimental_role'
IDENTIFICATION_METHOD = 'psi_mi:participant_identification_method'
# "unspecified role" says that no role was recorded.
UNSPECIFIED_ROLE = 'MI:0499'


def mi_terms_of(column: str) -> Callable[[Mapping[str, Any]], list[str]]:
    """Row function: the PSI-MI accessions in ``column``."""
    return lambda row: mi_terms(row.get(column))


def roles_of(column: str) -> Callable[[Mapping[str, Any]], list[str]]:
    """Row function: the recorded roles in ``column``."""
    return lambda row: [
        term for term in mi_terms(row.get(column)) if term != UNSPECIFIED_ROLE
    ]


def pmids_of(column: str) -> Callable[[Mapping[str, Any]], list[str]]:
    """Row function: the PubMed CURIEs in ``column``."""
    return lambda row: pmid_curies(row.get(column))


def taxa_of(column: str) -> Callable[[Mapping[str, Any]], list[str]]:
    """Row function: the NCBI Taxonomy CURIEs in ``column``."""
    return lambda row: taxon_curies(row.get(column))


def identifiers_of(
    *columns: str,
    namespaces: Mapping[str, Any] = MITAB_NAMESPACES,
) -> Callable[[Mapping[str, Any]], list[tuple[Any, str]]]:
    """Row function: the identifier pairs of the columns, for ``CV.from_pairs``."""
    return lambda row: _unique(
        pair
        for column in columns
        for pair in identifier_pairs(row.get(column), namespaces)
    )


def interaction_annotations(
    *,
    detection_method: str | None = None,
    interaction_type: str | None = None,
    source_database: str | None = None,
    publications: str | None = None,
    host_organism: str | None = None,
) -> list[CV]:
    """
    Annotations of an interaction from the standard columns a source has.

    Each argument is a column name: the detection methods become evidence
    types, the interaction types the original predicate and the source
    databases the supporting data sources, all as PSI-MI accessions; PubMed
    identifiers become publications and the host organism the taxon the
    interaction was observed in.
    """
    annotations = []
    if detection_method:
        annotations.append(
            CV(
                term=slots.has_evidence_of_type,
                value=mi_terms_of(detection_method),
            )
        )
    if interaction_type:
        annotations.append(
            CV(
                term=slots.original_predicate,
                value=mi_terms_of(interaction_type),
            )
        )
    if source_database:
        annotations.append(
            CV(
                term=slots.supporting_data_source,
                value=mi_terms_of(source_database),
            )
        )
    if publications:
        annotations.append(
            CV(term=slots.publications, value=pmids_of(publications))
        )
    if host_organism:
        annotations.append(
            CV(term=slots.in_taxon, value=taxa_of(host_organism))
        )
    return annotations


def participant_annotations(
    *,
    taxon: str | None = None,
    biological_role: str | None = None,
    experimental_role: str | None = None,
    identification_method: str | None = None,
) -> list[CV]:
    """
    Annotations of a participant from the standard columns a source has.

    Each argument is a column name. Roles and identification methods are
    PSI-MI accessions; an unspecified role is left out.
    """
    annotations = []
    if taxon:
        annotations.append(CV(term=slots.in_taxon, value=taxa_of(taxon)))
    if biological_role:
        annotations.append(
            CV(term=BIOLOGICAL_ROLE, value=roles_of(biological_role))
        )
    if experimental_role:
        annotations.append(
            CV(term=EXPERIMENTAL_ROLE, value=roles_of(experimental_role))
        )
    if identification_method:
        annotations.append(
            CV(
                term=IDENTIFICATION_METHOD,
                value=mi_terms_of(identification_method),
            )
        )
    return annotations
