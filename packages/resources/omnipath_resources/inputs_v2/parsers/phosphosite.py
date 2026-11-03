"""Read the PhosphoSitePlus kinase-substrate table and its BioPAX evidence file.

Source file (``Kinase_Substrate_Dataset``, tab separated). Three lines precede
the column names: the release date, the license statement and an empty line.

    050522
    Data extracted from PhosphoSitePlus(R), created by Cell Signaling ...

    GENE  KINASE  KIN_ACC_ID  KIN_ORGANISM  SUBSTRATE  SUB_GENE_ID  SUB_ACC_ID ...
    Pak2  PAK2    Q64303      rat           MEK1       170851       Q01986     ...

The second file, ``Kinase_substrates.owl`` (BioPAX level 3, RDF/XML), lists the
same kinase-substrate sites with the publications and methods behind them. One
``Catalysis`` has a controller (the kinase, a ``Protein`` with a
``ProteinReference`` that carries the accession), and a controlled
``BiochemicalReaction`` whose right-hand protein has one ``ModificationFeature``
(``ModificationFeature_S336p``: the residue, its position and the PSI-MOD
term). The reaction is the phosphorylation of the substrate site, so the
``Evidence`` elements on it are shared by every kinase of the site. An
``Evidence`` has a ``PublicationXref`` (a PubMed record or a PhosphoSitePlus
record without PubMed ID) and, for most, an ``evidenceCode`` with a PSI-MI
term. The file defines each object once and refers to it by ``rdf:resource``
elsewhere, so the elements have to be resolved by identifier.

The two files are joined on kinase accession, substrate accession, residue and
position (``iter_kinase_substrate_evidence``).
"""

from __future__ import annotations

from collections.abc import Collection, Generator, Iterable, Mapping
import itertools
import re
from types import SimpleNamespace
from typing import Any
import xml.etree.ElementTree as ET

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv, read_xml_root

# Columns of the table, in source order.
COLUMNS = (
    'GENE',
    'KINASE',
    'KIN_ACC_ID',
    'KIN_ORGANISM',
    'SUBSTRATE',
    'SUB_GENE_ID',
    'SUB_ACC_ID',
    'SUB_GENE',
    'SUB_ORGANISM',
    'SUB_MOD_RSD',
    'SITE_GRP_ID',
    'SITE_+/-7_AA',
    'DOMAIN',
    'IN_VIVO_RXN',
    'IN_VITRO_RXN',
    'CST_CAT#',
)
# The license statement on the second line of the file; the license term of the
# module rests on it.
_LICENSE_STATEMENT = (
    'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License'
)
# "S336": one-letter residue and its position in the substrate's sequence.
_SITE_REGEX = re.compile(r'^([A-Z])([1-9]\d*)$')
# The residue with seven residues on each side; modified residues are in lower
# case and "_" pads a window that runs past an end of the sequence.
_WINDOW_REGEX = re.compile(r'^[A-Za-z_]{15}$')
# The two reaction columns hold "X" or a blank.
_FLAGS = {'X': True, '': False}


def iter_kinase_substrate(
    opener, *, organisms: Collection[str], **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per row of the kinase-substrate table.

    The lines above the column names are read here and the table itself by the
    shared TSV reader. A file without the license statement or with other
    columns is an error, not an empty dataset. An organism label that is not
    in ``organisms``, a site, a sequence window or a reaction flag of another
    form is an error as well.

    Each dictionary has the source columns as written, and:

    - ``release``: the first line of the file, the release date (``050522``).
    - ``residue``, ``position``: one-letter residue and integer position read
      from ``SUB_MOD_RSD``.
    - ``in_vivo``, ``in_vitro``: the two reaction columns as booleans.
    """
    if not opener or not opener.result:
        return
    lines = iter(opener.result)
    preamble = []
    header = None
    for line in itertools.islice(lines, 10):
        if line.startswith('GENE\t'):
            header = line
            break
        preamble.append(line.strip())
    if header is None:
        raise ValueError('PhosphoSite table: no line with the column names')
    if header.rstrip('\r\n').split('\t') != list(COLUMNS):
        raise ValueError(f'PhosphoSite table has other columns: {header!r}')
    if not any(_LICENSE_STATEMENT in line for line in preamble):
        raise ValueError(
            'PhosphoSite table lacks the CC BY-NC-SA 3.0 statement; '
            'check the terms of this file before using it'
        )
    release = preamble[0] if preamble else None
    table = SimpleNamespace(result=itertools.chain([header], lines))
    count = 0
    for row in iter_tsv(table):
        site = _SITE_REGEX.match(row['SUB_MOD_RSD'] or '')
        if not site:
            raise ValueError(f'PhosphoSite site not understood: {row["SUB_MOD_RSD"]!r}')
        window = row['SITE_+/-7_AA'] or ''
        if not _WINDOW_REGEX.match(window) or window[7].upper() != site[1]:
            raise ValueError(
                f'PhosphoSite sequence window does not fit {site[0]}: {window!r}'
            )
        for column in ('KIN_ORGANISM', 'SUB_ORGANISM'):
            if row[column] not in organisms:
                raise ValueError(f'Unknown PhosphoSite organism: {row[column]!r}')
        row['release'] = release
        row['residue'] = site[1]
        row['position'] = int(site[2])
        for key, column in (('in_vivo', 'IN_VIVO_RXN'), ('in_vitro', 'IN_VITRO_RXN')):
            flag = (row[column] or '').strip()
            if flag not in _FLAGS:
                raise ValueError(f'Unknown PhosphoSite {column} value: {row[column]!r}')
            row[key] = _FLAGS[flag]
        count += 1
        yield row
    if not count:
        raise ValueError('PhosphoSite table has no rows')


# ================================ BioPAX evidence ============================

_BP = '{http://www.biopax.org/release/biopax-level3.owl#}'
_RDF = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
# Phosphorylated residues of the file: PSI-MOD term, its label, the residue.
_PHOSPHO_RESIDUES = {
    'MOD:00046': ('o-phospho-l-serine', 'S'),
    'MOD:00047': ('o-phospho-l-threonine', 'T'),
    'MOD:00048': ("o4'-phospho-l-tyrosine", 'Y'),
}
# "ModificationFeature_S336p": residue and position of the phosphorylated site.
_FEATURE_REGEX = re.compile(r'^ModificationFeature_([A-Z])([1-9]\d*)p$')
# The fraction of the file's catalysis records that must find a row of the
# table. The two files are one release, and 99.7% match; a release mismatch
# would leave far fewer.
MIN_MATCHED_FRACTION = 0.95


def _resolve(prop: ET.Element, index: Mapping[str, ET.Element]) -> ET.Element:
    """The object a property element holds: its child, or what it refers to."""
    reference = prop.get(_RDF + 'resource')
    if reference is not None:
        if not reference.startswith('#') or reference[1:] not in index:
            raise ValueError(f'PhosphoSite BioPAX: unresolved reference {reference!r}')
        return index[reference[1:]]
    children = list(prop)
    if len(children) != 1:
        raise ValueError(
            f'PhosphoSite BioPAX: <{prop.tag}> holds {len(children)} objects, expected 1'
        )
    return children[0]


def _one_prop(element: ET.Element, tag: str) -> ET.Element:
    """The only ``tag`` property of ``element``."""
    props = element.findall(_BP + tag)
    if len(props) != 1:
        raise ValueError(
            f'PhosphoSite BioPAX: {element.get(_RDF + "ID")!r} has {len(props)} '
            f'<{tag}> properties, expected 1'
        )
    return props[0]


def _one(element: ET.Element, tag: str, index: Mapping[str, ET.Element]) -> ET.Element:
    """The object of the only ``tag`` property of ``element``."""
    return _resolve(_one_prop(element, tag), index)


def _text(element: ET.Element, tag: str) -> str | None:
    child = element.find(_BP + tag)
    return child.text.strip() if child is not None and child.text else None


def _protein_reference(
    protein: ET.Element, index: Mapping[str, ET.Element]
) -> dict[str, str | None]:
    """Accession, database and organism label of a ``Protein``'s reference."""
    reference = _one(protein, 'entityReference', index)
    xrefs = [_resolve(x, index) for x in reference.findall(_BP + 'xref')]
    ids = [x for x in xrefs if x.tag == _BP + 'UnificationXref' and _text(x, 'id')]
    if len(ids) > 1:
        raise ValueError(
            f'PhosphoSite BioPAX: {reference.get(_RDF + "ID")!r} has several accessions'
        )
    organism = reference.find(_BP + 'organism')
    return {
        'accession': _text(ids[0], 'id') if ids else None,
        'database': _text(ids[0], 'db') if ids else None,
        'organism': _text(_resolve(organism, index), 'displayName')
        if organism is not None
        else None,
    }


def _site(reaction: ET.Element, index: Mapping[str, ET.Element]) -> tuple[str, int]:
    """Residue and position of the single phosphorylated site of a reaction."""
    right = _one(reaction, 'right', index)
    features = [_resolve(f, index) for f in right.findall(_BP + 'feature')]
    if len(features) != 1:
        raise ValueError(
            f'PhosphoSite BioPAX: {reaction.get(_RDF + "ID")!r} has '
            f'{len(features)} features, expected 1'
        )
    feature = features[0]
    name = _FEATURE_REGEX.match(feature.get(_RDF + 'ID', ''))
    location = _one(feature, 'featureLocation', index)
    kind = _one(feature, 'modificationType', index)
    xref = _one(kind, 'xref', index)
    term = _PHOSPHO_RESIDUES.get(_text(xref, 'id') or '')
    if not name or not term:
        raise ValueError(
            f'PhosphoSite BioPAX: site not understood: {feature.get(_RDF + "ID")!r}'
        )
    label = (_text(kind, 'term') or '').lower()
    if (
        _text(location, 'positionStatus') != 'EQUAL'
        or _text(location, 'sequencePosition') != name[2]
        or term != (label, name[1])
    ):
        raise ValueError(
            f'PhosphoSite BioPAX: residue, position and modification term of '
            f'{feature.get(_RDF + "ID")!r} disagree'
        )
    return name[1], int(name[2])


def _evidence(
    reaction: ET.Element,
    index: Mapping[str, ET.Element],
    methods: Mapping[str, str],
) -> tuple[list[str], list[str], int]:
    """PubMed IDs, PSI-MI method terms and PubMed-less records of a reaction.

    Each list is distinct and in file order. A publication record without a
    PubMed ID (``recid_...``, a PhosphoSitePlus curation record) contributes
    only to the count of the third value.
    """
    pmids: list[str] = []
    terms: list[str] = []
    records: set[str] = set()
    for prop in reaction.findall(_BP + 'evidence'):
        evidence = _resolve(prop, index)
        for xref_prop in evidence.findall(_BP + 'xref'):
            xref = _resolve(xref_prop, index)
            if xref.tag != _BP + 'PublicationXref':
                raise ValueError(
                    f'PhosphoSite BioPAX: evidence {evidence.get(_RDF + "ID")!r} '
                    f'has a {xref.tag} reference'
                )
            database, number = _text(xref, 'db'), _text(xref, 'id')
            if database is None and number is None:
                records.add(xref.get(_RDF + 'ID', ''))
            elif database != 'pubmed':
                raise ValueError(
                    f'PhosphoSite BioPAX: publication {xref.get(_RDF + "ID")!r} '
                    f'is from database {database!r}'
                )
            elif number and is_pmid(number) and number not in pmids:
                pmids.append(number)
        for code_prop in evidence.findall(_BP + 'evidenceCode'):
            code = _resolve(code_prop, index)
            xref = _one(code, 'xref', index)
            term = _text(xref, 'id')
            if (
                _text(xref, 'db') != 'PSI-MI'
                or term not in methods
                or (_text(code, 'term') or '').lower() != methods[term].lower()
            ):
                raise ValueError(
                    f'PhosphoSite BioPAX: evidence method not understood: '
                    f'{term!r} {_text(code, "term")!r}'
                )
            if term not in terms:
                terms.append(term)
    return pmids, terms, len(records)


def iter_catalysis_evidence(
    opener, *, methods: Mapping[str, str], **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per ``Catalysis`` of the BioPAX file.

    ``methods`` maps the PSI-MI accessions the file may use to their labels as
    the file writes them (case aside); another method, a database other than
    PubMed for a publication, a site that is not a phosphorylated serine,
    threonine or tyrosine with matching name, position and term, and a
    reference to an object the file does not define are errors.

    Fields of a dictionary:

    - ``catalysis``: the ``rdf:ID`` of the catalysis (``Catalysis_43850_457``).
    - ``kinase``, ``substrate``: accession of the controller and of the
      modified protein as the file writes it, or None where its reference has
      none (two kinases of the file); ``kinase_database``, ``substrate_database``
      (``UniProtKB``, ``RefSeq``, ``GenBank``, ``Ensembl``) and
      ``kinase_organism``, ``substrate_organism`` (label of the file).
    - ``residue``, ``position``: the phosphorylated residue and its position.
    - ``publications``: the distinct PubMed IDs of the reaction's evidence.
    - ``methods``: the distinct PSI-MI accessions of the evidence codes.
    - ``records_without_pmid``: publication records of the evidence that have
      no PubMed ID, which are not carried.
    """
    root = read_xml_root(opener)
    if root is None:
        raise ValueError('PhosphoSite BioPAX file is empty')
    index: dict[str, ET.Element] = {}
    for element in root.iter():
        identifier = element.get(_RDF + 'ID')
        if identifier is not None:
            if identifier in index:
                raise ValueError(f'PhosphoSite BioPAX: {identifier!r} defined twice')
            index[identifier] = element
    reactions: dict[str, tuple[tuple[str, int], tuple]] = {}
    records = []
    for catalysis in root.findall(_BP + 'Catalysis'):
        controller = _protein_reference(_one(catalysis, 'controller', index), index)
        reaction = _one(catalysis, 'controlled', index)
        identifier = reaction.get(_RDF + 'ID')
        if reaction.tag != _BP + 'BiochemicalReaction':
            raise ValueError(f'PhosphoSite BioPAX: {identifier!r} is not a reaction')
        if identifier not in reactions:
            # The reaction, its site and its evidence are shared by the
            # catalyses of the site.
            substrate = _protein_reference(_one(reaction, 'right', index), index)
            reactions[identifier] = (
                _site(reaction, index),
                substrate,
                _evidence(reaction, index, methods),
            )
        (residue, position), substrate, (pmids, terms, others) = reactions[identifier]
        records.append(
            {
                'catalysis': catalysis.get(_RDF + 'ID'),
                'kinase': controller['accession'],
                'kinase_database': controller['database'],
                'kinase_organism': controller['organism'],
                'substrate': substrate['accession'],
                'substrate_database': substrate['database'],
                'substrate_organism': substrate['organism'],
                'residue': residue,
                'position': position,
                'publications': pmids,
                'methods': terms,
                'records_without_pmid': others,
            }
        )
    # The tree is the bulk of the memory; it is not needed to yield the records.
    del root, index, reactions
    if not records:
        raise ValueError('PhosphoSite BioPAX file has no catalysis')
    yield from records


def _site_key(
    kinase: Any, substrate: Any, residue: Any, position: Any
) -> tuple | None:
    """The identity of an observation; None if kinase or substrate has no accession."""
    kinase, substrate = str(kinase or '').strip(), str(substrate or '').strip()
    return (kinase, substrate, residue, position) if kinase and substrate else None


def join_evidence(
    rows: Iterable[dict[str, Any]],
    records: Iterable[dict[str, Any]],
    *,
    min_matched: float = MIN_MATCHED_FRACTION,
) -> Generator[dict[str, Any], None, None]:
    """Add the BioPAX publications and methods to the table rows.

    A table row and a catalysis record are the same observation when kinase
    accession, substrate accession, residue and position are equal, as the two
    files write them. Each row gets:

    - ``publications``, ``evidence_methods``: those of its catalysis record
      (empty lists where it has none);
    - ``evidence_record``: True where a catalysis record matched.

    Records that match no row are not emitted; rows that match no record are
    emitted without evidence. A record or row with no kinase accession or no
    substrate accession matches nothing. A key that two records share, or
    fewer matching records than ``min_matched`` of all (a sign that the files
    are of different releases), is an error.
    """
    by_key: dict[tuple, dict[str, Any]] = {}
    total = 0
    for record in records:
        key = _site_key(
            record['kinase'], record['substrate'], record['residue'], record['position']
        )
        total += 1
        if key is None:
            continue
        if key in by_key:
            raise ValueError(f'PhosphoSite BioPAX: two catalysis records for {key}')
        by_key[key] = record
    if not total:
        raise ValueError('PhosphoSite BioPAX file has no catalysis records')
    matched = set()
    for row in rows:
        key = _site_key(
            row['KIN_ACC_ID'], row['SUB_ACC_ID'], row['residue'], row['position']
        )
        record = by_key.get(key) if key else None
        row['evidence_record'] = record is not None
        row['publications'] = list(record['publications']) if record else []
        row['evidence_methods'] = list(record['methods']) if record else []
        if record:
            matched.add(key)
        yield row
    if len(matched) < min_matched * total:
        raise ValueError(
            f'PhosphoSite: only {len(matched)} of {total} catalysis records '
            'of the BioPAX file match a row of the table; the files may be of '
            'different releases'
        )


def iter_kinase_substrate_evidence(
    opener,
    *,
    organisms: Collection[str],
    methods: Mapping[str, str],
    evidence_download,
    min_matched: float = MIN_MATCHED_FRACTION,
    force_refresh: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """The rows of ``iter_kinase_substrate`` with the evidence of the BioPAX file.

    ``evidence_download`` is the ``Download`` of ``Kinase_substrates.owl.gz``,
    opened here with the same ``force_refresh``. See ``join_evidence``.
    """
    if not opener or not opener.result:
        return
    records = iter_catalysis_evidence(
        evidence_download.open(force_refresh=force_refresh), methods=methods
    )
    yield from join_evidence(
        iter_kinase_substrate(opener, organisms=organisms),
        records,
        min_matched=min_matched,
    )
