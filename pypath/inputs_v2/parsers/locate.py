"""Read the LOCATE subcellular localization database (release 6, 2008).

LOCATE is distributed as one XML file for each species (``LOCATE_human_v6_
20081121.xml``, ``LOCATE_mouse_v6_20081121.xml``) of 600 to 700 MB, with one
``LOCATE_protein`` element for each protein. This reader streams the file and
keeps only the parts it needs of each protein:

- ``protein``: ``organism`` (``Human`` or ``Mouse``), ``class`` (the membrane
  topology class of LOCATE: ``cytoplasmic``, ``mtmp``, ``secretome``,
  ``typeI``, ``typeII``) and a ``source`` with an accession; ``xrefs`` lists
  accessions of other databases. The ones read are the UniProt accessions and
  entry names, the Entrez Gene id, the Ensembl, RefSeq and GenBank protein
  accessions, the HGNC id (digits, 28,689 human proteins) and the MGI id
  (``MGI ID``, digits; the module writes ``MGI:<digits>``).
- ``literature/reference``: a publication (``source`` with the PubMed
  identifier as ``accn``; 0 where it has none), the ``organism`` it studied, a
  ``notes`` text and a ``locations/location`` whose ``goid`` attribute lists
  the GO identifiers of its ``tier1`` and ``tier2`` location names, separated
  by ``;``. This is the manually curated localization of the database.
- ``externalannot/reference``: a localization taken from another database: its
  ``locations/location`` (with ``goid``), the ``source`` that names the
  database (``source_name``) and the record of the protein in it (``accn``),
  and sometimes a GO ``evidence`` code. Its ``pmid`` is always 0.
- ``scl_prediction/source``: a localization predicted by a ``method`` (CELLO,
  MultiLoc, Proteome_Analyst, WoLFPSORT, pTarget), with the ``location``, its
  ``goid`` (several, separated by ``;``, where the method predicts a
  combination such as ``cytoplasm_nucleus``) and the method's ``evaluation``
  score.

``experimental_data`` (the images of the experiments of LOCATE) is empty in
the release. A localization without a GO identifier (11,478 predictions of
``No prediction`` and some of chloroplast and vacuole, and 2 literature
references that name only a location) has no row, and neither has a prediction
with the score 0, which is no prediction of the location.
"""

from __future__ import annotations

from collections.abc import Generator, Iterator
import re
from typing import Any
import xml.etree.ElementTree as ET

LITERATURE = 'literature'
EXTERNAL = 'external annotation'
PREDICTION = 'prediction'
PROTEIN_CLASS = 'protein class'

_GO = re.compile(r'GO:\d{7}')
# A UniProt accession (6 or 10 characters), and an entry name.
_ACCESSION = re.compile(
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)
_ENTRY = re.compile(r'^[A-Z0-9]+_[A-Z0-9]+$')
# The HGNC and MGI ids of the xrefs are written as digits.
_DIGITS = re.compile(r'^\d+$')

# The `source_name` of the accessions that identify a protein, with the key of
# the row's `ids` that holds them.
_IDS = {
    'UniProtKB-SwissProt': 'uniprot',
    'UniProt/SPTrEMBL': 'uniprot',
    'Entrez Gene': 'entrez',
    'Ensembl-Peptide Human': 'ensp',
    'Ensembl-Peptide Mouse': 'ensp',
    'RefSeq Protein': 'refseq_protein',
    'Entrez Protein': 'genbank_protein',
    # Identifiers of the gene of the protein record, like the Entrez Gene id.
    'HGNC': 'hgnc',
    'MGI ID': 'mgi',
}


def _text(element: ET.Element | None, path: str) -> str | None:
    value = element.findtext(path) if element is not None else None
    value = value.strip() if value else None
    return value or None


def _gos(value: str | None) -> list[str]:
    """The GO identifiers of a ``goid`` value, separated by ``;``; others raise."""
    gos = []
    for part in (value or '').split(';'):
        part = part.strip()
        if not part:
            continue
        if not _GO.fullmatch(part):
            raise ValueError(f'Not a GO identifier: {part!r}')
        if part not in gos:
            gos.append(part)
    return gos


def _identifiers(protein: ET.Element, uid: str) -> dict[str, Any]:
    """The accessions of a protein: ``uid`` (LOCATE's own) and the keys of ``_IDS``.

    UniProt accessions and entry names are told apart by their form: the
    ``UniProt/SPTrEMBL`` accession is an entry name for most proteins. Order
    is kept and repeats dropped.
    """
    ids: dict[str, list[str]] = {
        'uniprot': [], 'uniprot_entry': [], 'entrez': [], 'ensp': [],
        'refseq_protein': [], 'genbank_protein': [], 'hgnc': [], 'mgi': [],
    }

    def add(source_name: str | None, accession: str | None) -> None:
        key = _IDS.get(source_name or '')
        if key is None or not accession:
            return
        accession = accession.strip()
        if key == 'uniprot' and not _ACCESSION.match(accession):
            if not _ENTRY.match(accession):
                raise ValueError(f'Not a UniProt accession or entry name: {accession!r}')
            key = 'uniprot_entry'
        if key in ('hgnc', 'mgi'):
            if not _DIGITS.match(accession):
                raise ValueError(f'Not an {source_name} id: {accession!r}')
            # The MGI id is written with its prefix, as the Matrisome module writes it.
            if key == 'mgi':
                accession = f'MGI:{accession}'
        if accession not in ids[key]:
            ids[key].append(accession)

    # The accession the protein element carries comes first, then the xrefs;
    # the Swiss-Prot accessions are before the TrEMBL ones. The GenBank protein
    # accessions of the xrefs (up to nine for a protein) are not kept: only the
    # one of the protein element is.
    source = protein.find('protein/source')
    add(_text(source, 'source_name'), _text(source, 'accn'))
    xrefs = [
        (_text(x, 'source/source_name'), _text(x, 'source/accn'))
        for x in protein.findall('xrefs/xref')
        if _text(x, 'source/source_name') != 'Entrez Protein'
    ]
    for wanted in ('UniProtKB-SwissProt', 'UniProt/SPTrEMBL'):
        for name, accession in xrefs:
            if name == wanted:
                add(name, accession)
    for name, accession in xrefs:
        if name not in ('UniProtKB-SwissProt', 'UniProt/SPTrEMBL'):
            add(name, accession)
    return {'locate': uid, **ids}


def _protein_fields(element: ET.Element) -> dict[str, Any]:
    organism = _text(element, 'protein/organism')
    if organism is None:
        raise ValueError(f'Protein {element.get("uid")} has no organism')
    return {
        'ids': _identifiers(element, element.get('uid')),
        'organism': organism,
        'class': _text(element, 'protein/class'),
    }


def _statements(element: ET.Element) -> Iterator[dict[str, Any]]:
    """The localization rows of a protein: literature, external annotation, prediction."""
    protein = _protein_fields(element)
    for reference in element.findall('literature/reference'):
        location = reference.find('locations/location')
        pmid = _text(reference, 'source/accn')
        for go in _gos(location.get('goid') if location is not None else None):
            yield {
                **protein,
                'kind': LITERATURE,
                'go': go,
                'pmid': None if pmid in (None, '0') else pmid,
                'study_organism': _text(reference, 'organism'),
                'notes': _text(reference, 'notes'),
            }
    for reference in element.findall('externalannot/reference'):
        location = reference.find('locations/location')
        # The first `source` of a reference is empty; the second names the database.
        sources = [s for s in reference.findall('source') if _text(s, 'source_name')]
        if len(sources) != 1:
            raise ValueError(f'An annotation of {element.get("uid")} has {len(sources)} sources')
        for go in _gos(location.get('goid') if location is not None else None):
            yield {
                **protein,
                'kind': EXTERNAL,
                'go': go,
                'database': _text(sources[0], 'source_name'),
                'record': _text(sources[0], 'accn'),
                'evidence_code': _text(reference, 'evidence'),
            }
    for source in element.findall('scl_prediction/source'):
        score = float(_text(source, 'evaluation'))
        if score == 0.0:
            continue
        for go in _gos(_text(source, 'goid')):
            yield {
                **protein,
                'kind': PREDICTION,
                'go': go,
                'method': _text(source, 'method'),
                'location': _text(source, 'location'),
                'score': _text(source, 'evaluation'),
            }
    if protein['class']:
        yield {**protein, 'kind': PROTEIN_CLASS}


def _proteins(opener) -> Iterator[ET.Element]:
    """The ``LOCATE_protein`` elements of the file, one at a time."""
    result = getattr(opener, 'result', None) if opener else None
    if isinstance(result, dict):
        result = next(iter(result.values()), None)
    if result is None:
        raise ValueError('LOCATE: nothing to read')
    events = ET.iterparse(result, events=('start', 'end'))
    root = None
    count = 0
    for event, element in events:
        if event == 'start' and root is None:
            root = element
            if root.tag != 'LOCATE_doc':
                raise ValueError(f'Not a LOCATE file: root element {root.tag!r}')
        elif event == 'end' and element.tag == 'LOCATE_protein':
            count += 1
            yield element
            # Release what is read: the file is 600 MB.
            element.clear()
            root.clear()
    if not count:
        raise ValueError('LOCATE: the file has no protein')


def iter_localizations(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One row for each localization statement of each protein, and its class.

    Common fields: ``ids`` (a dict: ``locate``, ``uniprot``, ``uniprot_entry``,
    ``entrez``, ``ensp``, ``refseq_protein``, ``genbank_protein``, ``hgnc``,
    ``mgi`` (``MGI:<digits>``); lists, except
    ``locate``), ``organism`` and ``class`` (None where the protein has none),
    ``kind``. For ``literature``, ``external annotation`` and ``prediction``:
    ``go``; for literature also ``pmid`` (None for 0), ``study_organism`` and
    ``notes``; for external annotation ``database``, ``record`` and
    ``evidence_code``; for prediction ``method``, ``location`` (the method's
    label) and ``score`` (as written). For ``protein class``: only the common
    fields.
    """
    for element in _proteins(opener):
        yield from _statements(element)


def iter_protein_classes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The classes of LOCATE's proteins, each once: ``name``."""
    seen: set[str] = set()
    for element in _proteins(opener):
        name = _text(element, 'protein/class')
        if name and name not in seen:
            seen.add(name)
            yield {'name': name}
