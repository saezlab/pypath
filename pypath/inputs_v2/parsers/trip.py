"""
Parser of the TRIP Database export and of the TRIP protein list.

The export is one CSV file with a five-line header and 32 columns, read with
the shared CSV reader (:func:`iter_csv`) under fixed column names. The site's
protein list is an HTML page with two tables of the same layout (TRP channels
and interacting proteins), read with the shared table reader
(:func:`iter_html_tables`); it gives the UniProt accession of each protein name
in human, mouse and rat. :func:`iter_trip_observations` joins the two: the
accession of a participant is looked up by its name and its species.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv, iter_html_tables, read_text

# The 32 columns of the export; its header spans five lines.
_COLUMNS = (
    'channel',
    'interactor',
    'method',
    # Screening
    'screening_channel_species',
    'screening_channel_region',
    'screening_interactor_species',
    'screening_interactor_tissue',
    'screening_sample_type',
    'screening_non_experimental',
    # In vitro validation
    'in_vitro_channel_species',
    'in_vitro_channel_region',
    'in_vitro_channel_expression',
    'in_vitro_interactor_species',
    'in_vitro_interactor_region',
    'in_vitro_interactor_expression',
    # In vivo validation
    'in_vivo_endogenous_cell',
    'in_vivo_overexpressed_cell',
    'in_vivo_channel_species',
    'in_vivo_channel_region',
    'in_vivo_interactor_species',
    'in_vivo_interactor_region',
    # Characterization
    'binding_channel_species',
    'binding_channel_region',
    'binding_interactor_species',
    'binding_interactor_region',
    'stoichiometry',
    'affinity',
    # Functional consequence
    'ptm',
    'trafficking',
    'activity',
    'reference',
    'overlap',
)
_HEADER_LINES = 5
# Cells that identify the layout, by header line and column.
_HEADER_ANCHORS = {
    0: {
        'screening_channel_species': 'Screening',
        'in_vitro_channel_species': 'Validation',
        'binding_channel_species': 'Characterization',
        'ptm': 'Functional consequence',
        'reference': 'Reference',
    },
    1: {
        'screening_non_experimental': 'Non-experimental screening',
        'in_vivo_endogenous_cell': 'In vivo validation',
        'stoichiometry': 'Stoichiometry',
        'affinity': 'Affinity (Kd)',
        'ptm': 'Post-translational modification',
        'trafficking': 'Subcellular trafficking',
        'activity': 'Activity',
    },
    4: {
        'channel': 'TRP channels',
        'interactor': 'Interactors',
        'method': 'Methods',
        'overlap': 'Overlap',
    },
}
_CATEGORY_COLUMNS = {
    'screening': _COLUMNS[3:9],
    'in_vitro_validation': _COLUMNS[9:15],
    'in_vivo_validation': _COLUMNS[15:21],
    'characterization': _COLUMNS[21:27],
    'functional_consequence': _COLUMNS[27:30],
}
_SPECIES_COLUMNS = {
    'channel': (
        'screening_channel_species',
        'in_vitro_channel_species',
        'in_vivo_channel_species',
        'binding_channel_species',
    ),
    'interactor': (
        'screening_interactor_species',
        'in_vitro_interactor_species',
        'in_vivo_interactor_species',
        'binding_interactor_species',
    ),
}

# Species labels of the species columns -> NCBI Taxonomy ID, which the protein
# list is keyed by. "Xenopus" names a genus, not a species.
SPECIES_TAXA = {
    'Human': '9606',
    'Mouse': '10090',
    'Rat': '10116',
    'Rabbit': '9986',
    'Bovine': '9913',
    'Chicken': '9031',
    'Pig': '9823',
    'Guinea pig': '10141',
    'Xenopus': None,
}
# Values of the species columns that name no species.
_NO_SPECIES = frozenset(
    {
        'Not specified',
        'Not specifed',
        'Not determined',
        'Not used',
        'Not used as a bait',
        'Endogenous protein',
    }
)
# Columns of the protein list after the name, by taxon.
_PROTEIN_LIST_TAXA = ('9606', '10090', '10116')
_UNIPROT_RE = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)
# "Phosphorylation (Tyr-273)": a modification and the channel residue.
_PTM_RE = re.compile(r'([A-Za-z-]+)(?: \(((?:[A-Za-z]{3})-(\d+))\))?')


def repair_greek(text: str) -> str:
    """Restore the Greek letters of TRIP protein names.

    The export and the later site pages write each Greek letter 0x210 code
    points too low ("CaMKIIƢ" for "CaMKIIβ"); older captures of the
    same pages have the Greek letters.
    """
    return ''.join(
        chr(ord(char) + 0x210) if 0x181 <= ord(char) <= 0x1B9 else char
        for char in text
    )


def protein_accessions(html: str) -> dict[str, dict[str, str]]:
    """UniProt accessions of the TRIP protein list, by name and taxon.

    The page has two tables of the same layout (TRP channels and interacting
    proteins): a name and its UniprotKB ID in human, mouse and rat.
    """
    accessions: dict[str, dict[str, str]] = {}
    for table in iter_html_tables(html):
        header = [cell['text'] for row in table['rows'] for cell in row if cell['header']]
        if header != ['Name', 'UniprotKB ID', 'Human', 'Mouse', 'Rat']:
            continue
        for row in table['rows']:
            if len(row) != 4 or any(cell['header'] for cell in row):
                continue
            name = repair_greek(row[0]['text'])
            by_taxon = {
                taxon: cell['text']
                for taxon, cell in zip(_PROTEIN_LIST_TAXA, row[1:])
                if _UNIPROT_RE.fullmatch(cell['text'])
            }
            if name in accessions and accessions[name] != by_taxon:
                raise ValueError(f'TRIP protein list names {name!r} twice')
            accessions[name] = by_taxon
    if not accessions:
        raise ValueError('TRIP: the protein list has no UniprotKB ID table')
    return accessions


def _species(row: dict[str, str], role: str) -> str | None:
    """The species of a participant, if the line names exactly one."""
    labels = set()
    for column in _SPECIES_COLUMNS[role]:
        label = row[column]
        if not label or label in _NO_SPECIES:
            continue
        if label not in SPECIES_TAXA:
            raise ValueError(f'Unknown TRIP species: {label!r} in {column}')
        labels.add(label)
    return labels.pop() if len(labels) == 1 else None


def _mirror(row: dict[str, str]) -> tuple[str, ...]:
    """The line as the other TRP channel's table lists it."""
    swapped = dict(row)
    for column in _COLUMNS:
        other = column.replace('channel', 'interactor')
        if other != column and other in row:
            swapped[column], swapped[other] = row[other], row[column]
    return tuple(swapped[column] for column in _COLUMNS)


def iter_trip_table(opener, **kwargs: Any) -> Generator[dict[str, str], None, None]:
    """Yield the data lines of the export, values stripped, keyed by column."""
    # The header lines are read as rows and checked below; a line with
    # another number of fields is an error.
    for number, row in enumerate(
        iter_csv(opener, fieldnames=_COLUMNS, strict=True, **kwargs), 1
    ):
        row = {column: value.strip() for column, value in row.items()}
        if number <= _HEADER_LINES:
            for column, expected in _HEADER_ANCHORS.get(number - 1, {}).items():
                if row[column] != expected:
                    raise ValueError(
                        f'TRIP: header line {number} has {row[column]!r} '
                        f'where {expected!r} was expected'
                    )
            continue
        row['line'] = str(number)
        yield row


def iter_trip_observations(
    opener,
    force_refresh: bool = False,
    proteins_download=None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per observation of the TRIP export.

    ``proteins_download`` is the download of the protein list page.

    Output fields, besides the 32 source columns (stripped, otherwise as
    written) and ``line``, the line number in the file:

    - ``channel_name``, ``interactor_name``: the names with their Greek
      letters restored (see :func:`repair_greek`).
    - ``categories``: the evidence categories whose columns are filled.
    - ``channel_species``, ``interactor_species``: the species label of the
      participant if the line names exactly one, else None.
    - ``channel_taxon``, ``interactor_taxon``: its NCBI Taxonomy ID.
    - ``channel_uniprot``, ``interactor_uniprot``: the UniProt accession the
      TRIP protein list gives for this name in this species (human, mouse or
      rat), else None. No accession is chosen when the species is not named.
    - ``ptm_type``, ``ptm_residue``: the two parts of the ``ptm`` column
      (``Phosphorylation`` and ``Tyr-273``); ``ptm_amino_acid`` and
      ``ptm_position``: the residue split into the three-letter name of the
      amino acid (``Tyr``) and the position as an integer (``273``). All are
      None where the cell has no residue.
    - ``mirrored``: True for the second listing of an observation between two
      TRP channels. TRIP lists these in the table of each channel (column
      ``overlap``); the second of two lines that differ only by the swap of
      channel and interactor is marked. Functional consequences are directed
      and never marked.
    """
    if proteins_download is None:
        raise ValueError('TRIP: the download of the protein list is needed')
    accessions = protein_accessions(
        read_text(proteins_download.open(force_refresh=force_refresh))
    )
    yield from _observations(iter_trip_table(opener), accessions)


def _observations(rows, accessions: dict[str, dict[str, str]]):
    awaited: Counter = Counter()
    for row in rows:
        record: dict[str, Any] = dict(row)
        record['categories'] = [
            category
            for category, columns in _CATEGORY_COLUMNS.items()
            if any(row[column] for column in columns)
        ]
        for role in ('channel', 'interactor'):
            name = repair_greek(row[role])
            species = _species(row, role)
            taxon = SPECIES_TAXA[species] if species else None
            record[f'{role}_name'] = name
            record[f'{role}_species'] = species
            record[f'{role}_taxon'] = taxon
            record[f'{role}_uniprot'] = accessions.get(name, {}).get(taxon)
        match = _PTM_RE.fullmatch(row['ptm']) if row['ptm'] else None
        record['ptm_type'] = match[1] if match else None
        record['ptm_residue'] = match[2] if match else None
        record['ptm_amino_acid'] = match[2][:3] if match and match[2] else None
        record['ptm_position'] = (
            int(match[3]) if match and match[3] else None
        )
        record['mirrored'] = False
        if row['overlap'] and 'functional_consequence' not in record['categories']:
            key = tuple(row[column] for column in _COLUMNS)
            if awaited[key]:
                awaited[key] -= 1
                record['mirrored'] = True
            else:
                awaited[_mirror(row)] += 1
        yield record
