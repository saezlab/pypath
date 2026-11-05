"""Raw parsers for BRENDA inputs_v2 datasets."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Generator
import re
from typing import Any, NamedTuple

from omnipath_resources.internals.ontology_schema import OntologyRelationship, OntologyTerm


ID = re.compile(r'^#([\d,]+)#')
REF_ID = re.compile(r'^\<(\d+)\>.*')
REF_PMID = re.compile(r'.*\{Pubmed:(\d+)\}')
EC_ID = re.compile(
    r'^(?P<ec>\d+\.\d+\.\d+\.(?:\d+|[nB]\d+))'
    r'(?:\s+\((?P<status>.*)\))?$'
)
EC_REFERENCE = re.compile(r'\bEC\s+(\d+\.\d+\.\d+\.(?:\d+|[nB]\d+))\b')
UNIPROT_IDENTIFIER = re.compile(
    r'\{([^{};]+);\s*source:\s*(?:UniProt|SwissProt|TrEMBL)\}',
    flags=re.IGNORECASE,
)
TRAILING_REFERENCE = re.compile(r'\s*<[\d,]+>\s*$')
TRAILING_REFERENCES = re.compile(r'\s*<([\d,\s]+)>\s*$')
REFERENCE_MARKER = re.compile(r'\s*<[\d,\s]+>')
NOTE_START = re.compile(r'\(\s*#[\d,]+#')
NOTE_SEPARATOR = re.compile(r';\s*(?=#[\d,]+#)')
ACCESSION_BLOCK = re.compile(r'\s*\{[^{}]*\}')


class Commentary(NamedTuple):
    """A ``#proteins# text <references>`` statement with its notes."""

    proteins: list[str]
    text: str
    notes: list[Commentary]
    references: list[str]


TOP_EC_CLASS_NAMES = {
    '1': 'Oxidoreductases',
    '2': 'Transferases',
    '3': 'Hydrolases',
    '4': 'Lyases',
    '5': 'Isomerases',
    '6': 'Ligases',
    '7': 'Translocases',
}


def iter_enzyme_ontology(opener, max_records: int | None = None, **_kwargs):
    records = []
    parents = set()
    for record in iter_brenda_records(opener):
        ec, _status = parse_ec_id(record.get('ID'))
        if not ec:
            continue
        records.append(record)
        parents.update(ec_parent_chain(ec))

    emitted = 0
    leaf_ids = {
        ec
        for record in records
        if (ec := parse_ec_id(record.get('ID'))[0])
    }
    for ec in sorted(parents - leaf_ids, key=ec_sort_key):
        yield term_to_row(synthetic_ec_term(ec))
        emitted += 1
        if max_records is not None and emitted >= max_records:
            return

    for record in records:
        term = ec_record_to_term(record)
        if term is None:
            continue
        yield term_to_row(term)
        emitted += 1
        if max_records is not None and emitted >= max_records:
            return


def iter_protein_enzyme_class_annotations(
    opener,
    max_records: int | None = None,
    **_kwargs,
):
    emitted = 0
    seen = set()
    for record in iter_brenda_records(opener, fields={'ID', 'PR'}):
        ec, status = parse_ec_id(record.get('ID'))
        if not ec or status:
            continue
        ec_term_id = ec_term_id_from_ec(ec)
        for protein_record in record.get('PR', []):
            for uniprot_id in extract_uniprot_ids(protein_record):
                key = (uniprot_id, ec_term_id)
                if key in seen:
                    continue
                seen.add(key)
                yield {
                    'UniProt': [uniprot_id],
                    'EC': ec,
                }
                emitted += 1
                if max_records is not None and emitted >= max_records:
                    return


def iter_molecular_forms(opener, max_records=None, **kwargs):
    """Keep each EN/PM observation with its explicit protein and citations.

    Official 2026.1 README: EN=engineering, PM=posttranslation modification;
    #protein# and <reference> delimiters are local to each EC record. Each
    protein keeps only the commentary and references the source assigns to it.
    Catalogue states are not broadcast to the separate kinetics datasets.
    """
    emitted = 0
    for record in iter_brenda_records(
        opener, fields={'ID', 'PR', 'RF', 'EN', 'PM'}
    ):
        ec, status = parse_ec_id(record.get('ID'))
        if not ec or status:
            continue
        proteins = {}
        for protein in record.get('PR', []):
            if parsed := parse_commentary(protein):
                organism = clean_text(ACCESSION_BLOCK.sub('', parsed.text))
                for pid in parsed.proteins:
                    proteins[pid] = (
                        protein,
                        organism or None,
                        [n.text for n in parsed.notes if pid in n.proteins],
                    )
        references = process_references(record)
        for section in ('EN', 'PM'):
            for observation in record.get(section, []):
                parsed = parse_commentary(observation)
                if not parsed:
                    continue
                noted = {ref for note in parsed.notes for ref in note.references}
                for pid in parsed.proteins:
                    own = [note for note in parsed.notes if pid in note.proteins]
                    own_notes = dedupe([note.text for note in own])
                    if not observation_texts(parsed.text, own_notes):
                        continue
                    own_refs = (
                        dedupe([ref for note in own for ref in note.references])
                        or [ref for ref in parsed.references if ref not in noted]
                        or parsed.references
                    )
                    protein, organism, protein_notes = proteins.get(
                        pid, ('', None, [])
                    )
                    accessions = extract_uniprot_ids(protein)
                    yield {
                        'EC': ec,
                        'protein_record_id': pid,
                        'UniProt': accessions[0]
                        if len(accessions) == 1
                        else None,
                        'organism': organism,
                        'protein_notes': protein_notes,
                        'source_accessions': accessions,
                        'section': section,
                        'descriptor': parsed.text,
                        'notes': own_notes,
                        'Refs': [
                            references[ref]
                            for ref in own_refs
                            if references.get(ref)
                        ],
                    }
                    emitted += 1
                    if max_records is not None and emitted >= max_records:
                        return


def parse_commentary(value: str) -> Commentary | None:
    """Split BRENDA ``#1,2# text (#1# note <3>; #2# note <4>) <3,4>`` lines.

    Markers are removed from all texts; protein and reference numbers are
    local to the EC record.
    """
    match = ID.match(value)
    if not match:
        return None
    pids = split_ref_ids(match[1])
    body = value[match.end():].strip()
    refs = []
    if trailing := TRAILING_REFERENCES.search(body):
        refs = split_ref_ids(trailing[1])
        body = body[: trailing.start()]
    notes = []
    if body.endswith(')') and (start := NOTE_START.search(body)):
        commentary = body[start.start() + 1 : -1]
        # Some source lines close the commentary twice.
        while commentary.endswith(')') and commentary.count(')') > commentary.count('('):
            commentary = commentary[:-1]
        for entry in NOTE_SEPARATOR.split(commentary):
            note = parse_commentary(entry.strip()) or Commentary(pids, entry, [], [])
            if text := clean_text(REFERENCE_MARKER.sub('', note.text)):
                notes.append(note._replace(text=text))
        body = body[: start.start()]
    return Commentary(pids, clean_text(REFERENCE_MARKER.sub('', body)), notes, refs)


def observation_texts(descriptor: str, notes: list[str]) -> list[str]:
    """One readable statement per note; ``more`` is BRENDA's free-text placeholder."""
    if descriptor.lower() in {'more', 'additional information'}:
        return notes
    return [f'{descriptor}: {note}' for note in notes] or [descriptor]


def term_record_to_term(row: dict[str, Any]) -> OntologyTerm | None:
    if not row.get('id'):
        return None

    replaced_by = row.get('replaced_by') or []
    relationships = [
        OntologyRelationship(type='replaced_by', target=target)
        for target in replaced_by
    ]

    return OntologyTerm(
        id=row['id'],
        name=row['name'],
        definition=row.get('definition') or None,
        synonyms=row.get('synonyms') or None,
        comments=row.get('comments') or None,
        is_a=row.get('is_a') or None,
        relationships=relationships or None,
        is_obsolete=str(row.get('is_obsolete', '')).lower() == 'true',
    )


def ec_term_id(row: dict[str, Any]) -> str:
    return ec_term_id_from_ec(row['EC'])


def process_references(record):
    refs = {}
    regexes = [
        REF_ID,
        REF_PMID,
    ]

    for entry in record.get('RF', []):
        res = [
            x.group(1) if (x := regex.match(entry)) else ''
            for regex in regexes
        ]
        refs[res[0]] = res[1]

    return refs


def term_to_row(term: OntologyTerm) -> dict[str, Any]:
    relationships = term.relationships or []
    return {
        'id': term.id,
        'name': term.name,
        'definition': term.definition or '',
        'synonyms': term.synonyms or [],
        'comments': term.comments or [],
        'is_a': term.is_a or [],
        'replaced_by': [
            relationship.target
            for relationship in relationships
            if relationship.type == 'replaced_by'
        ],
        'is_obsolete': 'true' if term.is_obsolete else 'false',
    }


def split_ref_ids(value: str) -> list[str]:
    return [item.strip() for item in value.split(',') if item.strip()]


def read_brenda_text(opener) -> str:
    if not opener or not getattr(opener, 'result', None):
        return ''

    keys = sorted(opener.result)
    handle = opener.result[keys[0]]
    if hasattr(handle, 'seek'):
        handle.seek(0)
    content = handle.read()
    return content.decode('utf-8', 'ignore') if isinstance(content, bytes) else str(content)


def iter_brenda_records(
    opener,
    *,
    fields: set[str] | None = None,
) -> Generator[dict[str, Any], None, None]:
    if not opener or not getattr(opener, 'result', None):
        return
    handle = opener.result[sorted(opener.result)[0]]
    if hasattr(handle, 'seek'):
        handle.seek(0)
    record = defaultdict(list)
    previous_key = None
    for raw_line in handle:
        line = raw_line.decode('utf-8') if isinstance(raw_line, bytes) else raw_line
        line = line.rstrip('\r\n')
        if line.startswith('///'):
            if record.get('ID'):
                yield dict(record)
            record = defaultdict(list)
            previous_key = None
        elif line[:1].isspace() and previous_key:
            if previous_key == 'ID':
                record['ID'] += ' ' + line.strip()
            else:
                record[previous_key][-1] += ' ' + line.strip()
        elif '\t' in line:
            key, value = line.split('\t', 1)
            previous_key = key if fields is None or key in fields else None
            if previous_key:
                if key == 'ID':
                    record[key] = value.strip()
                else:
                    record[key].append(value.strip())
    if record.get('ID'):
        yield dict(record)


def parse_ec_id(value: str | None) -> tuple[str | None, str | None]:
    if not value:
        return None, None
    match = EC_ID.match(value.strip().rstrip(','))
    if not match:
        return None, None
    return match.group('ec'), match.group('status')


def ec_term_id_from_ec(ec: str) -> str:
    return f'EC:{ec}'


def ec_parent(ec: str) -> str | None:
    parts = ec.split('.')
    if len(parts) <= 1:
        return None
    return '.'.join(parts[:-1])


def ec_parent_chain(ec: str) -> list[str]:
    parts = ec.split('.')
    return ['.'.join(parts[:i]) for i in range(1, len(parts))]


def ec_sort_key(ec: str) -> list[tuple[int, int | str]]:
    key = []
    for part in ec.split('.'):
        key.append((0, int(part)) if part.isdigit() else (1, part))
    return key


def first_nonempty(values: list[str] | None) -> str | None:
    for value in values or []:
        value = value.strip()
        if value:
            return value
    return None


def clean_text(value: str) -> str:
    return re.sub(r'\s+', ' ', value).strip()


def clean_reaction_definition(value: str) -> str:
    value = value.split(' (#', maxsplit=1)[0]
    value = TRAILING_REFERENCE.sub('', value)
    return clean_text(value)


def clean_synonym(value: str) -> str | None:
    text = TRAILING_REFERENCE.sub('', value).strip()
    if not text or text.startswith('#'):
        return None
    return clean_text(text)


def dedupe(values: list[str]) -> list[str]:
    seen = set()
    result = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            result.append(value)
    return result


def synthetic_ec_term(ec: str) -> OntologyTerm:
    parent = ec_parent(ec)
    return OntologyTerm(
        id=ec_term_id_from_ec(ec),
        name=TOP_EC_CLASS_NAMES.get(ec, f'EC {ec}'),
        is_a=[ec_term_id_from_ec(parent)] if parent else None,
        is_obsolete=False,
    )


def ec_record_to_term(record: dict[str, Any]) -> OntologyTerm | None:
    ec, status = parse_ec_id(record.get('ID'))
    if not ec:
        return None

    name = first_nonempty(record.get('RN')) or f'EC {ec}'
    systematic_name = first_nonempty(record.get('SN'))
    synonyms = []
    if systematic_name and systematic_name != name:
        synonyms.append(systematic_name)
    synonyms.extend(
        synonym
        for value in record.get('SY', [])
        if (synonym := clean_synonym(value)) and synonym != name
    )

    reactions = [
        reaction
        for value in record.get('RE', [])
        if (reaction := clean_reaction_definition(value))
    ]
    parent = ec_parent(ec)
    relationships = []
    if status and status.lower().startswith('transferred'):
        relationships.extend(
            OntologyRelationship(
                type='replaced_by',
                target=ec_term_id_from_ec(target),
            )
            for target in EC_REFERENCE.findall(status)
            if target != ec
        )

    return OntologyTerm(
        id=ec_term_id_from_ec(ec),
        name=name,
        definition='; '.join(dedupe(reactions)) or None,
        synonyms=dedupe(synonyms) or None,
        comments=[status] if status else None,
        is_a=[ec_term_id_from_ec(parent)] if parent else None,
        relationships=relationships or None,
        is_obsolete=bool(status),
    )


def extract_uniprot_ids(value: str) -> list[str]:
    identifiers = []
    for identifier in UNIPROT_IDENTIFIER.findall(value):
        identifiers.extend(
            item.strip()
            for item in identifier.split(' AND ')
            if item.strip()
        )
    return sorted(set(identifiers))


# ── regulation ───────────────────────────────────────────────────────────────

# The roles of a compound in BRENDA's AC, IN and CF fields.
REGULATION_ROLES = {'AC': 'activator', 'IN': 'inhibitor', 'CF': 'cofactor'}
# A Ki line: ``#1,2# 0.0051 {pyrazole} (notes) <refs>``; -999 means "no value".
KI_VALUE = re.compile(r'^#([\d,]+)#\s+(\S+)\s+\{(.+?)\}')
_SPACE = re.compile(r'\s+')


def _tsv_value(value: str) -> str:
    """A SPARQL TSV cell: a quoted literal is unquoted, an IRI or number stays."""
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        value = value[1:-1].replace('\\"', '"').replace('\\t', '\t').replace('\\\\', '\\')
    return value


def _iter_tsv(opener) -> Generator[list[str], None, None]:
    handle = opener.result if not isinstance(opener.result, dict) else opener.result[sorted(opener.result)[0]]
    lines = iter(handle)
    next(lines, None)  # the header names the query variables
    for line in lines:
        line = line.decode('utf-8') if isinstance(line, bytes) else line
        if line.strip():
            yield [_tsv_value(cell) for cell in line.rstrip('\r\n').split('\t')]


class BrendaCompounds:
    """BRENDA's own ligand table: names and synonyms to ligand IDs and structures.

    A BRENDA name finds its ligand when it names exactly one: a label first,
    then any name, then the same ignoring whitespace (the flat file wraps long
    names and gains spaces where it does). A name of several ligands finds none.
    """

    def __init__(self, names_opener, structures_opener) -> None:
        labels, names = defaultdict(set), defaultdict(set)
        squashed_labels, squashed = defaultdict(set), defaultdict(set)
        for ligand, kind, name in _iter_tsv(names_opener):
            key = _SPACE.sub('', name)
            names[name].add(ligand)
            squashed[key].add(ligand)
            if kind == 'label':
                labels[name].add(ligand)
                squashed_labels[key].add(ligand)
        self._lookups = (labels, names, squashed_labels, squashed)
        self.structures = {
            ligand: (inchikey.removeprefix('InChIKey='), inchi)
            for ligand, inchikey, inchi in _iter_tsv(structures_opener)
        }

    def ligand(self, name: str) -> str | None:
        for lookup, key in zip(self._lookups, (name, name, _SPACE.sub('', name), _SPACE.sub('', name))):
            found = lookup.get(key, ())
            if len(found) == 1:
                return next(iter(found))
            if found:
                return None
        return None


def _ki_values(record: dict) -> dict[tuple[str, str], list[str]]:
    """Ki values by protein and compound name; BRENDA's -999 is "no value"."""
    values = defaultdict(list)
    for line in record.get('KI', []):
        if not (match := KI_VALUE.match(line)) or match[2].startswith('-999'):
            continue
        compound = clean_text(match[3])
        for pid in split_ref_ids(match[1]):
            values[pid, compound].append(match[2])
    return values


def iter_regulation(
    opener,
    *,
    names_download,
    structures_download,
    max_records: int | None = None,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """One record per protein and compound of an AC, IN or CF statement.

    Each protein keeps only the notes and references the source assigns to it,
    as the molecular forms do. The compound is BRENDA's own ligand where its
    name finds one (``BrendaCompounds``); Ki values of the same protein and
    compound go with it.
    """
    if not opener or not getattr(opener, 'result', None):
        return
    force_refresh = kwargs.get('force_refresh', False)
    compounds = BrendaCompounds(
        names_download.open(force_refresh=force_refresh),
        structures_download.open(force_refresh=force_refresh),
    )
    emitted = 0
    for record in iter_brenda_records(
        opener, fields={'ID', 'PR', 'RF', 'AC', 'IN', 'CF', 'KI'}
    ):
        ec, status = parse_ec_id(record.get('ID'))
        if not ec or status:
            continue
        proteins = {}
        for protein in record.get('PR', []):
            if parsed := parse_commentary(protein):
                organism = clean_text(ACCESSION_BLOCK.sub('', parsed.text))
                accessions = extract_uniprot_ids(protein)
                for pid in parsed.proteins:
                    proteins[pid] = (organism or None, accessions)
        references = process_references(record)
        ki = _ki_values(record)
        for section, role in REGULATION_ROLES.items():
            for statement in record.get(section, []):
                parsed = parse_commentary(statement)
                if not parsed or parsed.text.lower() in {'more', 'additional information'}:
                    continue
                ligand = compounds.ligand(parsed.text)
                inchikey, inchi = compounds.structures.get(ligand, (None, None))
                noted = {ref for note in parsed.notes for ref in note.references}
                for pid in parsed.proteins:
                    own = [note for note in parsed.notes if pid in note.proteins]
                    own_refs = (
                        dedupe([ref for note in own for ref in note.references])
                        or [ref for ref in parsed.references if ref not in noted]
                        or parsed.references
                    )
                    organism, accessions = proteins.get(pid, (None, []))
                    yield {
                        'EC': ec,
                        'protein_record_id': pid,
                        'UniProt': accessions[0] if len(accessions) == 1 else None,
                        'organism': organism,
                        'role': role,
                        'compound_name': parsed.text,
                        'brenda_ligand': ligand,
                        'inchikey': inchikey,
                        'inchi': inchi,
                        'notes': dedupe([note.text for note in own]),
                        'ki': ki.get((pid, parsed.text), []) if role == 'inhibitor' else [],
                        'Refs': [references[ref] for ref in own_refs if references.get(ref)],
                    }
                    emitted += 1
                    if max_records is not None and emitted >= max_records:
                        return
