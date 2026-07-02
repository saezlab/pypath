"""Raw parsers for BRENDA inputs_v2 datasets."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Generator
import re
from typing import Any

from pypath.internals.ontology_schema import OntologyRelationship, OntologyTerm

RECORDS_ENABLED = {
    'ID',
    'PR',
    'AC',
    'IN',
    'CF',
    'KI',
    'KM',
    'RF',
}

ID = re.compile(r'^#([\d,]+)#')
PR_ORGANISM_NAME = re.compile(r"#\d+# '?([-\w\s\.\[\]]+[^\s#\{<\('])")
PR_IDENTIFIER = re.compile(r'.*\{(.*?)\;.*')
REFERENCE = re.compile(r'.*\<([\d,\s]+)\>$')
ROLE_COMPOUND = re.compile(r'^#[\d,]+# ([^#]+) [\(\<].*')
K_COMPOUND = re.compile(r'^#[\d,]+# ([-\d.e]+ \{.*\}).*')
DETAILS = re.compile(r'.*\((#.*\>)\).*')
REF_ID = re.compile(r'^\<(\d+)\>.*')
REF_PMID = re.compile(r'.*\{Pubmed:(\d+)\}')
EC_ID = re.compile(
    r'^(?P<ec>\d+\.\d+\.\d+\.(?:\d+|[nB]\d+))'
    r'(?:\s+\((?P<status>.*)\))?$'
)
EC_REFERENCE = re.compile(r'\bEC\s+(\d+\.\d+\.\d+\.(?:\d+|[nB]\d+))\b')
UNIPROT_IDENTIFIER = re.compile(
    r'\{([^{};]+);\s*source:\s*(?:UniProt|SwissProt)\}',
    flags=re.IGNORECASE,
)
TRAILING_REFERENCE = re.compile(r'\s*<[\d,]+>\s*$')

TOP_EC_CLASS_NAMES = {
    '1': 'Oxidoreductases',
    '2': 'Transferases',
    '3': 'Hydrolases',
    '4': 'Lyases',
    '5': 'Isomerases',
    '6': 'Ligases',
    '7': 'Translocases',
}

ROLES_MAPPER = {
    'AC': ('Activator', ROLE_COMPOUND),
    'IN': ('Inhibitor', ROLE_COMPOUND),
    'CF': ('Cofactor', ROLE_COMPOUND),
    'KI': ('Ki', K_COMPOUND),
    'KM': ('Km', K_COMPOUND),
}


def parser(opener, **_kwargs):
    for record in iter_brenda_records(opener, fields=RECORDS_ENABLED):
        ec, status = parse_ec_id(record.get('ID'))
        if not ec or status:
            continue

        parsed_record = defaultdict(list, record)
        parsed_record['ID'] = ec
        for row in process_record(parsed_record).values():
            yield normalize_processed_record(row)


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


def process_record(record):
    eid = record['ID']
    ref_dict = process_references(record)
    proc = defaultdict(lambda: defaultdict(set))

    for pr in record.get('PR', []):
        id_match = ID.match(pr)
        organism_match = PR_ORGANISM_NAME.match(pr)
        reference_match = REFERENCE.match(pr)
        if not id_match or not organism_match or not reference_match:
            continue

        pid = id_match.group(1)
        org = organism_match.group(1)
        upids = g.group(1) if (g := PR_IDENTIFIER.match(pr)) else ''
        upids = [i for i in upids.split(' AND ') if i]
        refs = split_ref_ids(reference_match.group(1))

        proc[pid]['EC'] = eid
        proc[pid]['UniProt'].update(upids)
        proc[pid]['#'] = pid
        proc[pid]['Organism'].add(org)
        proc[pid]['Refs'].update(refs)

    for key in ['AC', 'IN', 'CF', 'KI', 'KM']:
        proc = process_record_roles(proc, record, key)

    for pid in proc.keys():
        proc[pid]['Refs'] = [
            ref
            for i in proc[pid]['Refs']
            if (ref := ref_dict[i])
        ]

    return proc


def process_record_roles(proc, record, key):
    if key not in ROLES_MAPPER:
        return {'ERROR': '`key` not found in `ROLES_MAPPER` for processing'}

    new_key, comp_regex = ROLES_MAPPER[key]

    for role_record in record.get(key, []):
        compound_match = comp_regex.match(role_record)
        id_match = ID.match(role_record)
        reference_match = REFERENCE.match(role_record)
        if not compound_match or not id_match or not reference_match:
            continue

        compound = compound_match.group(1)
        aux = x.group(1) if (x := DETAILS.match(role_record)) else ''
        initial_pids = set(id_match.group(1).split(','))
        initial_refs = set(split_ref_ids(reference_match.group(1)))

        aux = re.sub(r'; #', r'||#', aux)

        for entry in aux.split('||'):
            if not entry:
                continue

            entry_id_match = ID.match(entry)
            entry_reference_match = REFERENCE.match(entry)
            if not entry_id_match or not entry_reference_match:
                continue

            pids = entry_id_match.group(1).split(',')
            refs = split_ref_ids(entry_reference_match.group(1))

            for pid in pids:
                proc[pid][new_key].add(compound)
                proc[pid]['Refs'].update(refs)
                initial_pids -= {pid}

            initial_refs -= set(refs)

        for pid in initial_pids:
            proc[pid][new_key].add(compound)
            proc[pid]['Refs'].update(initial_refs)

    return proc


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


def normalize_processed_record(row: dict[str, Any]) -> dict[str, Any]:
    return {
        key: sorted(value) if isinstance(value, set) else value
        for key, value in row.items()
    }


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
    text = read_brenda_text(opener)
    for entry in text.replace('\r\n', '\n').replace('\n\t', ' ').split('///'):
        record: dict[str, Any] = defaultdict(list)
        for line in entry.split('\n'):
            if not line or '\t' not in line:
                continue
            key, value = line.split('\t', maxsplit=1)
            if fields is not None and key not in fields:
                continue
            value = value.strip()
            if key == 'ID':
                record[key] = value
            else:
                record[key].append(value)

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
