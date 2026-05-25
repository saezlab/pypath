"""KEGG REST bulk-download parser for metabolic reactions, enzymes and compounds."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from collections.abc import Generator
from typing import Any


_DATA_CACHE: dict[str, dict[str, list[dict[str, Any]]]] = {}

_ORGANISM_TAXA = {
    'hsa': '9606',
    'mmu': '10090',
    'rno': '10116',
    'dre': '7955',
    'cel': '6239',
    'dme': '7227',
    'sce': '4932',
}

_ARROW_RE = re.compile(r'<=>|=>|<=')
_STOICH_RE = re.compile(r'^(?:>\s*)?(?:\d+|[a-zA-Z]|\([^)]+\))\s+')
_EQUATION_PART_RE = re.compile(
    r'^(?:(?P<stoich>\d+|[A-Za-z]|\([^)]+\))\s+)?(?P<compound>[CG]\d{5})$'
)


def _raw(
    opener,
    data_type: str,
    organism: str = 'mmu',
    force_refresh: bool = False,
    max_records: int | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield normalised KEGG records for the requested dataset.

    Parameters
    ----------
    opener:
        A pypath curl opener whose ``result`` attribute is a dict of
        file-like objects keyed by filename, matching the multi-file
        download pattern used by dlmachine.
    data_type:
        One of ``'reactions'``.  Additional types (e.g. ``'compounds'``)
        may be added here as the module evolves.
    organism:
        KEGG organism code, e.g. ``'mmu'`` (mouse) or ``'hsa'`` (human).
    force_refresh:
        Bypass the in-process cache and re-parse from the opener.
    max_records:
        If set, yield at most this many records (useful during development).
    """
    records = _load_records(opener, organism=organism, force_refresh=force_refresh)
    rows = records.get(data_type, [])

    if max_records is not None:
        rows = rows[:max_records]

    yield from rows


def _load_records(
    opener,
    organism: str = 'mmu',
    force_refresh: bool = False,
) -> dict[str, list[dict[str, Any]]]:

    cache_key = organism
    if _DATA_CACHE.get(cache_key) and not force_refresh:
        return _DATA_CACHE[cache_key]

    handles = _get_handles(opener)
    records = _build_records(handles, organism=organism)

    _DATA_CACHE[cache_key] = records
    return records


def _get_handles(opener) -> dict[str, Any]:
    """Extract named file handles from the opener's result."""
    if not opener or not opener.result:
        return {}
    if isinstance(opener.result, dict):
        return opener.result
    # Single-file opener — caller passed the wrong opener type; return empty
    return {}


def _iter_tsv(handle) -> Generator[tuple[str, str], None, None]:
    """Yield (col1, col2) pairs from a two-column KEGG TSV response."""
    if not handle:
        return
    if hasattr(handle, 'seek'):
        handle.seek(0)
    for raw_line in handle:
        if isinstance(raw_line, bytes):
            raw_line = raw_line.decode('utf-8', 'ignore')
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split('\t', 1)
        if len(parts) == 2:
            yield parts[0].strip(), parts[1].strip()


def _strip_prefix(value: str, prefix: str) -> str:
    """Remove a leading 'prefix:' string if present."""
    if value.startswith(prefix + ':'):
        return value[len(prefix) + 1:]
    return value


# ---------------------------------------------------------------------------
# Translation table builders
# ---------------------------------------------------------------------------

def _build_uniprot_map(handle) -> dict[str, list[str]]:
    """
    Build gene-to-UniProt mapping from /conv/uniprot/{organism}.

    Returns dict mapping bare gene ID (e.g. 'mmu:12385') to list of
    UniProt accessions.
    """
    result: dict[str, list[str]] = {}
    for gene_id, uniprot_raw in _iter_tsv(handle):
        uniprot = _strip_prefix(uniprot_raw, 'up')
        result.setdefault(gene_id, []).append(uniprot)
    return result


def _build_ec_map(handle) -> dict[str, list[str]]:
    """
    Build gene-to-EC mapping from /link/enzyme/{organism}.

    Returns dict mapping bare gene ID to list of EC number strings.
    """
    result: dict[str, list[str]] = {}
    for gene_id, ec_raw in _iter_tsv(handle):
        result.setdefault(gene_id, []).append(ec_raw)
    return result


def _build_rxn_map(handle) -> dict[str, list[str]]:
    """
    Build EC-to-reaction mapping from /link/reaction/enzyme.

    Returns dict mapping EC string to list of bare reaction IDs.
    """
    result: dict[str, list[str]] = {}
    for ec_id, rxn_raw in _iter_tsv(handle):
        rxn_id = _strip_prefix(rxn_raw, 'rn')
        result.setdefault(ec_id, []).append(rxn_id)
    return result


def _build_compound_map(
    list_handle,
    chebi_handle,
    pubchem_handle=None,
) -> dict[str, dict[str, Any]]:
    """
    Build compound name-to-identifiers mapping.

    Combines /list/compound (names) and /conv/chebi/compound (ChEBI IDs).

    Returns dict mapping each compound name to::

        {
            'kegg_id': 'cpd:C00001',
            'chebi':   'chebi:15377',
        }
    """
    # cpd:CXXXXX -> external IDs
    chebi_lookup: dict[str, str] = {}
    for cpd_raw, chebi_raw in _iter_tsv(chebi_handle):
        cpd_id = _strip_prefix(cpd_raw, 'cpd')
        chebi_id = chebi_raw  # keep full "chebi:XXXXX" form
        chebi_lookup[cpd_id] = chebi_id

    pubchem_lookup: dict[str, str] = {}
    for cpd_raw, pubchem_raw in _iter_tsv(pubchem_handle):
        cpd_id = _strip_prefix(cpd_raw, 'cpd')
        pubchem_lookup[cpd_id] = pubchem_raw

    # cpd:CXXXXX -> first (canonical) name
    result: dict[str, dict[str, Any]] = {}
    for cpd_raw, names_raw in _iter_tsv(list_handle):
        cpd_id = _strip_prefix(cpd_raw, 'cpd')
        for name in names_raw.split(';'):
            name = name.strip()
            if not name:
                continue
            if name not in result:
                result[name] = {
                    'kegg_id': f'cpd:{cpd_id}',
                    'chebi':   chebi_lookup.get(cpd_id),
                    'pubchem': pubchem_lookup.get(cpd_id),
                }
    return result


def _build_compound_id_map(
    list_handle,
    chebi_handle,
    pubchem_handle=None,
) -> dict[str, dict[str, Any]]:
    """Build compound ID-to-metadata mapping from KEGG list/conv endpoints."""
    chebi_lookup: dict[str, str] = {}
    for cpd_raw, chebi_raw in _iter_tsv(chebi_handle):
        chebi_lookup[_strip_prefix(cpd_raw, 'cpd')] = chebi_raw

    pubchem_lookup: dict[str, str] = {}
    for cpd_raw, pubchem_raw in _iter_tsv(pubchem_handle):
        pubchem_lookup[_strip_prefix(cpd_raw, 'cpd')] = pubchem_raw

    result: dict[str, dict[str, Any]] = {}
    for cpd_raw, names_raw in _iter_tsv(list_handle):
        cpd_id = _strip_prefix(cpd_raw, 'cpd')
        names = [name.strip() for name in names_raw.split(';') if name.strip()]
        result[cpd_id] = {
            'name': names[0] if names else '',
            'synonyms': names[1:],
            'chebi': chebi_lookup.get(cpd_id),
            'pubchem': pubchem_lookup.get(cpd_id),
        }
    return result


def _build_rxn_to_uniprot(
    uniprot_map: dict[str, list[str]],
    ec_map: dict[str, list[str]],
    rxn_map: dict[str, list[str]],
) -> dict[str, list[str]]:
    """
    Chain gene→EC→reaction to produce reaction→UniProt mapping.

    Returns dict mapping bare reaction ID to sorted list of UniProt accessions.
    """
    # EC → UniProt (via gene)
    ec_to_uniprot: dict[str, set[str]] = {}
    for gene_id, ec_list in ec_map.items():
        uniprots = uniprot_map.get(gene_id, [])
        for ec in ec_list:
            ec_to_uniprot.setdefault(ec, set()).update(uniprots)

    # Reaction → UniProt
    rxn_to_uniprot: dict[str, list[str]] = {}
    for ec, rxn_list in rxn_map.items():
        uniprots = ec_to_uniprot.get(ec, set())
        for rxn_id in rxn_list:
            rxn_to_uniprot.setdefault(rxn_id, set()).update(uniprots)  # type: ignore[arg-type]

    return {rxn_id: sorted(ups) for rxn_id, ups in rxn_to_uniprot.items()}


# ---------------------------------------------------------------------------
# Equation parser
# ---------------------------------------------------------------------------

def _strip_stoich(name: str) -> str:
    """Remove leading stoichiometric coefficients (2, n, (n+1), > 2, etc.)."""
    return _STOICH_RE.sub('', name).strip()


def _parse_equation(description: str) -> tuple[list[str], list[str]]:
    """
    Extract substrate and product name lists from a KEGG reaction description.

    KEGG descriptions follow the pattern::

        "Enzyme name; possibly more names; A + B <=> C + D"

    The equation is always after the last semicolon.  A greedy ``rsplit``
    on ``;`` correctly handles multi-semicolon descriptions without the
    non-greedy regex bug present in earlier versions of this pipeline.

    Returns
    -------
    tuple[list[str], list[str]]
        (substrates, products) — each a list of cleaned compound name strings.
    """
    equation = description.rsplit(';', 1)[-1].strip()

    match = _ARROW_RE.search(equation)
    if not match:
        return [], []

    sub_str  = equation[:match.start()].strip()
    prod_str = equation[match.end():].strip()

    substrates = [_strip_stoich(s.strip()) for s in sub_str.split(' + ') if s.strip()]
    products   = [_strip_stoich(p.strip()) for p in prod_str.split(' + ') if p.strip()]

    return substrates, products


def _parse_equation_part(part: str) -> dict[str, str] | None:
    """Parse one KEGG equation participant token such as ``2 C00031``."""
    match = _EQUATION_PART_RE.match(part.strip())
    if not match:
        return None
    return {
        'kegg_id': match.group('compound'),
        'stoichiometry': match.group('stoich') or '1',
    }


def _parse_kegg_equation(
    equation: str,
) -> tuple[str, list[dict[str, str]], list[dict[str, str]]]:
    """Parse a KEGG equation into direction, reactants, and products."""
    match = _ARROW_RE.search(equation)
    if not match:
        return '', [], []

    arrow = match.group(0)
    left = equation[:match.start()].strip()
    right = equation[match.end():].strip()

    left_parts = [
        parsed
        for token in left.split(' + ')
        if token.strip()
        for parsed in [_parse_equation_part(token)]
        if parsed
    ]
    right_parts = [
        parsed
        for token in right.split(' + ')
        if token.strip()
        for parsed in [_parse_equation_part(token)]
        if parsed
    ]

    if arrow == '<=>':
        return 'REVERSIBLE', left_parts, right_parts
    if arrow == '=>':
        return 'LEFT-TO-RIGHT', left_parts, right_parts
    if arrow == '<=':
        return 'RIGHT-TO-LEFT', right_parts, left_parts
    return '', left_parts, right_parts


def _resolve_compound(
    name: str,
    compound_map: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Look up a compound name and return its identifiers dict."""
    entry = compound_map.get(name)
    if entry:
        return {'name': name, 'kegg_id': entry['kegg_id'], 'chebi': entry['chebi']}
    return {'name': name, 'kegg_id': None, 'chebi': None}


def _iter_kegg_flat_entries(handle) -> Generator[dict[str, list[str]], None, None]:
    """Yield fixed-field KEGG flat-file entries as ``field -> lines`` dicts."""
    if not handle:
        return
    if hasattr(handle, 'seek'):
        handle.seek(0)

    fields: dict[str, list[str]] = {}
    current: str | None = None

    for raw_line in handle:
        if isinstance(raw_line, bytes):
            raw_line = raw_line.decode('utf-8', 'ignore')
        line = raw_line.rstrip('\n')
        if line == '///':
            if fields:
                yield fields
            fields = {}
            current = None
            continue
        if not line.strip():
            continue

        key = line[:12].strip()
        value = line[12:].strip()
        if key:
            current = key
            fields.setdefault(key, []).append(value)
        elif current:
            fields.setdefault(current, []).append(value)

    if fields:
        yield fields


def _entry_id(entry: dict[str, list[str]]) -> str:
    value = (entry.get('ENTRY') or [''])[0].split()
    return value[0] if value else ''


def _entry_values(entry: dict[str, list[str]], field: str) -> list[str]:
    return [value.strip() for value in entry.get(field, []) if value.strip()]


def _parse_names(entry: dict[str, list[str]]) -> list[str]:
    names: list[str] = []
    for line in _entry_values(entry, 'NAME'):
        names.extend(name.strip() for name in line.split(';') if name.strip())
    return names


def _parse_whitespace_values(entry: dict[str, list[str]], field: str) -> list[str]:
    values: list[str] = []
    for line in _entry_values(entry, field):
        values.extend(token.strip() for token in line.split() if token.strip())
    return values


def _parse_orthology(entry: dict[str, list[str]]) -> list[str]:
    ko_ids: list[str] = []
    for line in _entry_values(entry, 'ORTHOLOGY'):
        token = line.split(None, 1)[0]
        if token.startswith('K'):
            ko_ids.append(token)
    return ko_ids


def _parse_rclass(entry: dict[str, list[str]]) -> list[str]:
    rclasses: list[str] = []
    for line in _entry_values(entry, 'RCLASS'):
        token = line.split(None, 1)[0]
        if token.startswith('RC'):
            rclasses.append(token)
    return rclasses


def _parse_dblinks(entry: dict[str, list[str]]) -> dict[str, list[str]]:
    links: dict[str, list[str]] = {}
    active_db: str | None = None
    for line in _entry_values(entry, 'DBLINKS'):
        if ':' in line:
            db, values = line.split(':', 1)
            active_db = db.strip()
            links.setdefault(active_db, []).extend(values.split())
        elif active_db:
            links.setdefault(active_db, []).extend(line.split())
    return links


def _parse_pathways(entry: dict[str, list[str]]) -> list[dict[str, str]]:
    pathways: list[dict[str, str]] = []
    for line in _entry_values(entry, 'PATHWAY'):
        parts = line.split(None, 1)
        if not parts:
            continue
        pathway_id = parts[0].strip()
        pathway_name = parts[1].strip() if len(parts) > 1 else ''
        pathways.append({
            'pathway_id': pathway_id,
            'pathway_name': pathway_name,
        })
    return pathways


def _empty_pathway_record(
    pathway_id: str,
    pathway_name: str = '',
    taxon_id: str = '',
) -> dict[str, Any]:
    return {
        'pathway_id': pathway_id,
        'pathway_name': pathway_name,
        'reaction_ids': set(),
        'taxon_id': taxon_id,
        'protein_members': [],
        'small_molecule_members': [],
        'ortholog_members': [],
        'linked_pathway_members': [],
    }


def _read_handle_text(handle) -> str:
    if hasattr(handle, 'seek'):
        handle.seek(0)
    if hasattr(handle, 'read'):
        content = handle.read()
        return content.decode('utf-8') if isinstance(content, bytes) else str(content)
    return ''.join(
        chunk.decode('utf-8') if isinstance(chunk, bytes) else str(chunk)
        for chunk in handle
    )


def _strip_kegg_object_prefix(value: str) -> str:
    if ':' in value:
        return value.split(':', 1)[1]
    return value


def _graphics_name(entry_element: ET.Element) -> str:
    graphics = entry_element.find('graphics')
    return graphics.get('name', '').strip() if graphics is not None else ''


def _split_graphics_names(graphics_name: str) -> list[str]:
    return [
        name.strip()
        for name in graphics_name.replace('...', '').split(',')
        if name.strip()
    ]


def _kgml_pathway_id(root: ET.Element, fallback_id: str) -> str:
    name = root.get('name', '').strip()
    if name.startswith('path:'):
        return name[5:]
    return fallback_id


def _append_unique_member(
    members: list[dict[str, str]],
    member: dict[str, str],
    key_fields: tuple[str, ...],
) -> None:
    key = tuple(member.get(field, '') for field in key_fields)
    if not any(tuple(existing.get(field, '') for field in key_fields) == key for existing in members):
        members.append(member)


def _parse_kgml_documents(
    kgml_handles,
    *,
    compound_id_map: dict[str, dict[str, Any]],
    uniprot_map: dict[str, list[str]],
    taxon_id: str,
) -> dict[str, dict[str, Any]]:
    """Parse KGML documents into rich pathway records."""
    if not kgml_handles:
        return {}

    if isinstance(kgml_handles, dict):
        items = kgml_handles.items()
    else:
        items = [('unknown', kgml_handles)]

    pathways: dict[str, dict[str, Any]] = {}

    for fallback_id, handle in items:
        text = _read_handle_text(handle).strip()
        if not text:
            continue
        try:
            root = ET.fromstring(text)
        except ET.ParseError:
            continue

        pathway_id = _kgml_pathway_id(root, fallback_id)
        pathway = pathways.setdefault(
            pathway_id,
            _empty_pathway_record(
                pathway_id,
                root.get('title', '').strip(),
                taxon_id,
            ),
        )
        if not pathway['pathway_name']:
            pathway['pathway_name'] = root.get('title', '').strip()
        pathway['taxon_id'] = taxon_id

        entry_by_id: dict[str, ET.Element] = {}
        for entry in root.findall('entry'):
            entry_id = entry.get('id', '').strip()
            if entry_id:
                entry_by_id[entry_id] = entry

            entry_type = entry.get('type', '').strip()
            names = [name.strip() for name in entry.get('name', '').split() if name.strip()]
            reaction_ids = [
                _strip_kegg_object_prefix(reaction)
                for reaction in entry.get('reaction', '').split()
                if reaction.strip()
            ]
            graphics_name = _graphics_name(entry)
            graphics_names = _split_graphics_names(graphics_name)

            for reaction_id in reaction_ids:
                pathway['reaction_ids'].add(reaction_id)

            if entry_type == 'gene':
                for index, gene_id in enumerate(names):
                    if ':' not in gene_id:
                        continue
                    uniprots = uniprot_map.get(gene_id) or ['']
                    gene_name = (
                        graphics_names[index]
                        if index < len(graphics_names)
                        else graphics_name
                    )
                    for uniprot in uniprots:
                        _append_unique_member(
                            pathway['protein_members'],
                            {
                                'kegg_id': gene_id,
                                'entrez_id': _strip_kegg_object_prefix(gene_id),
                                'uniprot_id': uniprot,
                                'name': gene_name,
                                'reaction_ids': ';'.join(reaction_ids),
                            },
                            ('kegg_id', 'uniprot_id'),
                        )

            elif entry_type in {'compound', 'glycan', 'drug'}:
                for name in names:
                    member_id = _strip_kegg_object_prefix(name)
                    metadata = compound_id_map.get(member_id, {})
                    _append_unique_member(
                        pathway['small_molecule_members'],
                        {
                            'kegg_id': member_id,
                            'name': str(metadata.get('name') or graphics_name or member_id),
                            'chebi_id': str(metadata.get('chebi') or ''),
                            'pubchem_id': str(metadata.get('pubchem') or ''),
                        },
                        ('kegg_id',),
                    )

            elif entry_type == 'ortholog':
                for name in names:
                    member_id = name if ':' in name else f'ko:{name}'
                    _append_unique_member(
                        pathway['ortholog_members'],
                        {
                            'kegg_id': member_id,
                            'name': graphics_name or member_id,
                        },
                        ('kegg_id',),
                    )

            elif entry_type == 'map':
                for name in names:
                    linked_id = _strip_kegg_object_prefix(name)
                    _append_unique_member(
                        pathway['linked_pathway_members'],
                        {
                            'pathway_id': linked_id,
                            'pathway_name': graphics_name or linked_id,
                        },
                        ('pathway_id',),
                    )

        for reaction in root.findall('reaction'):
            for reaction_id in reaction.get('name', '').split():
                pathway['reaction_ids'].add(_strip_kegg_object_prefix(reaction_id))

    return pathways


def _join_member_field(
    members: list[dict[str, str]],
    field: str,
) -> str:
    return '||'.join(member.get(field, '') for member in members)


def _pathway_record_to_row(pathway: dict[str, Any]) -> dict[str, str]:
    protein_members = pathway.get('protein_members', [])
    small_molecule_members = pathway.get('small_molecule_members', [])
    ortholog_members = pathway.get('ortholog_members', [])
    linked_pathway_members = pathway.get('linked_pathway_members', [])

    return {
        'pathway_id': pathway['pathway_id'],
        'pathway_name': pathway.get('pathway_name', ''),
        'reaction_ids': ';'.join(sorted(pathway.get('reaction_ids', set()))),
        'taxon_id': pathway.get('taxon_id', ''),
        'protein_member_kegg_ids': _join_member_field(protein_members, 'kegg_id'),
        'protein_member_entrez_ids': _join_member_field(protein_members, 'entrez_id'),
        'protein_member_uniprot_ids': _join_member_field(protein_members, 'uniprot_id'),
        'protein_member_names': _join_member_field(protein_members, 'name'),
        'protein_member_reaction_ids': _join_member_field(protein_members, 'reaction_ids'),
        'small_molecule_member_kegg_ids': _join_member_field(small_molecule_members, 'kegg_id'),
        'small_molecule_member_names': _join_member_field(small_molecule_members, 'name'),
        'small_molecule_member_chebi_ids': _join_member_field(small_molecule_members, 'chebi_id'),
        'small_molecule_member_pubchem_ids': _join_member_field(small_molecule_members, 'pubchem_id'),
        'ortholog_member_kegg_ids': _join_member_field(ortholog_members, 'kegg_id'),
        'ortholog_member_names': _join_member_field(ortholog_members, 'name'),
        'linked_pathway_ids': _join_member_field(linked_pathway_members, 'pathway_id'),
        'linked_pathway_names': _join_member_field(linked_pathway_members, 'pathway_name'),
    }


def _participant_fields(
    participants: list[dict[str, str]],
    compound_id_map: dict[str, dict[str, Any]],
) -> dict[str, str]:
    kegg_ids: list[str] = []
    names: list[str] = []
    chebis: list[str] = []
    pubchems: list[str] = []
    stoichiometries: list[str] = []

    for participant in participants:
        cpd_id = participant['kegg_id']
        metadata = compound_id_map.get(cpd_id, {})
        kegg_ids.append(cpd_id)
        names.append(str(metadata.get('name') or ''))
        chebis.append(str(metadata.get('chebi') or ''))
        pubchems.append(str(metadata.get('pubchem') or ''))
        stoichiometries.append(participant['stoichiometry'])

    return {
        'kegg_id': '||'.join(kegg_ids),
        'name': '||'.join(names),
        'chebi': '||'.join(chebis),
        'pubchem': '||'.join(pubchems),
        'stoichiometry': '||'.join(stoichiometries),
    }


# ---------------------------------------------------------------------------
# Record builder
# ---------------------------------------------------------------------------

def _build_records(
    handles: dict[str, Any],
    organism: str = 'mmu',
) -> dict[str, list[dict[str, Any]]]:
    """
    Parse all KEGG bulk-download handles and return a dict of record lists.

    Expected keys in ``handles``
    ----------------------------
    ``'conv_uniprot'``
        Response from ``/conv/uniprot/{organism}``
    ``'link_enzyme'``
        Response from ``/link/enzyme/{organism}``
    ``'link_reaction'``
        Response from ``/link/reaction/enzyme``
    ``'list_compound'``
        Response from ``/list/compound``
    ``'conv_chebi'``
        Response from ``/conv/chebi/compound``
    ``'conv_pubchem'``
        Response from ``/conv/pubchem/compound``
    ``'list_reaction'``
        Response from ``/list/reaction``
    ``'get_reaction'``
        Combined response from batched ``/get/<reaction_ids>``
    ``'kgml_pathways'``
        Mapping of KEGG pathway ID to KGML file handle.
    """
    uniprot_map  = _build_uniprot_map(handles.get('conv_uniprot'))
    ec_map       = _build_ec_map(handles.get('link_enzyme'))
    rxn_map      = _build_rxn_map(handles.get('link_reaction'))
    compound_id_map = _build_compound_id_map(
        handles.get('list_compound'),
        handles.get('conv_chebi'),
        handles.get('conv_pubchem'),
    )
    rxn_to_up    = _build_rxn_to_uniprot(uniprot_map, ec_map, rxn_map)

    reactions: list[dict[str, Any]] = []
    pathways_by_id: dict[str, dict[str, Any]] = {}
    taxon_id = _ORGANISM_TAXA.get(organism, '')

    for entry in _iter_kegg_flat_entries(handles.get('get_reaction')):
        rxn_id = _entry_id(entry)
        if not rxn_id:
            continue

        names = _parse_names(entry)
        definition = ' '.join(_entry_values(entry, 'DEFINITION'))
        equation = ' '.join(_entry_values(entry, 'EQUATION'))
        direction, reactants, products = _parse_kegg_equation(equation)
        reactant_fields = _participant_fields(reactants, compound_id_map)
        product_fields = _participant_fields(products, compound_id_map)
        dblinks = _parse_dblinks(entry)
        pathways = _parse_pathways(entry)
        uniprot_ids = rxn_to_up.get(rxn_id, [])

        for pathway in pathways:
            pathway_id = pathway['pathway_id']
            pathway_record = pathways_by_id.setdefault(
                pathway_id,
                _empty_pathway_record(pathway_id, pathway['pathway_name'], taxon_id),
            )
            if not pathway_record['pathway_name'] and pathway['pathway_name']:
                pathway_record['pathway_name'] = pathway['pathway_name']
            pathway_record['reaction_ids'].add(rxn_id)

        reactions.append({
            'reaction_id': rxn_id,
            'reaction_name': ';'.join(names),
            'reaction_definition': definition,
            'reaction_equation': equation,
            'conversion_direction': direction,
            'ec_numbers': ';'.join(_parse_whitespace_values(entry, 'ENZYME')),
            'ko_ids': ';'.join(f'ko:{ko_id}' for ko_id in _parse_orthology(entry)),
            'rclass_ids': ';'.join(f'rc:{rclass_id}' for rclass_id in _parse_rclass(entry)),
            'rhea_ids': ';'.join(dblinks.get('RHEA', [])),
            'pathway_ids': ';'.join(pathway['pathway_id'] for pathway in pathways),
            'pathway_names': ';'.join(pathway['pathway_name'] for pathway in pathways),
            'reactant_kegg_id': reactant_fields['kegg_id'],
            'reactant_name': reactant_fields['name'],
            'reactant_chebi': reactant_fields['chebi'],
            'reactant_pubchem': reactant_fields['pubchem'],
            'reactant_stoichiometry': reactant_fields['stoichiometry'],
            'product_kegg_id': product_fields['kegg_id'],
            'product_name': product_fields['name'],
            'product_chebi': product_fields['chebi'],
            'product_pubchem': product_fields['pubchem'],
            'product_stoichiometry': product_fields['stoichiometry'],
            'uniprot_ids': ';'.join(uniprot_ids),
            'taxon_id': taxon_id,
        })

    if handles.get('get_reaction'):
        kgml_pathways = _parse_kgml_documents(
            handles.get('kgml_pathways'),
            compound_id_map=compound_id_map,
            uniprot_map=uniprot_map,
            taxon_id=taxon_id,
        )
        for pathway_id, kgml_pathway in kgml_pathways.items():
            pathway_record = pathways_by_id.setdefault(pathway_id, kgml_pathway)
            if not pathway_record['pathway_name'] and kgml_pathway.get('pathway_name'):
                pathway_record['pathway_name'] = kgml_pathway['pathway_name']
            pathway_record['taxon_id'] = kgml_pathway.get('taxon_id', pathway_record.get('taxon_id', ''))
            pathway_record['reaction_ids'].update(kgml_pathway.get('reaction_ids', set()))
            for field in (
                'protein_members',
                'small_molecule_members',
                'ortholog_members',
                'linked_pathway_members',
            ):
                for member in kgml_pathway.get(field, []):
                    key_fields = (
                        ('kegg_id', 'uniprot_id')
                        if field == 'protein_members'
                        else ('pathway_id',)
                        if field == 'linked_pathway_members'
                        else ('kegg_id',)
                    )
                    _append_unique_member(pathway_record[field], member, key_fields)

        pathways = [
            _pathway_record_to_row(pathway)
            for pathway in sorted(pathways_by_id.values(), key=lambda item: item['pathway_id'])
        ]
        return {'reactions': reactions, 'pathways': pathways}

    compound_map = _build_compound_map(
        handles.get('list_compound'),
        handles.get('conv_chebi'),
        handles.get('conv_pubchem'),
    )

    for rxn_raw, description in _iter_tsv(handles.get('list_reaction')):
        rxn_id = _strip_prefix(rxn_raw, 'rn')

        uniprot_ids = rxn_to_up.get(rxn_id)
        if not uniprot_ids:
            continue

        substrates, products = _parse_equation(description)

        reactions.append({
            'Reaction':  rxn_id,
            'UniProt':   ';'.join(uniprot_ids),
            'Substrate': [_resolve_compound(s, compound_map) for s in substrates],
            'Product':   [_resolve_compound(p, compound_map) for p in products],
        })

    return {'reactions': reactions}
