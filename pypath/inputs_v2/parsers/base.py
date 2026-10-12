"""
Common raw data parsers for standard file formats.

These parsers handle standard formats: delimited text, JSON, JSONL, GMT, SQLite,
Parquet, R data frames, spreadsheets, the tables and links of HTML pages, XML and
PSI-MI XML. A resource needs its own parser only for a format these do not read.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable, Sequence
import csv
from html.parser import HTMLParser
import io
import itertools
import json
from pathlib import Path
import re
import sqlite3
from typing import Any
import xml.etree.ElementTree as ET

try:
    import duckdb
except ImportError:  # pragma: no cover - optional dependency
    duckdb = None


def _first_handle(opener) -> Any | None:
    """Extract the first file handle from an opener result."""
    if not opener or not opener.result:
        return None
    if isinstance(opener.result, dict):
        return next(iter(opener.result.values()), None)
    return opener.result


def content_text(content: Any, encoding: str = 'utf-8') -> str:
    """Text of a string, bytes, a file handle or a line iterator."""
    if content is None:
        return ''
    if hasattr(content, 'read'):
        content = content.read()
    if isinstance(content, bytes):
        return content.decode(encoding)
    if isinstance(content, str):
        return content
    return ''.join(
        chunk.decode(encoding) if isinstance(chunk, bytes) else str(chunk) for chunk in content
    )


def read_text(opener, encoding: str = 'utf-8') -> str:
    """The whole content of a download as text; the first member of an archive."""
    return content_text(_first_handle(opener), encoding=encoding)


def _text_lines(handle: Any, encoding: str = 'utf-8') -> Iterable[str]:
    """Lines of a text handle, of whole content, or of a binary handle (decoded)."""
    if isinstance(handle, bytes):
        handle = handle.decode(encoding)
    if isinstance(handle, str):
        return io.StringIO(handle.lstrip('\ufeff'), newline='')
    if isinstance(handle, io.TextIOBase):
        return handle
    return (line.decode(encoding) if isinstance(line, bytes) else line for line in handle)


def _require_columns(names: Sequence[str], required: Sequence[str] | None) -> None:
    """Raise if the header lacks a column the caller relies on."""
    missing = [name for name in required or () if name not in names]
    if missing:
        raise ValueError(f'Table lacks expected columns: {", ".join(missing)}')


def _raw(
    opener,
    delimiter: str = ',',
    fieldnames: Sequence[str] | None = None,
    skip_rows: int = 0,
    strict: bool = False,
    quotechar: str | None = '"',
    allow_repeated: bool = False,
    encoding: str = 'utf-8',
    required: Sequence[str] | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parse delimited text files.

    Args:
        opener: File opener from download_and_open
        delimiter: Field delimiter (default: ',')
        fieldnames: Column names, for a file without a header row, or in place
            of a header row that ``skip_rows`` skips
        skip_rows: Number of lines to skip at the start of the file
        strict: Raise ``ValueError`` for a row whose number of fields is not the
            number of columns, and skip rows that hold only blank fields
        quotechar: Character that quotes a field; None reads quotes as data
        allow_repeated: Accept a header that repeats a column name; the last
            column of that name is the one returned
        encoding: Encoding of a download opened in binary mode
        required: Column names the header must have; a missing one, or a table
            with no rows, raises ``ValueError``

    Yields:
        Dictionary for each row. Unless ``strict``, a short row has None for
        the columns it lacks and a long one has its extra fields as a list
        under the key None, as with ``csv.DictReader``.

    A header that repeats a column name raises ``ValueError``, because only one
    of those columns could be returned: name the columns with ``fieldnames``.
    """
    handle = _first_handle(opener)
    if not handle:
        return
    lines = _text_lines(handle, encoding)
    if skip_rows:
        lines = itertools.islice(lines, skip_rows, None)
    quoting = {'quotechar': quotechar} if quotechar else {'quoting': csv.QUOTE_NONE}
    rows = csv.reader(lines, delimiter=delimiter, **quoting)
    names = list(fieldnames) if fieldnames is not None else next(rows, None)
    if names is None:
        return
    header = ' '.join(names).lstrip('\ufeff ').lower()
    if header.startswith(('<!doctype html', '<html', '<head', '<body')):
        raise ValueError('Expected delimited source records, received an HTML document')
    repeated = sorted({name for name in names if name and names.count(name) > 1})
    if repeated and not allow_repeated:
        raise ValueError(
            f'Header repeats column names {repeated}; pass fieldnames, '
            'or allow_repeated=True to keep the last column of each name'
        )
    _require_columns(names, required)
    width = len(names)
    empty = required is not None
    for row in rows:
        if strict:
            if not any(field.strip() for field in row):
                continue
            if len(row) != width:
                raise ValueError(
                    f'Line {skip_rows + rows.line_num} has {len(row)} fields, expected {width}'
                )
        elif not row:
            continue
        record = dict(zip(names, row))
        if len(row) > width:
            record[None] = row[width:]
        else:
            for name in names[len(row):]:
                record[name] = None
        empty = False
        yield record
    if empty:
        raise ValueError('Table has no rows')


def iter_csv(opener, delimiter: str = ',', **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse CSV files with configurable delimiter; the options are those of ``_raw``."""
    yield from _raw(opener, delimiter=delimiter, **_kwargs)


def iter_tsv(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse TSV (tab-separated values) files; the options are those of ``_raw``."""
    yield from _raw(opener, delimiter='\t', **_kwargs)


def iter_semicolon(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse semicolon-delimited files; the options are those of ``_raw``."""
    yield from _raw(opener, delimiter=';', **_kwargs)


def iter_gmt(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Sets of a GMT (Gene Matrix Transposed) file, which has no header.

    Yields one dict per line: ``name`` (first column, as written),
    ``description`` (second column) and ``members`` (the remaining non-empty
    columns, in file order). Empty lines are skipped; a line with fewer than
    two columns is an error.
    """
    handle = _first_handle(opener)
    if not handle:
        return
    for number, line in enumerate(_text_lines(handle), start=1):
        line = line.rstrip('\r\n')
        if not line.strip():
            continue
        if line.lstrip('\ufeff ').lower().startswith(('<!doctype html', '<html')):
            raise ValueError('Expected a GMT file, received an HTML document')
        name, tab, rest = line.partition('\t')
        if not tab:
            raise ValueError(f'GMT line {number} has no description column: {line!r}')
        description, *members = rest.split('\t')
        yield {
            'name': name,
            'description': description,
            'members': [member for member in members if member.strip()],
        }


def iter_json(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    Parse JSON files.

    Yields individual records if the JSON is a list, otherwise yields the whole object.
    """
    handle = _first_handle(opener)
    if not handle:
        return
    data = json.load(handle)
    if isinstance(data, list):
        yield from data
    else:
        yield data


def iter_jsonl(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse JSONL (JSON Lines) files - one JSON object per line."""
    handle = _first_handle(opener)
    if not handle:
        return
    for line in handle:
        if line.strip():
            yield json.loads(line)


def read_bytes(opener) -> bytes:
    """The whole content of a binary download."""
    handle = _first_handle(opener)
    if handle is None:
        return b''
    # A small download can arrive as its whole content instead of a handle.
    if isinstance(handle, (bytes, bytearray)):
        return bytes(handle)
    if isinstance(handle, str):
        return handle.encode('utf-8')
    if hasattr(handle, 'read'):
        data = handle.read()
        return data if isinstance(data, bytes) else data.encode('utf-8')
    return b''.join(
        chunk if isinstance(chunk, bytes) else str(chunk).encode('utf-8') for chunk in handle
    )


def _cell(value: Any) -> Any:
    """A table cell as a plain Python value; missing (R ``NA``, NaN, empty) is None."""
    if hasattr(value, 'item'):
        value = value.item()
    if value is None or value != value:  # NaN
        return None
    if isinstance(value, str):
        value = value.strip()
        return value or None
    return value


def iter_rdata(opener, name: str | None = None, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Rows of an R data frame in an ``.rds`` or ``.rda``/``.RData`` file.

    ``name`` is the object to read from a file that holds several; a file with one
    object needs none. Missing values are None.
    """
    import warnings

    import rdata

    data = read_bytes(opener)
    if not data:
        return
    with warnings.catch_warnings():
        # An .rds file has no header saying so, and tibbles convert as data frames.
        warnings.filterwarnings('ignore', category=UserWarning, module='rdata')
        converted = rdata.conversion.convert(rdata.parser.parse_data(data))
    if isinstance(converted, dict):
        if name is None and len(converted) != 1:
            raise ValueError(f'R data file holds {sorted(converted)}; name the object to read')
        converted = converted[name] if name else next(iter(converted.values()))
    if not hasattr(converted, 'to_dict'):
        raise ValueError(f'Expected an R data frame, found {type(converted).__name__}')
    for row in converted.to_dict(orient='records'):
        yield {str(column): _cell(value) for column, value in row.items()}


def iter_excel(
    opener,
    sheet: str | int = 0,
    header_row: int = 1,
    fieldnames: Sequence[str] | None = None,
    required: Sequence[str] | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Rows of a worksheet in an ``.xlsx`` or ``.xls`` file.

    ``sheet`` is a sheet name or position; ``header_row`` is the row (from 1) that
    holds the column names, and rows above it are skipped. Empty cells are None and
    empty rows are skipped.

    ``fieldnames`` names the columns itself: for a header that repeats a name or
    spans several rows (``header_row`` is then its last row), or for a sheet without
    one (``header_row=0``). A header that repeats a name is an error otherwise,
    because only one of those columns could be returned.
    
    ``required`` lists column names the header must have; a missing one, or a
    sheet with no rows, raises ``ValueError``.
    """
    import io

    data = read_bytes(opener)
    if not data:
        return
    if data[:2] == b'PK':  # .xlsx is a zip archive
        import openpyxl

        book = openpyxl.load_workbook(io.BytesIO(data), read_only=True, data_only=True)
        table = book[sheet] if isinstance(sheet, str) else book.worksheets[sheet]
        rows = table.iter_rows(values_only=True)
    else:
        import xlrd

        book = xlrd.open_workbook(file_contents=data)
        table = book.sheet_by_name(sheet) if isinstance(sheet, str) else book.sheet_by_index(sheet)
        rows = (table.row_values(index) for index in range(table.nrows))
    rows = itertools.islice(rows, max(header_row - 1, 0), None)
    header = [str(name).strip() if name is not None else '' for name in (next(rows, ()) if header_row else ())]
    if fieldnames is not None:
        header = list(fieldnames)
    repeated = sorted({name for name in header if name and header.count(name) > 1})
    if repeated:
        raise ValueError(f'Worksheet header repeats column names {repeated}; pass fieldnames')
    _require_columns(header, required)
    empty = required is not None
    for row in rows:
        record = {name: _cell(value) for name, value in zip(header, row) if name}
        if any(value is not None for value in record.values()):
            empty = False
            yield record
    if empty:
        raise ValueError('Worksheet has no rows')


def iter_parquet(
    opener=None,
    path: Path | str | None = None,
    query: str | None = None,
    batch_size: int = 100_000,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Iterate over rows from a Parquet file or dataset using DuckDB."""
    if duckdb is None:
        raise ImportError("duckdb is required to iterate Parquet inputs.")

    if path is None:
        handle = _first_handle(opener)
        path = getattr(handle, 'name', None)

    if path is None:
        return
    path = Path(path)
    parquet_path = str(path / '**' / '*.parquet') if path.is_dir() else str(path)

    connection = duckdb.connect(':memory:')
    try:
        if query is None:
            query = "SELECT * FROM read_parquet(?)"
            cursor = connection.execute(query, [parquet_path])
        else:
            cursor = connection.execute(query, [parquet_path])

        columns = [desc[0] for desc in cursor.description]
        while rows := cursor.fetchmany(batch_size):
            for row in rows:
                yield dict(zip(columns, row))
    finally:
        connection.close()


def iter_sqlite(
    opener,
    table_name: str | None = None,
    sqlite_path: Path | None = None,
    db_rel_path: str | None = None,
    query: str | None = None,
    **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """
    Iterates over a table or query in a SQLite database.
    If the database is within an archive (e.g., .tar.gz), it extracts it
    to a local cache path first.

    Args:
        opener: The opener object for handling database connections.
        table_name: The name of the table to retrieve data from (if query not provided).
        sqlite_path: The local path to the cached SQLite database.
        db_rel_path: Optional relative path to the database file within the archive.
        query: Optional custom SQL query to execute.
        max_records: Optional maximum number of records to retrieve.

    Yields:
        Dictionary representing a row from the table or query.
    """

    if sqlite_path and not sqlite_path.exists():

        if not opener or not opener.result:

            return

        _extract_sqlite_from_opener(opener, sqlite_path, db_rel_path)

    results = _iter_table(sqlite_path, table_name, query=query)

    yield from results


def _extract_sqlite_from_opener(opener, target_path: Path, db_rel_path: str | None) -> None:
    """Extracts a SQLite database file from the opener to a local path."""
    if isinstance(opener.result, dict):
        if not db_rel_path:
            raise ValueError("Relative path for file in archive is required for extracted archive results.")

        if db_rel_path not in opener.result:
            raise FileNotFoundError(f"Could not find {db_rel_path} in archive.")

        db_handle = opener.result[db_rel_path]
    else:
        db_handle = opener.result

    print(f"Extracting SQLite database to: {target_path}")
    target_path.parent.mkdir(parents=True, exist_ok=True)

    with open(target_path, 'wb') as fp:
        while chunk := db_handle.read(1024 ** 2):
            fp.write(chunk)


def _iter_table(
    sqlite_path: Path,
    table_name: str | None = None,
    query: str | None = None,
) -> Generator[dict[str, Any], None, None]:
    """Pure generator for yielding rows from a SQLite table or query."""
    connection = None
    try:
        connection = sqlite3.connect(sqlite_path)
        connection.row_factory = sqlite3.Row
        cursor = connection.cursor()

        if not query:
            if not table_name:
                raise ValueError("Either table_name or query must be provided.")
            query = f"SELECT * FROM {table_name};"

        cursor.execute(query)

        for row in cursor:
            yield dict(row)

    except sqlite3.Error as e:
        print(f"SQLite error during iteration: {e}")
        raise
    finally:
        if connection:
            connection.close()


# ================================ HTML pages =================================


def _collapse(text: str) -> str:
    """Text with its runs of whitespace collapsed to one space."""
    return ' '.join(text.split())


class _HtmlTables(HTMLParser):
    """Collect tables in the order their ``<table>`` tags open."""

    _SKIPPED = frozenset({'script', 'style'})

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.tables: list[dict[str, Any]] = []
        # Tables that are open now, the innermost last.
        self._open: list[dict[str, Any]] = []
        self._skipping = 0

    def _cell(self) -> dict[str, Any] | None:
        return self._open[-1].get('_cell') if self._open else None

    @staticmethod
    def _end_node(cell: dict[str, Any], end_link: bool = False) -> None:
        """Close the text node being read; a link closes only at its own end."""
        link = cell['_link']
        if link is not None and not end_link:
            return
        text = ''.join(cell['_node'])
        cell['_node'] = []
        if link is not None:
            link['text'] = _collapse(text)
            cell['content'].append(link)
            cell['_link'] = None
        elif text.strip():
            cell['content'].append(text)

    def _close_cell(self, table: dict[str, Any]) -> None:
        cell = table.pop('_cell', None)
        if cell is None:
            return
        self._end_node(cell, end_link=True)
        del cell['_node'], cell['_link']
        cell['text'] = _collapse(''.join(cell.pop('_parts')))
        table['_row'].append(cell)

    def _close_row(self, table: dict[str, Any]) -> None:
        self._close_cell(table)
        row = table.pop('_row', None)
        if row:
            table['rows'].append(row)

    def handle_starttag(self, tag, attrs):
        attributes = {key: value or '' for key, value in attrs}
        cell = self._cell()
        if cell is not None:
            self._end_node(cell, end_link=tag == 'a')
        if tag in self._SKIPPED:
            self._skipping += 1
        elif tag == 'table':
            table = {'attrs': attributes, 'rows': []}
            self.tables.append(table)
            self._open.append(table)
        elif self._open:
            table = self._open[-1]
            if tag == 'tr':
                self._close_row(table)
                table['_row'] = []
            elif tag in ('td', 'th'):
                self._close_cell(table)
                table.setdefault('_row', [])
                table['_cell'] = {
                    'header': tag == 'th',
                    'attrs': attributes,
                    'links': [],
                    'content': [],
                    '_parts': [],
                    '_node': [],
                    '_link': None,
                }
            elif cell is not None:
                if tag == 'a' and attributes.get('href'):
                    cell['links'].append(attributes['href'])
                    cell['_link'] = {'href': attributes['href'], 'text': '', 'attrs': attributes}
                elif tag == 'br':
                    cell['_parts'].append(' ')

    def handle_endtag(self, tag):
        cell = self._cell()
        if cell is not None:
            self._end_node(cell, end_link=tag == 'a')
        if tag in self._SKIPPED:
            self._skipping = max(0, self._skipping - 1)
        elif not self._open:
            return
        elif tag == 'table':
            self._close_row(self._open.pop())
        elif tag == 'tr':
            self._close_row(self._open[-1])
        elif tag in ('td', 'th'):
            self._close_cell(self._open[-1])

    def handle_data(self, data):
        cell = self._cell()
        if self._skipping or cell is None:
            return
        cell['_parts'].append(data)
        cell['_node'].append(data)

    def handle_comment(self, data):
        cell = self._cell()
        if cell is not None:
            self._end_node(cell)

    def close(self):
        super().close()
        while self._open:
            self._close_row(self._open.pop())


def iter_html_tables(html: str) -> Generator[dict[str, Any], None, None]:
    """Tables of an HTML document, or of the markup of one element, in document order.

    Each table is ``{'attrs': {...}, 'rows': [[cell, ...], ...]}``. ``attrs``
    are the attributes of the ``<table>`` tag. A cell is a dictionary:

    - ``header``: True for ``<th>``, False for ``<td>``;
    - ``attrs``: the cell's attributes (``colspan``, ``title``, ...);
    - ``text``: the cell's text, with runs of whitespace collapsed;
    - ``links``: the ``href`` of every link in the cell, in order;
    - ``content``: the cell's text nodes and links in document order. A text
      node is a string as written. A link is ``{'href': ..., 'text': ...,
      'attrs': {...}}`` and its text is not repeated as a text node. This
      keeps apart pieces of a cell that ``text`` runs together.

    Spanning cells are not repeated: a row holds the cells written in its
    ``<tr>``. A table inside a cell is yielded on its own, after the table
    that contains it, and its content is not part of the outer cell. Rows
    without cells are left out, and scripts and styles are not text.
    """
    parser = _HtmlTables()
    parser.feed(html)
    parser.close()
    yield from parser.tables


class _HtmlLinks(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.links: list[dict[str, Any]] = []
        self._link: dict[str, Any] | None = None

    def _close(self) -> None:
        if self._link is not None:
            self._link['text'] = _collapse(''.join(self._link.pop('_parts')))
            self.links.append(self._link)
            self._link = None

    def handle_starttag(self, tag, attrs):
        if tag != 'a':
            return
        self._close()
        attributes = {key: value or '' for key, value in attrs}
        if attributes.get('href'):
            self._link = {
                'href': attributes['href'],
                'attrs': attributes,
                '_parts': [],
            }

    def handle_endtag(self, tag):
        if tag == 'a':
            self._close()

    def handle_data(self, data):
        if self._link is not None:
            self._link['_parts'].append(data)

    def close(self):
        super().close()
        self._close()


def iter_html_links(html: str) -> Generator[dict[str, Any], None, None]:
    """Links (``<a href=...>``) of an HTML document, in document order.

    Each link is ``{'href': ..., 'text': ..., 'attrs': {...}}``; ``text`` is
    the link's text with runs of whitespace collapsed.
    """
    parser = _HtmlLinks()
    parser.feed(html)
    parser.close()
    yield from parser.links


# ==================================== XML ====================================


def read_xml_root(opener) -> ET.Element | None:
    """Root element of an XML download; None if there is no download.

    ``Opener.result`` is a file handle, bytes, a string, an iterator of lines
    or, for an archive, a mapping of member names to one of these (the first
    member is read).
    """
    result = getattr(opener, 'result', None) if opener else None
    if isinstance(result, dict):
        result = next(iter(result.values()), None)
    if result is None:
        return None
    if hasattr(result, 'read'):
        return ET.parse(result).getroot()
    if isinstance(result, (bytes, str)):
        return ET.fromstring(result)
    parser = ET.XMLParser()
    for chunk in result:
        parser.feed(chunk)
    return parser.close()


# ================================ PSI-MI XML =================================

# The two levels write most things differently; the helpers they share take
# ``level1`` where the level 1 reader collapses whitespace inside a value and
# the level 2.5 reader only strips it.

# A PSI-MI accession, written with or without its prefix (``MI:0018``, ``0018``).
_MI_RE = re.compile(r'(?:MI:)?(\d{4})')
_PSIMI_LEVELS = {'1': '1', '2': '2.5'}


def _xml_text(element: ET.Element | None, path: str, level1: bool = False) -> str | None:
    """Text of a descendant, stripped; None if it is absent or empty."""
    if element is None:
        return None
    value = element.findtext(path) or ''
    return (_collapse(value) if level1 else value.strip()) or None


def _xml_attribute(element: ET.Element, name: str) -> str | None:
    """A stripped attribute value; None if it is absent or empty."""
    return (element.get(name) or '').strip() or None


def _psimi_names(element: ET.Element | None, level1: bool) -> dict[str, Any]:
    return {
        'label': _xml_text(element, '{*}names/{*}shortLabel', level1),
        'name': _xml_text(element, '{*}names/{*}fullName', level1),
    }


def _psimi_xrefs(element: ET.Element | None, level1: bool) -> list[dict[str, Any]]:
    """Cross-references of an element's ``<xref>``, the primary one first."""
    xref = element.find('{*}xref') if element is not None else None
    if xref is None:
        return []
    refs = [ref for ref in xref if ref.tag.endswith(('primaryRef', 'secondaryRef'))]
    if level1:
        return [
            {
                'primary': ref.tag.endswith('primaryRef'),
                'db': _xml_attribute(ref, 'db'),
                'id': _xml_attribute(ref, 'id'),
                'secondary': _xml_attribute(ref, 'secondary'),
                'version': _xml_attribute(ref, 'version'),
            }
            for ref in refs
        ]
    return [
        {
            'primary': ref.tag.endswith('primaryRef'),
            'db': (ref.get('db') or '').strip(),
            'db_ac': ref.get('dbAc'),
            'id': (ref.get('id') or '').strip(),
            'ref_type': ref.get('refType'),
            'ref_type_ac': ref.get('refTypeAc'),
        }
        for ref in refs
    ]


def _psimi_mi(xrefs: list[dict[str, Any]], level1: bool) -> str | None:
    """PSI-MI accession of a term, from its primary reference.

    Level 2.5 gives the ``id`` of the ``psi-mi`` reference as written; level 1
    files write the database as ``mi`` or ``psi-mi`` and the accession with or
    without its prefix, so it is returned as ``MI:0018``.
    """
    for ref in xrefs:
        database = (ref['db'] or '').lower()
        if not ref['primary']:
            continue
        if not level1 and database == 'psi-mi':
            return ref['id']
        if level1 and database in ('mi', 'psi-mi'):
            match = _MI_RE.fullmatch(ref['id'] or '')
            return f'MI:{match.group(1)}' if match else None
    return None


def _psimi_term(element: ET.Element | None, level1: bool) -> dict[str, Any] | None:
    """A controlled vocabulary element: its labels and its PSI-MI accession."""
    if element is None:
        return None
    xrefs = _psimi_xrefs(element, level1)
    return {
        **_psimi_names(element, level1),
        'mi': _psimi_mi(xrefs, level1),
        'xrefs': xrefs,
    }


def _psimi_organism(element: ET.Element | None, level1: bool) -> dict[str, Any] | None:
    if element is None:
        return None
    return {'taxon': _xml_attribute(element, 'ncbiTaxId'), **_psimi_names(element, level1)}


def _psimi_confidence(element: ET.Element | None) -> dict[str, Any] | None:
    if element is None:
        return None
    return {
        'unit': _xml_attribute(element, 'unit'),
        'value': _xml_attribute(element, 'value'),
    }


def _psimi_boolean(element: ET.Element, path: str) -> bool | None:
    """An ``xs:boolean`` element; None if it is absent or not a boolean."""
    return {'true': True, '1': True, 'false': False, '0': False}.get(
        (_xml_text(element, path, True) or '').lower()
    )


def _psimi_referenced(
    parent: ET.Element,
    inline_tag: str,
    reference_tag: str,
    catalogue: dict[str, dict[str, Any]],
    build,
    what: str,
) -> list[dict[str, Any]]:
    """Children given inline (expanded form) or by reference (compact form).

    A reference is the element's text (level 2.5) or its ``ref`` attribute
    (level 1).
    """
    items = []
    for child in parent:
        tag = child.tag.rpartition('}')[2]
        if tag == inline_tag:
            items.append(build(child))
        elif tag == reference_tag:
            key = (child.get('ref') or child.text or '').strip()
            if key not in catalogue:
                raise ValueError(f'PSI-MI XML refers to unknown {what} {key!r}')
            items.append(catalogue[key])
    return items


def _psimi_experiment(element: ET.Element, level1: bool) -> dict[str, Any]:
    if level1:
        host = _psimi_organism(element.find('{*}hostOrganism'), level1)
        hosts = [host] if host else []
    else:
        hosts = [
            _psimi_organism(host, level1)
            for host in element.iterfind('{*}hostOrganismList/{*}hostOrganism')
        ]
    detection, identification = (
        ('interactionDetection', 'participantDetection')
        if level1
        else ('interactionDetectionMethod', 'participantIdentificationMethod')
    )
    experiment = {
        'id': element.get('id'),
        **_psimi_names(element, level1),
        'bibliography': _psimi_xrefs(element.find('{*}bibref'), level1),
        'xrefs': _psimi_xrefs(element, level1),
        'host_organisms': hosts,
        'detection_method': _psimi_term(element.find('{*}' + detection), level1),
        'participant_identification_method': _psimi_term(
            element.find('{*}' + identification), level1
        ),
    }
    if level1:
        experiment['confidence'] = _psimi_confidence(element.find('{*}confidence'))
    return experiment


def _psimi_interactor(element: ET.Element, level1: bool) -> dict[str, Any]:
    if level1:
        return {
            'id': element.get('id'),
            **_psimi_names(element, level1),
            'xrefs': _psimi_xrefs(element, level1),
            'organism': _psimi_organism(element.find('{*}organism'), level1),
        }
    return {
        'id': element.get('id'),
        **_psimi_names(element, level1),
        'aliases': [
            {
                'value': alias.text.strip(),
                'type': alias.get('type'),
                'type_ac': alias.get('typeAc'),
            }
            for alias in element.iterfind('{*}names/{*}alias')
            if alias.text and alias.text.strip()
        ],
        'xrefs': _psimi_xrefs(element, level1),
        'type': _psimi_term(element.find('{*}interactorType'), level1),
        'organism': _psimi_organism(element.find('{*}organism'), level1),
    }


def _psimi_participant(
    element: ET.Element,
    interactors: dict[str, dict[str, Any]],
    level1: bool,
) -> dict[str, Any]:
    inline = 'proteinInteractor' if level1 else 'interactor'
    found = _psimi_referenced(
        element,
        inline,
        inline + 'Ref',
        interactors,
        lambda child: _psimi_interactor(child, level1),
        'interactor',
    )
    interactor = found[0] if found else None
    if level1:
        return {
            'interactor': interactor,
            'role': _xml_text(element, '{*}role', level1),
            'confidence': _psimi_confidence(element.find('{*}confidence')),
            'is_tagged': _psimi_boolean(element, '{*}isTaggedProtein'),
            'is_overexpressed': _psimi_boolean(element, '{*}isOverexpressedProtein'),
        }
    return {
        'id': element.get('id'),
        'interactor': interactor,
        'biological_role': _psimi_term(element.find('{*}biologicalRole'), level1),
        'experimental_roles': [
            _psimi_term(role, level1)
            for role in element.iterfind('{*}experimentalRoleList/{*}experimentalRole')
        ],
    }


def iter_psimi_xml(
    xml: bytes | str,
    level: str | None = None,
) -> Generator[dict[str, Any], None, None]:
    """Interactions of a PSI-MI XML document of level 1 or level 2.5.

    Yields one dict per ``<interaction>``. The reader is chosen by the ``level``
    attribute of the ``<entrySet>`` (``1``, or ``2`` for PSI-MI XML 2.5); any
    other level raises ``ValueError``, as does a root that is not an
    ``<entrySet>``. ``level`` (``'1'`` or ``'2.5'``) states the level the
    caller expects and makes another level an error, because the two levels
    do not yield the same keys.

    Experiments and interactors are read where the document has them: inside
    each interaction (expanded form), or once in the entry and referred to by
    the interactions (compact form). Values are the source's own, stripped of
    surrounding whitespace; absent elements and attributes give None or an
    empty list, and nothing is inferred.

    Both levels:

    - ``entry_index``: position of the ``<entry>`` in the document, from 1.
    - ``source``: the entry's ``<source>``: ``label``, ``name`` and ``xrefs``.
    - ``label``, ``name``, ``xrefs``: of the interaction.
    - ``interaction_types``: the interaction's ``<interactionType>`` terms.
    - ``experiments``: the interaction's experiments, each with ``id``,
      ``label``, ``name``, ``bibliography`` (cross-references of ``<bibref>``),
      ``xrefs``, ``host_organisms``, ``detection_method`` and
      ``participant_identification_method``.
    - ``participants``: in document order, each with ``interactor`` (``id``,
      ``label``, ``name``, ``xrefs``, ``organism``).
    - A controlled vocabulary term is ``{'label', 'name', 'mi', 'xrefs'}`` and
      an organism ``{'taxon', 'label', 'name'}``.

    Level 2.5 only:

    - ``id`` of the interaction and ``availability``, the texts of the entry's
      ``<availability>`` elements.
    - A participant has ``id``, ``biological_role`` and ``experimental_roles``;
      its interactor has ``aliases`` (``value``, ``type``, ``type_ac``) and
      ``type``.
    - A cross-reference is ``{'primary', 'db', 'db_ac', 'id', 'ref_type',
      'ref_type_ac'}``; ``db`` and ``id`` are empty strings when absent.
    - ``mi`` of a term is the ``id`` of its ``psi-mi`` reference, as written.

    Level 1 only:

    - ``index``: position of the interaction in the document, from 1 (level 1
      interactions have no identifier of their own).
    - ``confidence`` of the interaction, of an experiment and of a participant:
      ``{'unit', 'value'}``.
    - ``source`` also has ``bibliography``, ``release`` and ``release_date``.
    - A participant has ``role`` (the text of ``<role>``: bait, prey, neutral
      or unspecified), ``is_tagged`` and ``is_overexpressed``.
    - A cross-reference is ``{'primary', 'db', 'id', 'secondary', 'version'}``.
    - ``mi`` of a term is normalised to ``MI:0018``, and whitespace inside a
      label or name is collapsed.
    """
    root = ET.fromstring(xml)
    if not root.tag.endswith('entrySet'):
        raise ValueError(f'Expected a PSI-MI XML <entrySet>, found <{root.tag}>')
    found = _PSIMI_LEVELS.get(root.get('level'))
    if found is None:
        raise ValueError(f'Unsupported PSI-MI XML level {root.get("level")!r}')
    if level is not None and level != found:
        raise ValueError(f'Expected PSI-MI XML level {level}, found level {found}')
    level1 = found == '1'
    index = 0
    for entry_index, entry in enumerate(root.iterfind('{*}entry'), 1):
        source = entry.find('{*}source')
        source_record = {
            **_psimi_names(source, level1),
            'xrefs': _psimi_xrefs(source, level1),
        }
        if level1:
            source_record.update(
                bibliography=_psimi_xrefs(
                    source.find('{*}bibref') if source is not None else None, level1
                ),
                release=source.get('release') if source is not None else None,
                release_date=source.get('releaseDate') if source is not None else None,
            )
        availability = [
            _collapse(item.text)
            for item in entry.iterfind('{*}availabilityList/{*}availability')
            if item.text and item.text.strip()
        ]
        experiments = {
            element.get('id'): _psimi_experiment(element, level1)
            for element in entry.iterfind('{*}experimentList/{*}experimentDescription')
        }
        interactors = {
            element.get('id'): _psimi_interactor(element, level1)
            for element in entry.iterfind(
                '{*}interactorList/{*}' + ('proteinInteractor' if level1 else 'interactor')
            )
        }
        for interaction in entry.iterfind('{*}interactionList/{*}interaction'):
            index += 1
            experiment_list = interaction.find('{*}experimentList')
            common = {
                **_psimi_names(interaction, level1),
                'xrefs': _psimi_xrefs(interaction, level1),
                'interaction_types': [
                    _psimi_term(element, level1)
                    for element in interaction.iterfind('{*}interactionType')
                ],
            }
            experiments_found = (
                _psimi_referenced(
                    experiment_list,
                    'experimentDescription',
                    'experimentRef',
                    experiments,
                    lambda child: _psimi_experiment(child, level1),
                    'experiment',
                )
                if experiment_list is not None
                else []
            )
            participants = [
                _psimi_participant(participant, interactors, level1)
                for participant in interaction.iterfind(
                    '{*}participantList/{*}'
                    + ('proteinParticipant' if level1 else 'participant')
                )
            ]
            if level1:
                yield {
                    'entry_index': entry_index,
                    'index': index,
                    'source': source_record,
                    **common,
                    'confidence': _psimi_confidence(interaction.find('{*}confidence')),
                    'experiments': experiments_found,
                    'participants': participants,
                }
            else:
                yield {
                    'entry_index': entry_index,
                    'source': source_record,
                    'availability': availability,
                    'id': interaction.get('id'),
                    **common,
                    'experiments': experiments_found,
                    'participants': participants,
                }
