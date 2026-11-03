"""
Readers for the flat files of 3did, the database of three-dimensional
interacting domains.

3did lists, for pairs of Pfam domains (``3did_flat``) and for a Pfam domain
and a linear motif (``3did_dmi_flat``), the structures of the Protein Data
Bank in which the two are in contact. The files name the PDB entry, the chains
and the residue ranges of the two partners in the numbering of the PDB file,
and no protein: nothing in them says which UniProt protein a chain is. The
readers therefore give one record per structural instance, with the domains as
Pfam entries.

``3did_flat``, the domain-domain file::

    #=ID  <name 1>  <name 2>   (<Pfam accession 1>@Pfam  <Pfam accession 2>@Pfam)
    #=3D  <PDB id>  <chain>:<start>-<end>  <chain>:<start>-<end>  <score>  <Z-score>  [<topology>]
    <residue 1>  <residue 2>  <position 1>  <position 2>  <contact type>
    ...
    //

``3did_dmi_flat``, the domain-motif file::

    #=ID  <domain name>  (Pfam)  <motif name>  (<motif source>)
    #=PT  <motif pattern> (<date>)
    #=3D  <PDB id>  <domain chain>:<start>-<end>  <motif chain>:<start>-<end>  <peptide>  <contacts>  <topology>
    //

The residue contact lines below a structural instance are not read.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
import re
from typing import Any

__all__ = ['iter_domain_domain', 'iter_domain_motif']

_ACCESSION = re.compile(r'\s*\(?(PF\d{5})\.(\d+)@Pfam\)?\s*')
_REGION = re.compile(r'([^:\s]+):(-?\d+[A-Za-z]?)-(-?\d+[A-Za-z]?)')
_PDB_ID = re.compile(r'[0-9][a-z0-9]{3}')
_PATTERN = re.compile(r'(\S+) \(([A-Za-z]{3} \d{4})\)')


def _lines(opener) -> Iterable[str]:
    if not opener or not opener.result:
        raise ValueError('3did file is empty')
    handle = opener.result
    if isinstance(handle, dict):
        handle = next(iter(handle.values()))
    if isinstance(handle, str):
        return handle.splitlines()
    return handle


def _region(cell: str, where: str) -> dict[str, Any]:
    match = _REGION.fullmatch(cell)
    if not match:
        raise ValueError(f'3did {where}: malformed chain range {cell!r}')
    return {'chain': match.group(1), 'start': match.group(2), 'end': match.group(3)}


def _number(cell: str, where: str) -> float:
    try:
        return float(cell)
    except ValueError:
        raise ValueError(f'3did {where}: {cell!r} is not a number') from None


def iter_domain_domain(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The structural instances of ``3did_flat``, one dict each.

    Keys:

    - ``name_1``, ``name_2``: the Pfam names of the two domains (``SH2``);
    - ``pfam_1``, ``pfam_2``: their Pfam accessions without the version
      (``PF00017``), and ``pfam_version_1``, ``pfam_version_2`` the version
      numbers the file writes after the dot;
    - ``pdb_id``: the PDB entry;
    - ``region_1``, ``region_2``: ``chain``, ``start`` and ``end`` of each
      domain in the PDB file, as written (the numbering of the PDB file, which
      can carry an insertion code);
    - ``score``, ``z_score``: the InterPreTS score and Z-score (floats);
    - ``topology``: the last field of the line, ``n:m``, None where a line
      has none (24,854 of 1,038,439 lines).

    A file that is empty, a structural instance before the first domain pair,
    a line of another width, a malformed range or accession, and a file with
    no instance raise ``ValueError``.
    """
    pair: dict[str, Any] | None = None
    found = False
    for number, line in enumerate(_lines(opener), 1):
        if not line.startswith('#='):
            continue
        fields = line.rstrip('\n').split('\t')
        tag = fields[0]
        if tag == '#=ID':
            if len(fields) != 5:
                raise ValueError(f'3did line {number}: {len(fields)} fields in a domain pair')
            accessions = [_ACCESSION.fullmatch(fields[3]), _ACCESSION.fullmatch(fields[4])]
            if not all(accessions) or not fields[1] or not fields[2]:
                raise ValueError(f'3did line {number}: malformed domain pair {line!r}')
            pair = {
                'name_1': fields[1],
                'name_2': fields[2],
                'pfam_1': accessions[0].group(1),
                'pfam_version_1': accessions[0].group(2),
                'pfam_2': accessions[1].group(1),
                'pfam_version_2': accessions[1].group(2),
            }
        elif tag == '#=3D':
            if pair is None:
                raise ValueError(f'3did line {number}: structural instance without a domain pair')
            if len(fields) not in (6, 7) or not _PDB_ID.fullmatch(fields[1]):
                raise ValueError(f'3did line {number}: malformed structural instance {line!r}')
            found = True
            yield {
                **pair,
                'pdb_id': fields[1],
                'region_1': _region(fields[2], f'line {number}'),
                'region_2': _region(fields[3], f'line {number}'),
                'score': _number(fields[4], f'line {number}'),
                'z_score': _number(fields[5], f'line {number}'),
                'topology': fields[6] if len(fields) == 7 else None,
            }
        else:
            raise ValueError(f'3did line {number}: unknown record {tag!r}')
    if not found:
        raise ValueError('3did file has no structural instances')


def iter_domain_motif(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The structural instances of ``3did_dmi_flat``, one dict each.

    Keys:

    - ``domain``: the Pfam name of the domain (the file gives no accession);
    - ``motif``: 3did's name of the motif (``PCNA_C_LIG_0-0``, the name of the
      domain that binds it, a motif type and a number);
    - ``motif_source``: the database the motif comes from (``PLoS_CB_2010``);
    - ``pattern``, ``pattern_date``: the regular expression of the motif and
      the date the file gives with it (``Feb 2025``);
    - ``pdb_id``: the PDB entry;
    - ``domain_region``, ``motif_region``: ``chain``, ``start`` and ``end`` in
      the PDB file, as written;
    - ``peptide``: the sequence of the motif instance in the structure;
    - ``contacts``: the number of contextual contacts (int);
    - ``topology``: the number of the binding topology (int).

    A file that is empty, an instance before the first motif, a line of
    another width, a malformed range or number, and a file with no instance
    raise ``ValueError``.
    """
    motif: dict[str, Any] | None = None
    found = False
    for number, line in enumerate(_lines(opener), 1):
        if line.startswith('//'):
            motif = None
            continue
        if not line.startswith('#='):
            continue
        fields = line.rstrip('\n').split('\t')
        tag = fields[0]
        if tag == '#=ID':
            if len(fields) != 5 or fields[2] != '(Pfam)':
                raise ValueError(f'3did line {number}: malformed domain-motif pair {line!r}')
            motif = {
                'domain': fields[1],
                'motif': fields[3],
                'motif_source': fields[4].strip('()'),
                'pattern': None,
                'pattern_date': None,
            }
        elif tag == '#=PT':
            match = _PATTERN.fullmatch(fields[1]) if len(fields) == 2 else None
            if motif is None or not match:
                raise ValueError(f'3did line {number}: malformed motif pattern {line!r}')
            motif['pattern'], motif['pattern_date'] = match.groups()
        elif tag == '#=3D':
            if motif is None:
                raise ValueError(f'3did line {number}: structural instance without a motif')
            if len(fields) != 7 or not _PDB_ID.fullmatch(fields[1]):
                raise ValueError(f'3did line {number}: malformed structural instance {line!r}')
            found = True
            yield {
                **motif,
                'pdb_id': fields[1],
                'domain_region': _region(fields[2], f'line {number}'),
                'motif_region': _region(fields[3], f'line {number}'),
                'peptide': fields[4],
                'contacts': int(_number(fields[5], f'line {number}')),
                'topology': int(_number(fields[6], f'line {number}')),
            }
        else:
            raise ValueError(f'3did line {number}: unknown record {tag!r}')
    if not found:
        raise ValueError('3did file has no structural instances')
