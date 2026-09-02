"""Tests for the PMI-DB raw parser."""

from __future__ import annotations

from io import StringIO
from types import SimpleNamespace

from pypath.inputs_v2.pmidb import parser


def _opener(text: str) -> SimpleNamespace:
    return SimpleNamespace(result = StringIO(text))


def test_parser_accepts_download_control_arguments():
    rows = list(
        parser(
            _opener('C00002\tHMDB00538\tATP\turl\tformula\n'),
            header = ['kegg', 'hmdb', 'name', 'link', 'formula'],
            force_refresh = True,
        )
    )

    assert rows == [
        {
            'kegg': 'C00002',
            'hmdb': 'HMDB00538',
            'name': 'ATP',
            'link': 'url',
            'formula': 'formula',
        }
    ]
