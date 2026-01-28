#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
GitHub and GitLab repository utilities for standard GEM access.
"""

from __future__ import annotations

from collections.abc import Generator
from urllib.parse import quote as url_quote
import json

import pypath.share.curl as curl
import pypath.resources.urls as urls

from ._common import _log
from ._records import MetatlasGem

__all__ = [
    'metatlas_gems',
    'git_raw_file',
]


def metatlas_gems() -> Generator[MetatlasGem, None, None]:
    """
    Lists standard GEMs from the Metabolic Atlas validation index.

    Standard GEMs are community-curated models following a standardized
    format, hosted on GitHub or GitLab.

    Yields:
        MetatlasGem named tuples with repository metadata.
    """

    _log('Downloading standard GEM index from Metabolic Atlas.')

    index = _parse_gem_index()

    if not index:
        return

    for host, repos in index.items():
        for repo in repos:
            # Extract GEM name from repo path (last part)
            gem_name = repo.split('/')[-1]
            metadata = _gem_metadata(gem_name)

            if host == 'github':
                git_url = f'https://github.com/{repo}'
            elif host == 'gitlab':
                git_url = f'https://gitlab.com/{repo}'
            else:
                git_url = ''

            yield MetatlasGem(
                model_name=gem_name,
                git_host=host,
                git_repo=repo,
                git_url=git_url,
                latest_version=metadata.get('version', ''),
                commits=metadata.get('commits'),
                contributors=metadata.get('contributors'),
            )


def _parse_gem_index() -> dict:
    """
    Parses the standard-GEM validation index.json.

    Returns:
        Dictionary mapping GEM names to their repository info.
    """

    url = urls.urls['metatlas']['gems_index']
    c = curl.Curl(url, silent=False, large=False)

    if c.result is None:
        _log('Failed to download GEM index.')
        return {}

    return json.loads(c.result)


def _gem_metadata(gem_name: str) -> dict:
    """
    Retrieves metadata for a specific GEM from validation results.

    Args:
        gem_name: Name of the GEM (e.g., 'Human-GEM').

    Returns:
        Dictionary with version, commits, contributors info.
    """

    url = urls.urls['metatlas']['gem_metadata'] % gem_name
    c = curl.Curl(url, silent=True, large=False)

    if c.result is None:
        return {}

    try:
        data = json.loads(c.result)
        # Structure is {gem_name: {metadata: {...}, releases: [...]}}
        gem_data = data.get(gem_name, {})
        metadata = gem_data.get('metadata', {})

        # Extract version from releases if available
        releases = gem_data.get('releases', [])
        version = ''
        if releases:
            # Get the first release's tag name if available
            first_release = releases[0] if releases else {}
            for key in first_release.keys():
                if key not in ('standard-GEM', 'test_results'):
                    version = key
                    break

        return {
            'version': version,
            'commits': metadata.get('commits'),
            'contributors': metadata.get('contributors'),
        }
    except json.JSONDecodeError:
        return {}


def git_raw_file(
        host: str,
        repo: str,
        ref: str,
        path: str,
) -> str | None:
    """
    Downloads a raw file from GitHub or GitLab.

    Args:
        host: 'github' or 'gitlab'.
        repo: Repository path (e.g., 'SysBioChalmers/Human-GEM').
        ref: Git reference (branch, tag, or commit hash).
        path: Path to file within repository.

    Returns:
        File content as string, or None if download failed.
    """

    if host == 'github':
        url = urls.urls['metatlas']['github_raw'] % (repo, ref, path)
    elif host == 'gitlab':
        url = urls.urls['metatlas']['gitlab_raw'] % (repo, ref, path)
    else:
        _log(f'Unknown git host: {host}')
        return None

    c = curl.Curl(url, silent=False, large=True)

    if c.fileobj is None:
        _log(f'Failed to download {path} from {repo}.')
        return None

    # Read the content from the file object
    c.fileobj.seek(0)
    return c.fileobj.read()


def git_repo_tree(
        host: str,
        repo: str,
        ref: str,
) -> list[dict]:
    """
    Lists files in a repository at a specific ref.

    Args:
        host: 'github' or 'gitlab'.
        repo: Repository path (e.g., 'SysBioChalmers/Human-GEM').
        ref: Git reference (branch, tag, or commit hash).

    Returns:
        List of file dictionaries with path and type info.
    """

    if host == 'github':
        url = urls.urls['metatlas']['github_tree'] % (repo, ref)
    elif host == 'gitlab':
        encoded_repo = url_quote(repo, safe='')
        url = urls.urls['metatlas']['gitlab_tree'] % (encoded_repo, ref)
    else:
        _log(f'Unknown git host: {host}')
        return []

    c = curl.Curl(url, silent=False, large=False)

    if c.result is None:
        _log(f'Failed to get tree for {repo}@{ref}.')
        return []

    data = json.loads(c.result)

    if host == 'github':
        return data.get('tree', [])
    else:
        return data
