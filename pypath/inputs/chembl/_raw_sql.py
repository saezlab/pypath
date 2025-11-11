#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

from typing import NamedTuple, Sequence, TypeAlias
from collections.abc import Generator
from pathlib import Path
import requests
import tarfile
import sqlite3
import logging
import shutil

from pypath.share import curl, cache
from pypath.resources import urls


logger = logging.getLogger(__name__)

# Finds the latest version string in the README file
RELEASE_PREFIX = "* Release:"
class VersionInfo(NamedTuple):
    """A pair of format version and regular version."""

    fmt_version: str
    version: str
    
class VersionPathPair(NamedTuple):
    """A pair of a version and path."""

    version: str
    path: Path

class _VersionFlavorsHelper(NamedTuple):
    """A pair of format version and regular version."""

    fmt_version: str
    version: str


#: A hint for a version, which can either be an integer, string, or float (for minor versions)
VersionHint: TypeAlias = str | int | float | VersionInfo
def download_sqlite(
    version: VersionHint | None = None,
    *,
    prefix: Sequence[str] | None = None,
    return_version: bool = False,
) -> Path | VersionPathPair:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory where pypath should store the downloaded file.
    :param return_version: Should the version get returned? Turn this to true if you're
        looking up the latest version and want to reduce redundant code.

    :returns: If ``return_version`` is true, return a pair of the version and the local
        file path to the downloaded ``*.tar.gz`` file. Otherwise, just return the path.
    """
    return _download_helper(
        suffix="_sqlite.tar.gz",
        version=version,
        prefix=prefix,
        return_version=return_version,
    )


def _download_helper(
    suffix: str,
    version: VersionHint | None = None,
    prefix: Sequence[str] | None = None,
    *,
    return_version: bool,
    filename_repeats_version: bool = True,
) -> Path | VersionPathPair:
    """Ensure the latest ChEMBL file with the given suffix is downloaded.

    Args:
        suffix: The suffix of the file
        version: The version number of ChEMBL to get. If none specified, uses
            :func:`latest` to look up the latest.
        prefix: The directory that pypath should store the downloaded file.
        return_version: Should the version get returned? Turn this to true if you're
            looking up the latest version and want to reduce redundant code.
        param filename_repeats_version: True if filename contains ``chembl_<version>`` in
            the beginning. Set to false to allow downloading arbitrarily named files.

    Returns: 
        If ``return_version`` is true, return a pair of the version and the local
            file path to the downloaded file. Otherwise, just return the path.

    Raises: 
        ValueError: If file could not be downloaded
    """
    version_info = _get_version_info(version)

    if filename_repeats_version:
        filename = f"chembl_{version_info.fmt_version}{suffix}"
    else:
        filename = suffix

    # Use decided prefix to create a clear clean cache path for Curl.curl
    if prefix:
        cache_path = Path(*prefix) / version_info.version / filename
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        # full path for curl.Curl
        cache_location = str(cache_path)
    else:
        # If no prefix, let curl.Curl use its default hash-based cache
        cache_location = True

    # Creates possible URLs to download ChEMBL file
    release_url = urls.urls['chembl']['releases']
    base = f"{release_url}/chembl_{version_info.fmt_version}"
    chembl_urls = [
        f"{base}/{filename}",
        f"{base}/archived/{filename}",
    ]
    for chembl_url in chembl_urls:
        try:
            curl_obj = curl.Curl(chembl_url, large=True, silent=False, cache=cache_location, process=False)
            path = Path(curl_obj.cache_file_name)
            print(path)
        except OSError:
            continue
        if return_version:
            return VersionPathPair(version_info.version, path)
        else:
            return path

    # String of all attempted URLs for error message
    urls_fmt = "\n".join(f"   - {url}" for url in chembl_urls)

    # Get cache directory for error message
    cache_dir = cache.get_cachedir()
    raise ValueError(f"""\

[ChEMBL v{version_info.fmt_version}] could not ensure {filename}

1. It wasn't found in the pypath cache (typically at {cache_dir}):

2. It couldn't be downloaded from any of the following URLs:
{urls_fmt}
    """)

def _get_version_info(
    version: VersionHint | None,
    prefix: Sequence[str] | None = None,
) -> VersionInfo:
    """Get version info for the given version."""

    # If version is already VersionInfo, return it
    if isinstance(version, VersionInfo):
        return version

    flavor = _ensure_version_helper(version)

    return VersionInfo(flavor.fmt_version, flavor.version)

def latest(*, full: bool = False) -> str |VersionInfo:
    """Get the latest version of ChEMBL as a string.

    Returns: The latest version string of ChEMBL

    Raises: 
        ValueError: If the latest README can not be parsed
    """
    latest_readme_url = urls.urls['chembl']['latest_readme']
    res = requests.get(latest_readme_url, timeout=5)
    res.raise_for_status()
    for line_binary in res.iter_lines(decode_unicode=True):
        line: str = line_binary.decode("utf8")
        if line.startswith(RELEASE_PREFIX):
            line = line.removeprefix(RELEASE_PREFIX)
            line = line.strip()
            line = line.removeprefix("chembl_")
            if not full:
                return line
            return _get_version_info(line)

    raise ValueError("could not find latest ChEMBL version")

def _ensure_version_helper(version: VersionHint | None) -> _VersionFlavorsHelper:
    """
    Helper function to ensure a version is of type `_VersionFlavorsHelper`.

    If version is None, it is replaced with the latest version.
    If version is an int, it is converted to a string and left padded with a zero.
    If version is a str, it is stripped of leading zeros and any decimal points are replaced with underscores.
    If version is a float, it is converted to a string and any decimal points are replaced with underscores.

    Args:
        version (VersionHint | None): The version to ensure.
    Returns: 
        `_VersionFlavorsHelper`: object with the canonicalized version number
    Raises:
      TypeError: If version is not of type int, str, float or None
    """
    if version is None:
        version = latest()
    if isinstance(version, int):
        # versions 1-9 are left padded with a zero
        fmt_version = f"{version:02}"
        version = str(version)
    elif isinstance(version, str):
        # remove all leading zeros
        version = version.lstrip("0")

        # for versions 22.1 and 24.1, it's important to canonicalize the version number
        # for versions < 10 it's important to left pad with a zero
        fmt_version = version.replace(".", "_").zfill(2)
    elif isinstance(version, float):
        version = str(version)
        fmt_version = version.replace(".", "_")
    else:
        raise TypeError(f"invalid type for version: {version}")

    return _VersionFlavorsHelper(fmt_version, version)

def versions(
    *, full: bool = False
) -> list[str] | list[VersionInfo]:
    """
    Get all versions of ChEMBL. Allows the user to find all the versions available.
    The user can pick a specific version for downloading.
    
    Returns:
        list[str] | list[VersionInfo]: List of version strings or VersionInfo objects.
    """
    latest_version_info = latest(full=True)
    rv = [str(i).zfill(2) for i in range(1, int(latest_version_info.version) + 1)]
    # Side version in ChEMBL
    rv.extend(["22_1", "24_1"])
    rv = sorted(rv, reverse=True)
    if not full:
        return rv
    return [_get_version_info(version) for version in rv]

def download_extract_sqlite(
    version: VersionHint | None = None,
    *,
    prefix: Sequence[str] | None = None,
    return_version: bool = False,
    retain: bool = False,
) -> Path | VersionPathPair:
    """Ensure the latest ChEMBL SQLite dump is downloaded and extracted.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param return_version: Should the version get returned? Turn this to true if you're
        looking up the latest version and want to reduce redundant code.
    :param retain: If true, keeps the original archive.

    :returns: If ``return_version`` is true, return a pair of the version and the local
        file path to the downloaded ChEMBLSQLite database file. Otherwise, just return
        the path.

    :raises FileNotFoundError: If no database file could be found in the extracted
        directories
    """

    # Gets the version info
    version_info = _get_version_info(version)
    version_string = version_info.version
    fmt_version = version_info.fmt_version

    # Get the directory to store the data
    if prefix:
        base_dir = Path(*prefix)
    else:
        base_dir = Path.home() / "pypath" / ".chembl_data"

    # Determine the database file path
    version_dir = base_dir / version_string
    db_filename = f"chembl_{fmt_version}.db"
    rv = version_dir / db_filename

    # If the database file doesn't exist, download and extract it
    if not rv.is_file():

        # Ensure version directory exists
        version_dir.mkdir(parents=True, exist_ok=True)

        # Download the archive
        tar_path = download_sqlite(version=version_info, prefix=prefix, return_version=False)

        # Extract the database file from the archive
        with tarfile.open(tar_path, mode="r:gz") as tar_file:
            tar_info = _get_tar_info(tar_file)
            if tar_info is None:
                raise FileNotFoundError("could not find a .db file in the ChEMBL archive")
            logger.info("unarchiving %s to %s", tar_path, rv)
            with tar_file.extractfile(tar_info) as source, open(rv, "wb") as dest:
                shutil.copyfileobj(source, dest)

        # Delete the original archive if not needed and user owns file location
        if not retain and prefix:
            logger.info("deleting original archive %s", tar_path)
            tar_path.unlink()

    # Return the database file path
    if return_version:
        return VersionPathPair(version_info.version, rv)
    else:
        return rv # return the .db file path
    
def _get_tar_info(tar_file: tarfile.TarFile) -> tarfile.TarInfo | None:
    """Walk an archive and find a file with the ``.db`` extension."""
    for tar_info in tar_file:
        if tar_info.name.endswith(".db"):
            return tar_info
    return None

def query_sql_data(
    table_name: str,
    chembl_version: str = "latest"
) -> Generator[dict, None, None]:
    """
    Retrieves all data from a specified table in the ChEMBL database.

    This function uses chembl_downloader to fetch a specific version of the
    ChEMBL SQLite database, queries the specified table, and yields each
    row as a dictionary.

    Args:
        table_name (str): The name of the table to retrieve data from
                          (e.g., 'molecule_dictionary', 'activities').
        chembl_version (str): The ChEMBL version to use for reproducibility.
                              Defaults to 'latest'.

    Yields:
        dict: A dictionary representing a row from the table.
    """
    # Ensure the specified (or latest) version of the database is downloaded
    if chembl_version == "latest":
        # Let the library determine the latest version
        sqlite_path = download_sqlite()
    else:
        sqlite_path = download_sqlite(version=chembl_version)

    # Connect to the database
    with chembl_downloader.connect(sqlite_path) as connection:
        # Use sqlite3.Row as the row_factory to get dict-like rows
        connection.row_factory = sqlite3.Row
        cursor = connection.cursor()

        # Construct and execute the SQL query
        query = f"SELECT * FROM {table_name};"
        cursor.execute(query)

        # Yield each row as a dictionary
        for row in cursor:
            yield dict(row)