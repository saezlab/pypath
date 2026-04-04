#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
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

from typing import Literal

import lxml.etree as etree

from pypath.inputs.hmdb.schema.common import XMLNS
import pypath.resources.urls as urls
from pypath.share.downloads import dm, _resolve_data_dir


def hmdb_xml(dataset: Literal['metabolites']) -> etree.iterparse:
    """
    Download and open the XML file of a dataset from the HMDB.

    If a pre-downloaded plain XML file exists at
    ``<cachedir>/hmdb/hmdb_{dataset}.xml`` it is used directly, bypassing
    the download.  This is the recommended workaround when Cloudflare blocks
    automated requests to hmdb.ca: download the file manually in a browser
    and place it at that path.
    """

    RECORD_TAGNAMES = {
        'metabolites': 'metabolite',
        'proteins': 'protein',
    }

    tag = f'{XMLNS}{RECORD_TAGNAMES[dataset]}'

    # Use pre-downloaded plain XML if available (Cloudflare bypass).
    xml_path = _resolve_data_dir() / 'hmdb' / f'hmdb_{dataset}.xml'

    if xml_path.exists():
        return etree.iterparse(str(xml_path), tag=tag)

    url = urls.urls['hmdb'][dataset]

    # Download file using download manager
    file_path = dm.download(
        url,
        filename=f'hmdb_{dataset}.zip',
        subfolder='hmdb',
    )

    # Open the extracted XML file from the zip
    import zipfile
    with zipfile.ZipFile(file_path) as zf:
        xml_file = zf.open(f'hmdb_{dataset}.xml')
        return etree.iterparse(xml_file, tag=tag)
