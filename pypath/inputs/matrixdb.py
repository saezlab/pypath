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

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.uniprot_db as uniprot_db


def matrixdb_interactions(organism = 9606):
    """
    Downloads and processes MatrixDB interactions in MITAB format.
    Returns a list of MatrixdbInteraction namedtuples.
    """
    
    MatrixdbInteraction = collections.namedtuple(
        'MatrixdbInteraction',
        [
            'id_a',
            'id_b', 
            'alt_ids_a',
            'alt_ids_b',
            'aliases_a',
            'aliases_b',
            'detection_method',
            'author',
            'pmids',
            'taxid_a',
            'taxid_b',
            'interaction_type',
            'source_db',
            'interaction_id',
            'confidence'
        ]
    )

    url = urls.urls['matrixdb']['url']
    
    # Handle ZIP file format for new MatrixDB version
    if url.endswith('.zip'):
        c = curl.Curl(url, silent = False, large = True)
        zipfile = curl.FileOpener(
            c.fileobj.name,
            compr = 'zip',
            files_needed = ['matrixdb_CORE.tab'],
            large = True,
        )
        data = zipfile.result
        fileobj = data['matrixdb_CORE.tab']
        content = fileobj.read()
        if isinstance(content, bytes):
            content = content.decode('utf-8')
        f = content.split('\n')
        c.close()
    else:
        c = curl.Curl(url, silent = False, large = True)
        f = c.result
    
    result = []
    
    # Skip header line
    header_skipped = False
    
    for line in f:
        if not header_skipped:
            header_skipped = True
            continue
            
        # Handle both file iterator and list from ZIP
        if isinstance(line, str):
            line = line.replace('\n', '').replace('\r', '')
        else:
            line = line.decode('utf-8') if isinstance(line, bytes) else str(line)
            line = line.replace('\n', '').replace('\r', '')
        
        if not line.strip():
            continue
            
        fields = line.split('\t')
        
        # Skip if not enough fields for basic MITAB format
        if len(fields) < 15:
            continue
            
        # Extract organism info (taxid fields are at positions 9 and 10)
        try:
            taxid_a = fields[9] if fields[9] != '-' else None
            taxid_b = fields[10] if fields[10] != '-' else None
            
            # Extract numeric taxid for organism filtering
            spec_a = 0
            spec_b = 0
            if taxid_a and 'taxid:' in taxid_a:
                spec_a = int(taxid_a.split(':')[1].split('(')[0])
            if taxid_b and 'taxid:' in taxid_b:
                spec_b = int(taxid_b.split(':')[1].split('(')[0])
        except (IndexError, ValueError):
            spec_a = spec_b = 0
            taxid_a = taxid_b = None

        # Apply organism filter
        if organism is not None and not (spec_a == organism or spec_b == organism):
            continue
            
        # Create MITAB interaction record
        interaction = MatrixdbInteraction(
            id_a=fields[0],
            id_b=fields[1],
            alt_ids_a=fields[2] if len(fields) > 2 else '-',
            alt_ids_b=fields[3] if len(fields) > 3 else '-',
            aliases_a=fields[4] if len(fields) > 4 else '-',
            aliases_b=fields[5] if len(fields) > 5 else '-',
            detection_method=fields[6] if len(fields) > 6 else '-',
            author=fields[7] if len(fields) > 7 else '-',
            pmids=fields[8] if len(fields) > 8 else '-',
            taxid_a=taxid_a,
            taxid_b=taxid_b,
            interaction_type=fields[11] if len(fields) > 11 else '-',
            source_db=fields[12] if len(fields) > 12 else '-',
            interaction_id=fields[13] if len(fields) > 13 else '-',
            confidence=fields[14] if len(fields) > 14 else '-'
        )
        
        result.append(interaction)

    # Only close if it's a file object (not a list from ZIP)
    if hasattr(f, 'close'):
        f.close()

    return result


def _matrixdb_protein_list(category, organism = 9606):
    """
    Returns a set of proteins annotated by MatrixDB.

    :arg str category:
        The protein annotation category. Possible values: `ecm`, `membrane`
        or `secreted`.
    """

    url = urls.urls['matrixdb']['%s_proteins' % category]
    c = curl.Curl(url, silent = False, large = True)

    # Handle missing files (only ecm_proteins.csv exists currently)
    if c.result is None or c.status != 0:
        import pypath.share.session as session
        _logger = session.Logger(name = 'matrixdb_input')
        _logger._log(
            f'MatrixDB {category} proteins file not available. '
            f'URL: {url}, Status: {c.status}'
        )
        return set()

    proteins = set()

    # header row
    _ = next(c.result)

    for l in c.result:
        if not l:
            continue

        # CSV format, not tab-delimited
        fields = l.strip().replace('"', '').split(',')
        if len(fields) > 0 and fields[0]:
            proteins.add(fields[0])

    proteins = mapping.map_names(proteins, 'uniprot', 'uniprot')

    if organism:

        uniprots = uniprot_db.all_uniprots(
            organism = organism,
            swissprot = True,
        )
        proteins = proteins & set(uniprots)

    return proteins


def matrixdb_membrane_proteins(organism = 9606):
    """
    Returns a set of membrane protein UniProt IDs retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('membrane', organism = organism)


def matrixdb_secreted_proteins(organism = 9606):
    """
    Returns a set of secreted protein UniProt IDs retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('secreted', organism = organism)


def matrixdb_ecm_proteins(organism = 9606):
    """
    Returns a set of ECM (extracellular matrix) protein UniProt IDs
    retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('ecm', organism = organism)


def matrixdb_annotations(organism = 9606):

    MatrixdbAnnotation = collections.namedtuple(
        'MatrixdbAnnotation',
        ('mainclass',),
    )
    annot = collections.defaultdict(set)

    for cls in ('membrane', 'secreted', 'ecm'):
        cls_annot = MatrixdbAnnotation(mainclass = cls)

        method = globals()['matrixdb_%s_proteins' % cls]

        for uniprot in method(organism = organism):
            annot[uniprot].add(cls_annot)

    return dict(annot)
