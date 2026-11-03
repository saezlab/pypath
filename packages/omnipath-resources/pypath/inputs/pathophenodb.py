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

import os
import pickle
import hashlib
import collections

import pypath.resources.urls as urls
import pypath.share.cache as cache
import pypath.share.session as session
import pypath.share.curl as curl

_logger = session.Logger(name = 'pathophenodb_input')
_log = _logger._log


DiseasePathogen = collections.namedtuple(
    'DiseasePathogen',
    (
        'disease_id',
        'disease',
        'pathogen_taxid',
        'pathogen',
        'evidence'
    ),
)


def disease_pathogen_interactions():
    """
    Retrieves disease pathogen relationships from PathoPhenoDb.

    Returns:
        Disease-pathogen relationships as a list of tuples.
    """

    url = urls.urls['pathophenodb']['url']
    urlmd5 = hashlib.md5(url.encode()).hexdigest()
    cache_path = os.path.join(cache.get_cachedir(), urlmd5)

    if os.path.exists(cache_path):

        try:

            _log(f'Loading from cache: `{cache_path}`.')

            with open(cache_path, 'rb') as fp:

                return pickle.load(fp)

        except:

            _log(
                f'Failed to load from `{cache_path}`, '
                'falling back to download.'
            )

    try:

        from rdflib import Graph

    except ModuleNotFoundError:

        _logger._console(
            'No module `rdflib` is available. '
            'Please install it to access PathoPhenoDB: '
            'pip install rdflib'
        )
        _log('Returning empty result!')
        return []

    # Download the RDF file
    c = curl.Curl(url, large=True, silent=False)
    
    try:
        # Parse the RDF data manually due to malformed triples
        import re
        
        rdf_data = c.result
        if hasattr(rdf_data, '__iter__') and not isinstance(rdf_data, (str, bytes)):
            rdf_data = ''.join(rdf_data)
        
        # Parse relationships manually with optimized data structures
        result = set()
        
        # Optimized data structures for multiple relationships per disease
        disease_relations = {}  # disease_id -> [rel_ids]
        relation_evidence = {}  # rel_id -> evidence_id  
        relation_pathogen = {}  # rel_id -> pathogen_id
        labels = {}  # entity_id -> label
        
        # Parse each line
        for line in rdf_data.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            # Extract triples manually
            # Disease SIO_000255 relationship
            if 'SIO_000255' in line:
                match = re.search(r'<([^>]+)>\s+<[^>]*SIO_000255>\s+<([^>]+)>', line)
                if match:
                    disease_id, rel_id = match.groups()
                    if disease_id not in disease_relations:
                        disease_relations[disease_id] = []
                    disease_relations[disease_id].append(rel_id)
            
            # Relationship RO_0002558 evidence
            elif 'RO_0002558' in line:
                match = re.search(r'<([^>]+)>\s+<[^>]*RO_0002558>\s+<([^>]+)>', line)
                if match:
                    rel_id, evidence_id = match.groups()
                    relation_evidence[rel_id] = evidence_id
            
            # Relationship RO_0002556 pathogen
            elif 'RO_0002556' in line:
                match = re.search(r'<([^>]+)>\s+<[^>]*RO_0002556>\s+<([^>]+)>', line)
                if match:
                    rel_id, pathogen_id = match.groups()
                    relation_pathogen[rel_id] = pathogen_id
            
            # Extract labels
            elif 'rdf-schema#label' in line:
                match = re.search(r'<([^>]+)>\s+<[^>]*rdf-schema#label>\s+"([^"]+)"', line)
                if match:
                    entity_id, label = match.groups()
                    labels[entity_id] = label
        
        # Build result tuples for all disease-pathogen combinations
        for disease_id, rel_ids in disease_relations.items():
            disease_label = labels.get(disease_id, 'Unknown')
            disease_id_str = disease_id.split('/')[-1].replace('_', ':')
            
            for rel_id in rel_ids:
                if rel_id in relation_evidence and rel_id in relation_pathogen:
                    pathogen_id = relation_pathogen[rel_id]
                    evidence_id = relation_evidence[rel_id]
                    
                    pathogen_label = labels.get(pathogen_id, 'Unknown')
                    evidence_label = labels.get(evidence_id, 'Unknown')
                    
                    # Extract pathogen taxid
                    pathogen_taxid = pathogen_id.split('/')[-1].split('_')[-1]
                    
                    pair = DiseasePathogen(
                        disease_id=disease_id_str,
                        disease=disease_label,
                        pathogen_taxid=pathogen_taxid,
                        pathogen=pathogen_label,
                        evidence=evidence_label,
                    )
                    
                    result.add(pair)
        
        result = list(result)
        
        # Cache the result
        with open(cache_path, 'wb') as fp:
            _log(f'Saving to cache: `{cache_path}`.')
            pickle.dump(obj=result, file=fp)
        
        return result
        
    except Exception as e:
        _logger._console(f'Failed to parse PathoPhenoDB RDF data: {e}')
        _log('Returning empty result!')
        return []
