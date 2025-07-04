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

from future.utils import iteritems

import csv
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common
import pypath.inputs.science as science


def get_proteinatlas(normal = True, pathology = True, cancer = True):

    result = {
        'normal':    collections.defaultdict(dict),
        'pathology': collections.defaultdict(dict),
    }


    def line(l):

        return l.strip('\n\r').split('\t')


    if normal:

        c = curl.Curl(urls.urls['proteinatlas']['normal'],
                    silent = False, large = True)
        fp = list(c.result.values())[0]
        hdr = line(fp.readline().strip())

        for l in fp:
            l = line(l)
            
            if len(l) < 7:  # Skip incomplete lines
                continue

            uniprots = mapping.map_name(l[1], 'genesymbol', 'uniprot')
            # New format: Gene, Gene name, Tissue, IHC tissue name, Cell type, Level, Reliability
            tissue = '%s:%s:%s' % (l[2], l[3], l[4])  # Tissue:IHC tissue name:Cell type

            for u in uniprots:
                result['normal'][tissue][u] = (l[5], l[6].strip())  # Level, Reliability


    if cancer or pathology:

        c = curl.Curl(urls.urls['proteinatlas']['pathology'],
                    silent = False, large = True)
        fp = list(c.result.values())[0]
        hdr = line(fp.readline())

        for l in fp:

            l = line(l)
            
            if len(l) < 7:  # Skip incomplete lines
                continue
                
            uniprots = mapping.map_name(l[1], 'genesymbol', 'uniprot')
            tissue   = l[2]  # Cancer type

            # New format: Gene, Gene name, Cancer, High, Medium, Low, Not detected
            values = {
                'High': int(l[3]) if l[3].isdigit() else 0,
                'Medium': int(l[4]) if l[4].isdigit() else 0,
                'Low': int(l[5]) if l[5].isdigit() else 0,
                'Not detected': int(l[6]) if l[6].isdigit() else 0,
            }

            for u in uniprots:

                result['pathology'][tissue][u] = values

    return dict((k, dict(v)) for k, v in iteritems(result))


def proteinatlas_annotations(normal = True, pathology = True, cancer = True):

    LEVELS = ('Not detected', 'Low', 'Medium', 'High')


    ProteinatlasAnnotation = collections.namedtuple(
        'ProtainatlasAnnotation',
        [
            'organ',
            'tissue',
            'cell_type',
            'level',
            'status',
            'n_not_detected',
            'n_low',
            'n_medium',
            'n_high',
            'prognostic',
            'favourable',
            'score',
            'pathology',
        ],
    )
    ProteinatlasAnnotation.__new__.__defaults__ = (
        (None,) * 5 + (False, False, None, False)
    )


    def n_or_none(ex, key):

        return ex[key] if key in ex else None


    data = get_proteinatlas(
        normal = normal,
        pathology = pathology,
        cancer = cancer,
    )

    result = collections.defaultdict(set)

    if normal:

        for tissue, gex in iteritems(data['normal']):

            organ = tissue
            cell_type = ''

            if ':' in tissue:
                parts = tissue.split(':')
                if len(parts) >= 3:
                    organ, tissue, cell_type = parts[0], parts[1], parts[2]
                elif len(parts) == 2:
                    organ, tissue = parts[0], parts[1]
                else:
                    organ = parts[0]

            organ = organ.strip()
            tissue = tissue.strip()
            cell_type = cell_type.strip()

            for uniprot, ex in iteritems(gex):
                uniprots = mapping.map_name(uniprot, 'uniprot', 'uniprot')

                for _uniprot in uniprots:
                    result[_uniprot].add(
                        ProteinatlasAnnotation(
                            organ = organ,
                            tissue = tissue,
                            cell_type = cell_type,
                            level = ex[0],
                            status = ex[1],
                        )
                    )

    if pathology or cancer:
        for condition, gex in iteritems(data['pathology']):
            for uniprot, ex in iteritems(gex):
                try:
                    effect, score = next(
                        i for i in iteritems(ex) if i[0] not in LEVELS
                    )
                    prognostic = not effect.startswith('unprognostic')
                    favourable = not effect.endswith('unfavourable')

                except StopIteration:
                    prognostic, favourable, score = None, None, None

                uniprots = mapping.map_name(uniprot, 'uniprot', 'uniprot')

                for _uniprot in uniprots:
                    result[_uniprot].add(
                        ProteinatlasAnnotation(
                            organ = condition,
                            tissue = condition,
                            cell_type = None,
                            level = max(
                                (i for i in iteritems(ex) if i[0] in LEVELS),
                                key = lambda i: i[1],
                                default = (None,),
                            )[0],
                            status = None,
                            n_not_detected = n_or_none(ex, 'Not detected'),
                            n_low = n_or_none(ex, 'Low'),
                            n_medium = n_or_none(ex, 'Medium'),
                            n_high = n_or_none(ex, 'High'),
                            prognostic = prognostic,
                            favourable = favourable,
                            score = score,
                            pathology = True,
                        )
                    )

    return dict(result)


def proteinatlas_interactions():
    """
    Downloads and processes interaction consensus data from Human Protein Atlas.
    """

    ProteinatlasInteraction = collections.namedtuple(
        'ProteinatlasInteraction',
        [
            'uniprot_a',
            'uniprot_b',
            'gene_a',
            'gene_b',
            'interaction_type',
            'evidence',
        ],
    )

    url = urls.urls['proteinatlas']['interactions']

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'r',
    )
    
    # Find the interaction consensus file
    file_key = None
    for key in c.result.keys():
        if 'interaction_consensus' in key:
            file_key = key
            break
    
    if not file_key:
        raise ValueError('Could not find interaction_consensus.tsv in downloaded files')
    
    reader = csv.DictReader(
        c.result[file_key],
        delimiter = '\t',
    )

    result = []

    for rec in reader:
        result.append(ProteinatlasInteraction(
            uniprot_a = rec.get('Uniprot_A', ''),
            uniprot_b = rec.get('Uniprot_B', ''), 
            gene_a = rec.get('Gene_A', ''),
            gene_b = rec.get('Gene_B', ''),
            interaction_type = rec.get('Interaction_type', ''),
            evidence = rec.get('Evidence', ''),
        ))

    return result


def proteinatlas_subcellular_annotations():

    ProteinatlasSubcellularAnnotation = collections.namedtuple(
        'ProteinatlasSubcellularAnnotation',
        [
            'location',
            'status',
            'reliability',
            'main_location',
            'additional_location',
        ],
    )

    url = urls.urls['proteinatlas']['subcell']

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'r',
    )
    
    # Find the correct TSV file in the zip
    file_key = None
    for key in c.result.keys():
        if 'subcellular_location.tsv' in key:
            file_key = key
            break
    
    if not file_key:
        raise ValueError('Could not find subcellular_location.tsv in downloaded files')
    
    reader = csv.DictReader(
        c.result[file_key],
        delimiter = '\t',
    )

    result = collections.defaultdict(set)

    for rec in reader:
        uniprots = mapping.map_name(rec['Gene name'], 'genesymbol', 'uniprot')

        for uniprot in uniprots:
            # Handle the new format with main and additional locations
            main_locations = rec.get('Main location', '').split(';') if rec.get('Main location') else []
            additional_locations = rec.get('Additional location', '').split(';') if rec.get('Additional location') else []
            
            all_locations = main_locations + additional_locations
            
            for location in all_locations:
                if location.strip():
                    result[uniprot].add(ProteinatlasSubcellularAnnotation(
                        location = location.strip(),
                        status = 'main' if location in main_locations else 'additional',
                        reliability = rec.get('Reliability', ''),
                        main_location = rec.get('Main location', ''),
                        additional_location = rec.get('Additional location', ''),
                    ))
            
            # Also handle Enhanced, Supported, Approved, Uncertain columns if they exist
            for status in ('Enhanced', 'Supported', 'Approved', 'Uncertain'):
                if rec.get(status):
                    for location in rec[status].split(';'):
                        if location.strip():
                            result[uniprot].add(ProteinatlasSubcellularAnnotation(
                                location = location.strip(),
                                status = status.lower(),
                                reliability = rec.get('Reliability', ''),
                                main_location = rec.get('Main location', ''),
                                additional_location = rec.get('Additional location', ''),
                            ))

    return dict(result)


def proteinatlas_interactions():
    """
    Downloads and processes interaction consensus data from Human Protein Atlas.
    """

    ProteinatlasInteraction = collections.namedtuple(
        'ProteinatlasInteraction',
        [
            'uniprot_a',
            'uniprot_b',
            'gene_a',
            'gene_b',
            'interaction_type',
            'evidence',
        ],
    )

    url = urls.urls['proteinatlas']['interactions']

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'r',
    )
    
    # Find the interaction consensus file
    file_key = None
    for key in c.result.keys():
        if 'interaction_consensus' in key:
            file_key = key
            break
    
    if not file_key:
        raise ValueError('Could not find interaction_consensus.tsv in downloaded files')
    
    reader = csv.DictReader(
        c.result[file_key],
        delimiter = '\t',
    )

    result = []

    for rec in reader:
        result.append(ProteinatlasInteraction(
            uniprot_a = rec.get('Uniprot_A', ''),
            uniprot_b = rec.get('Uniprot_B', ''), 
            gene_a = rec.get('Gene_A', ''),
            gene_b = rec.get('Gene_B', ''),
            interaction_type = rec.get('Interaction_type', ''),
            evidence = rec.get('Evidence', ''),
        ))

    return result


def proteinatlas_secretome_annotations():

    ProteinatlasSecretomeAnnotation = collections.namedtuple(
        'ProteinatlasSecretomeAnnotation',
        [
            'mainclass',
            'secreted',
        ],
    )

    # Use local rescued file if available, otherwise try to download
    url = urls.urls['proteinatlas'].get('secretome_rescued', urls.urls['proteinatlas']['secretome'])
    
    if 'secretome_rescued' in urls.urls['proteinatlas']:
        path = url  # Use local file path directly
    else:
        path = science.science_download(url)
        
    reader = inputs_common.read_xls(path)[1:]
    result = collections.defaultdict(set)

    for rec in reader:
        # Skip if not enough columns
        if len(rec) < 4:
            continue
            
        # Column 2 should contain UniProt IDs
        uniprot_column = rec[2] if rec[2] else ''
        
        for uniprot_original in uniprot_column.split(','):
            uniprot_original = uniprot_original.strip()
            if not uniprot_original:
                continue
                
            uniprots = mapping.map_name(
                uniprot_original,
                'uniprot',
                'uniprot',
            )

            for uniprot in uniprots:
                result[uniprot].add(ProteinatlasSecretomeAnnotation(
                    mainclass = rec[3] if len(rec) > 3 and rec[3] else '',
                    secreted = 'secreted' in rec[3].lower() if len(rec) > 3 and rec[3] else False,
                ))

    return dict(result)


def proteinatlas_interactions():
    """
    Downloads and processes interaction consensus data from Human Protein Atlas.
    """

    ProteinatlasInteraction = collections.namedtuple(
        'ProteinatlasInteraction',
        [
            'uniprot_a',
            'uniprot_b',
            'gene_a',
            'gene_b',
            'interaction_type',
            'evidence',
        ],
    )

    url = urls.urls['proteinatlas']['interactions']

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'r',
    )
    
    # Find the interaction consensus file
    file_key = None
    for key in c.result.keys():
        if 'interaction_consensus' in key:
            file_key = key
            break
    
    if not file_key:
        raise ValueError('Could not find interaction_consensus.tsv in downloaded files')
    
    reader = csv.DictReader(
        c.result[file_key],
        delimiter = '\t',
    )

    result = []

    for rec in reader:
        result.append(ProteinatlasInteraction(
            uniprot_a = rec.get('Uniprot_A', ''),
            uniprot_b = rec.get('Uniprot_B', ''), 
            gene_a = rec.get('Gene_A', ''),
            gene_b = rec.get('Gene_B', ''),
            interaction_type = rec.get('Interaction_type', ''),
            evidence = rec.get('Evidence', ''),
        ))

    return result
