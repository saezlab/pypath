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

from past.builtins import xrange, range

import collections

from lxml import etree

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.uniprot_db as uniprot_db
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping


def locate_localizations(
        organism = 9606,
        literature = True,
        external = True,
        predictions = False,
    ):

    record = collections.namedtuple(
        'LocateAnnotation',
        ('source', 'location', 'cls', 'pmid', 'score'),
    )
    record.__new__.__defaults__ = (None, None, None)

    organism_uniprots = set(
        uniprot_db.all_uniprots(organism = organism, swissprot = True)
    )

    organism_str = taxonomy.taxids[organism]
    url = urls.urls['locate']['url_rescued'] % organism_str
    fname = url.split('/')[-1][:-4]

    c = curl.Curl(
        url,
        large = True,
        default_mode = 'rb',
        silent = False,
        files_needed = [fname],
    )
    c.result[fname]

    parser = etree.iterparse(c.result[fname], events = ('start', 'end'))

    result = collections.defaultdict(set)
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == 'LOCATE_protein':

            tag_protein = elem.find('protein')
            this_uniprot = None
            this_uniprots = None
            this_entrez  = None
            this_organism = (
                tag_protein.find('organism').text
                    if tag_protein is not None else
                None
            )
            this_class = (
                tag_protein.find('class').text
                    if tag_protein is not None else
                None
            )

            xrefs = elem.find('xrefs')

            if xrefs is None:
                continue

            for xref in xrefs.findall('xref'):
                src = xref.find('source')
                src_name = src.find('source_name').text

                if src_name == 'UniProtKB-SwissProt':
                    this_uniprot = src.find('accn').text

                if src_name == 'Entrez Gene':
                    this_entrez = src.find('accn').text

                if src_name == 'UniProt/SPTrEMBL' and this_uniprot is None:
                    this_uniprot = src.find('accn').text

            # if we don't know what it is, does not make sense to proceed
            if this_uniprot is None and this_entrez is None:
                continue

            if this_uniprot:
                this_uniprots = mapping.map_name(
                    this_uniprot,
                    'uniprot',
                    'uniprot',
                    ncbi_tax_id = organism,
                )

            if not this_uniprots and this_entrez:
                this_uniprots = mapping.map_name(
                    this_entrez,
                    'entrez',
                    'uniprot',
                    ncbi_tax_id = organism,
                )

            this_uniprots = set(this_uniprots) & organism_uniprots

            # if we don't know what it is, does not make sense to proceed
            if not this_uniprots:
                continue

            if external:
                # External database annotations
                extannot = elem.find('externalannot')

                if extannot is not None:
                    for extannotref in extannot.findall('reference'):
                        sources = []

                        for src in extannotref.findall('source'):
                            src_name = src.find('source_name')

                            if src_name is not None:
                                sources.append(src_name.text)

                        sources = ';'.join(sources) if sources else None
                        locations =  extannotref.find('locations')

                        if locations is not None:
                            for location in locations.findall('location'):
                                for loc in location.iterchildren():
                                    if loc.tag[:4] == 'tier':
                                        this_loc = loc.text.lower().split(',')

                                        for uniprot in this_uniprots:
                                            for _loc in this_loc:
                                                result[uniprot].add(record(
                                                    source = sources,
                                                    location = _loc.strip(),
                                                    cls = this_class,
                                                    score = None,
                                                ))

            if predictions:
                # Predictions
                sclpred = elem.find('scl_prediction')

                if sclpred is not None:
                    for sclpred_src in sclpred.findall('source'):
                        score = float(sclpred_src.find('evaluation').text)

                        if score == 0.0:
                            continue

                        this_src = sclpred_src.find('method').text
                        this_loc = sclpred_src.find('location').text.lower()

                        if this_loc == 'no prediction':
                            continue

                        for uniprot in this_uniprots:
                            result[uniprot].add(record(
                                source = this_src,
                                location = this_loc,
                                cls = this_class,
                                score = score,
                            ))

            if literature:
                # Literature curation
                lit = elem.find('literature')

                if lit is not None:

                    for litref in lit.findall('reference'):

                        locs = set()

                        for lloc in (
                            litref.find('locations').findall('location')
                        ):

                            for loc in lloc.iterchildren():
                                if loc.tag[:4] == 'tier':
                                    locs.add(loc.text.lower())

                        pmid = litref.find('source')
                        pmid = (
                            None if pmid is None else pmid.find('accn').text
                        )

                        for loc in locs:

                            for uniprot in this_uniprots:

                                result[uniprot].add(
                                    record(
                                        source = 'literature',
                                        location = loc,
                                        pmid = pmid,
                                        cls = this_class,
                                        score = None,
                                    )
                                )

        used_elements.append(elem)

        # removing used elements to keep memory low
        if len(used_elements) > 1000:
            for _ in xrange(500):
                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return dict(result)
