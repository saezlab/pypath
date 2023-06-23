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

from typing import Union

import re
import collections

import pypath.utils.mapping as mapping
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.formats.obo as obo


def hpo_annotations() -> dict[str, set[tuple]]:
    """
    Human Phenotype Ontology annotations.

    Returns
        Dict of proteins as keys and sets of HPO annotations as values.
    """

    url = urls.urls['hpo']['gene']
    c = curl.Curl(url, large = True, silent = False)
    _ = next(c.result)

    fields = ('entrez_gene_id','entrez_gene_symbol','hpo_id')

    HpoAnnotation = collections.namedtuple(
        'HpoAnnotation',
        fields,defaults = ("",) * len(fields)
    )

    result = collections.defaultdict(set)

    for r in c.result:

        r = r.strip().split('\t')

        uniprots = mapping.map_name(r[0], 'entrez', 'uniprot')

        for uniprot in uniprots:

            result[uniprot].add(
                HpoAnnotation(
                    entrez_gene_id = r[0],
                    entrez_gene_symbol = r[1],
                    hpo_id = r[2],
                )
            )

    return dict(result)


def hpo_terms() -> dict[str, str]:
    """
    Human Phenotype Ontology accession to term mapping.
    """

    return hpo_ontology()['terms']


def hpo_diseases() -> dict[str, set[tuple]]:
    """
    HPO term-disease relationships from Human Phenotype Ontology.

    Returns
        A set of disease records for each HPO term.
    """

    url = urls.urls['hpo']['disease']
    c = curl.Curl(url, large = True, silent = False)

    HpoDisease = collections.namedtuple(
        'HpoDisease',
        (
            'omim',
            'name',
            'pmid',
            'qualifier',
            'evidence',
            'onset',
            'frequency',
            'sex',
            'modifier',
            'aspect',
        ),
    )

    result = collections.defaultdict(set)

    for r in c.result:

        if r[0] == '#': continue

        r = r.split('\t')

        pmid = re.sub('^PMID:', '', r[4]) if r[4][:4] == 'PMID' else None

        result[r[3]].add(
            HpoDisease(
                omim = r[0],
                name = r[1],
                pmid = pmid,
                qualifier = r[2] or None,
                evidence = r[5] or None,
                onset = r[6] or None,
                frequency = r[7] or None,
                sex = r[8] or None,
                modifier = r[9] or None,
                aspect = r[10],
            )
        )

    return dict(result)


def hpo_ontology() -> dict[str, dict[str, Union[str, set[str]]]]:
    """
    Ontology data from HPO.

    Returns
        Five dictionaries with term names, term definitions, parents in the
        ontology tree, term synonyms and cross references to other databases.
        The dicts "terms" and "defs" are one-to-one, while "parents",
        "synonyms" and "xrefs" are one-to-many mappings, the keys are always
        HPO terms.
    """

    url = urls.urls['hpo']['ontology']
    reader = obo.Obo(url)

    result = {
        'terms': {},
        'defs': {},
        'parents': collections.defaultdict(set),
        'synonyms': collections.defaultdict(set),
        'xrefs': collections.defaultdict(set),
    }

    for r in reader:

        if r.stanza != 'Term': continue

        if (
            r.name is None or
            r.name.value == 'obsolete' or
            r.attrs.get('is_obsolete')
        ):
            continue

        term = r.id.value

        name = (
            (r.name.value, r.name.modifiers)
                if r.name.modifiers else
            r.name.value
        )

        if isinstance(name, tuple): name = ' '.join(n for n in name if n)

        result['terms'][term] = name

        result['defs'][term] = r.definition.value if r.definition else None

        for key, obokey in (
            ('parents', 'is_a'),
            ('synonyms', 'synonym'),
            ('xrefs', 'xref'),
        ):

            proc = (
                lambda x: tuple(x.split(':'))
                    if key == 'xrefs' else
                lambda x: x
            )

            for x in r.attrs.get(obokey, ()):
                y = proc(x.value)
                result[key][term].update(
                    {
                        y(x.value)
                        if type(y) != tuple else
                        y
                    }
                )

    return {k: dict(v) for k, v in result.items()}
