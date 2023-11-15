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

import re
import itertools
import collections
import bs4
import warnings

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.progress as progress
import pypath.share.common as common
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera
import pypath.core.entity as entity


KeggPathway = collections.namedtuple(
    'KeggPathway',
    ['pathway'],
)


def kegg_interactions():
    """
    Downloads and processes KEGG Pathways.
    Returns list of interactions.
    """

    positive_terms = {'activation', 'expression'}
    negative_terms = {'inhibition', 'repression'}
    transc_terms = {'expression', 'repression'}
    mechanism_terms = {
        'phosphorylation',
        'binding/association',
        'dissociation',
        'ubiquitination',
        'dephosphorylation',
        'glycosylation',
        'state change',
        'methylation',
    }
    direct_terms = {'indirect effect'}

    KeggInteraction = collections.namedtuple(
        'KeggInteraction',
        [
            'id_a',
            'id_b',
            'effect',
            'pathway',
            'mechanism',
            'is_direct',
            'transcriptional',
        ],
    )

    rehsa = re.compile(r'.*(hsa[0-9]+).*')
    req_hdrs = [
        'Referer: http://www.genome.jp/kegg-bin/show_pathway'
        '?map=hsa04710&show_description=show'
    ]
    hsa_list = []
    interactions = []

    c = curl.Curl(urls.urls['kegg_pws']['list_url'], silent = True)
    htmllst = c.result
    lstsoup = bs4.BeautifulSoup(htmllst, 'html.parser')

    for a in lstsoup.find_all('a', href = True):
        m = rehsa.match(a['href'])

        if m:
            hsa_list.append((m.groups(0)[0], a.text))

    prg = progress.Progress(
        len(hsa_list), 'Processing KEGG Pathways', 1, percent = False
    )

    for hsa, pw in hsa_list:

        prg.step()
        c = curl.Curl(
            urls.urls['kegg_pws']['kgml_url_2'] % hsa,
            silent = True,
            req_headers = req_hdrs
        )
        kgml = c.result

        with warnings.catch_warnings():

            warnings.simplefilter('ignore')
            kgmlsoup = bs4.BeautifulSoup(kgml, 'html.parser')

        entries = {}

        for ent in kgmlsoup.find_all('entry'):
            gr = ent.find('graphics')

            if gr and 'name' in gr.attrs:
                entries[ent.attrs['id']] = [
                    n.strip()
                    for n in gr.attrs['name'].replace('...', '').split(',')
                ]

        uentries = dict([(eid, common.unique_list(
            common.flat_list([
                mapping.map_name(
                    gn, 'genesymbol', 'uniprot', strict = True) for gn in gns
            ]))) for eid, gns in iteritems(entries)])

        for rel in kgmlsoup.find_all('relation'):

            subtypes = {st.attrs['name'] for st in rel.find_all('subtype')}

            if (
                rel.attrs['entry1'] in uentries and
                rel.attrs['entry2'] in uentries and
                subtypes
            ):

                is_direct = 'indirect effect' not in subtypes
                effect = (
                    'inhibition'
                        if negative_terms & subtypes else
                    'activation'
                        if positive_terms & subtypes else
                    'unknown'
                )
                mechanism = ';'.join(mechanism_terms & subtypes)
                transcriptional = bool(transc_terms & subtypes)

                for u1 in uentries[rel.attrs['entry1']]:

                    for u2 in uentries[rel.attrs['entry2']]:

                        interactions.append(
                            KeggInteraction(
                                id_a = u1,
                                id_b = u2,
                                effect = effect,
                                pathway = pw,
                                mechanism = mechanism,
                                is_direct = is_direct,
                                transcriptional = transcriptional,
                            )
                        )

    prg.terminate()

    return common.unique_list(interactions)


def kegg_pathways():

    data = kegg_interactions()
    pws = common.unique_list(map(lambda i: i[3], data))
    proteins_pws = dict(map(lambda pw: (pw, set([])), pws))
    interactions_pws = dict(map(lambda pw: (pw, set([])), pws))

    for rec in data:

        u1, u2, eff, pw = rec[:4]
        proteins_pws[pw].add(u1)
        proteins_pws[pw].add(u2)
        interactions_pws[pw].add((u1, u2))

    return proteins_pws, interactions_pws


def kegg_pathway_annotations():

    result = collections.defaultdict(set)

    proteins, interactions = kegg_pathways()

    for pathway, uniprots in iteritems(proteins):
        record = KeggPathway(pathway = pathway)

        for uniprot in uniprots:
            result[uniprot].add(record)

    return dict(result)


def kegg_pathway_annotations_pathwaycommons():

    result = collections.defaultdict(set)

    url = urls.urls['kegg_pws']['pw_annot']
    c = curl.Curl(url, large = True, silent = False)

    for row in c.result:

        row = row.strip().split('\t')
        name = row[1].split(';', maxsplit = 1)[0]
        name = name.split(':', maxsplit = 1)[1].strip()
        uniprots = row[2:]

        annot = KeggPathway(pathway = name)

        for uniprot in uniprots:

            result[uniprot].add(annot)

    return dict(result)


def kegg_medicus(max_entity_variations = 10):
    """
    Retrieves and preprocesses the KEGG MEDICUS database. Returns a set of
    raw interaction records (with the original identifiers and some further
    attributes). Nested complexes and protein families are flattened which
    means each interacting pair is either a single protein or a protein
    complex. Then the combination of all variants of each interacting partner
    yields a separate record. E.g. if a family of 3 proteins interacts with
    a protein complex where one of the members can be 2 alternative proteins
    then this interaction yields 6 records.

    max_entity_variations : int
        In KEGG MEDICUS many molecular entities are protein families or
        families of often large and nested protein complexes. By this option
        you can limit largest number of variants a single entity might yield,
        so you won't end up with one complex yielding hundreds of
        combinatiorial variants.
    """


    reentity = re.compile(r'[,\+\(\)]|\w+')
    renminus2 = re.compile(r'\(n(?:-2)?\)')
    renetref = re.compile(r'\[(N|nt)\d{5}\]')


    KeggMedicusRawInteraction = collections.namedtuple(
        'KeggMedicusRawInteraction',
        [
            'id_a',
            'id_b',
            'name_a',
            'name_b',
            'effect',
            'itype',
            'pw_type',
            'type_a',
            'type_b',
            'network_id',
        ],
    )


    i_code = {
        '->': ('post_translational', 'stimulation'),
        '=>': ('transcriptional', 'stimulation'),
        '//': ('post_translational', 'missing'),
        '-|': ('post_translational', 'inhibition'),
        '=|': ('transcriptional', 'inhibition'),
        '--': ('post_translational', 'undirected'),
        '>>': ('post_translational', 'enzyme_enzyme'),
        '==': ('post_translational', 'missing'),
    }


    def process_entity(e):

        if isinstance(e, str):

            e = renminus2.sub('', e)
            e = reentity.findall(e)

        sub = 0
        stack = []
        cplex = False

        for it in e:

            if it == ',':

                continue

            elif it == ')':

                sub -= 1

                if not sub:

                    stack.append(process_entity(this_stack))

                else:

                    this_stack.append(it)

            elif sub:

                this_stack.append(it)

                if it == '(':

                    sub += 1

            elif it == '(':

                if not sub:

                    this_stack = []

                sub += 1

            elif it == '+':

                cplex = True

            else:

                stack.append(it)

        if cplex:

            stack = tuple(stack)

        return stack


    def flatten_entity(e):

        flat = []

        if isinstance(e, str):

            flat.append(e)

        elif isinstance(e, tuple):

            flat.extend(
                itertools.product(*(
                    (c,)
                        if isinstance(c, str) else
                    (flatten_entity(c),)
                        if isinstance(c, tuple) else
                    c
                    for c in e
                ))
            )

        elif isinstance(e, list):

            flat.extend(itertools.chain(*(flatten_entity(c) for c in e)))

        if any(
            any(isinstance(c, list) for c in flate)
            for flate in flat
        ):

            flat = list(
                itertools.chain(*(
                    flatten_entity(flate) for flate in flat
                ))
            )

        flat = [flatten_nested_complex(flate) for flate in flat]

        return flat


    def flatten_nested_complex(cplex):

        if is_nested_complex(cplex):

            cplex = tuple(
                member
                for members in cplex
                for member in (
                    members
                        if isinstance(members, tuple) else
                    (members,)
                )
            )

            if is_nested_complex(cplex):

                cplex = flatten_nested_complex(cplex)

        return cplex


    def is_nested_complex(cplex):

        return (
            isinstance(cplex, tuple) and
            any(isinstance(member, tuple) for member in cplex)
        )


    def get_interactions(connections, enames, pw_type, network_id):

        entities = dict(
            (
                i,
                flatten_entity(process_entity(connections[i]))
            )
            for i in range(0, len(connections), 2)
        )

        for i in range(0, len(connections) - 1, 2):

            itype, effect = i_code[connections[i + 1]]

            if (
                len(entities[i]) > max_entity_variations or
                len(entities[i + 2]) > max_entity_variations
            ):

                continue

            for id_a, id_b in itertools.product(entities[i], entities[i + 2]):

                name_a, type_a = get_name_type(id_a, enames)
                name_b, type_b = get_name_type(id_b, enames)

                yield KeggMedicusRawInteraction(
                    id_a = id_a,
                    id_b = id_b,
                    name_a = name_a,
                    name_b = name_b,
                    effect = effect,
                    itype = itype,
                    pw_type = pw_type,
                    type_a = type_a,
                    type_b = type_b,
                    network_id = network_id,
                )


    def get_name_type(_id, enames):

        return (
            tuple(zip(*(_get_name_type(i, enames) for i in _id)))
                if isinstance(_id, tuple) else
            _get_name_type(_id, enames)
        )


    def _get_name_type(_id, enames):

        if _id not in enames:

            dbget = kegg_dbget(_id)

            if not dbget:

                name, entity_type = (None, None)

            else:

                name = (
                    dbget['Name'][-1]
                        if isinstance(dbget['Name'], list) else
                    dbget['Name']
                )
                entity_type = dbget['Type'].lower()

            enames[_id] = (name, entity_type)

        return enames[_id]


    recollect = re.compile(r'^(GENE|PERTURBANT|VARIANT|METABOLITE)')
    recon = re.compile(r'(->|--|//|-\||=>|>>|=\||==)')
    rewrongspace = re.compile(r'(\d+) (?=\d+)')
    result = set()
    url = urls.urls['kegg_pws']['medicus']
    c = curl.Curl(url, silent = False, large = True)
    enames = {}
    collecting = None

    for row in c.result:

        begin_coll = recollect.match(row)

        if begin_coll:

            collecting = begin_coll.group()
            row = row.split(maxsplit = 1)[-1]

        if collecting:

            if not begin_coll and row[0] != ' ':

                collecting = None
                continue

            if collecting == 'GENE':

                row = row.split(';')[0]

            id_name = row.split(maxsplit = 1)

            if len(id_name) == 2:

                _id, name = id_name

            else:

                _id = id_name[0]
                dbget = kegg_dbget(_id)

                name = (
                    dbget['Name']
                        if 'Name' in dbget else
                    dbget['Composition']
                )

                if isinstance(name, list):

                    name = name[-1]

            enames[_id] = (name.strip(), collecting.lower())

    c.fileobj.seek(0)

    for row in c.fileobj:

        if row.startswith('ENTRY'):

            pw_type = None
            collecting = None
            network_id = row.split()[1]

        elif row.startswith('TYPE'):

            pw_type = row.strip().split()[-1].lower()

        elif row.startswith('  EXPANDED'):

            connections = renetref.sub('', row)
            connections = recon.sub(' \g<1> ', connections)
            connections = rewrongspace.sub('\g<1>,', connections)
            connections = connections.split()[1:]

        elif row.startswith('///'):

            result.update(
                set(get_interactions(
                    connections,
                    enames,
                    pw_type,
                    network_id
                ))
            )

    return result


def kegg_medicus_interactions(max_entity_variations = 10, complexes = False):
    """
    Retrieves and preprocesses human protein-protein and transcriptional
    regulatory interactions from the KEGG MEDICUS database. Optionally
    it returns protein complexes instead of interactions.

    max_entity_variations : int
        In KEGG MEDICUS many molecular entities are protein families or
        families of often large and nested protein complexes. By this option
        you can limit largest number of variants a single entity might yield,
        so you won't end up with one complex yielding hundreds of
        combinatiorial variants.
    complexes : bool
        Return a set of protein complexes instead of a list of molecular
        interactions.
    """

    KeggMedicusInteraction = collections.namedtuple(
        'KeggMedicusInteraction',
        [
            'id_a',
            'id_b',
            'entity_type_a',
            'entity_type_b',
            'interaction_type',
            'effect',
        ]
    )


    result = []
    cplexes = {}


    def process_complex(ids, symbols, types):

        if ids not in cplexes:

            if not all(t == 'gene' for t in types):

                cplexes[ids] = set()

            uniprots = [
                process_protein(id_, symbol)
                for id_, symbol in zip(ids, symbols)
            ]

            this_cplexes = {
                intera.Complex(
                    components = components,
                    sources = 'KEGG-MEDICUS',
                )
                for components in itertools.product(*uniprots)
            }

            cplexes[ids] = this_cplexes

        return cplexes[ids]


    def process_protein(id_, symbol):

        return (
            mapping.map_name(id_, 'entrez', 'uniprot') or
            mapping.map_name(id_, 'genesymbol', 'uniprot')
        )


    def process_partner(ids, symbols, types = None):

        return (
            process_protein(ids, symbols)
                if isinstance(ids, str) else
            process_complex(ids, symbols, types)
        )


    for rec in kegg_medicus(max_entity_variations = max_entity_variations):

        for id_a, id_b in itertools.product(
            process_partner(rec.id_a, rec.name_a, rec.type_a),
            process_partner(rec.id_b, rec.name_b, rec.type_b),
        ):

            if not complexes:

                result.append(
                    KeggMedicusInteraction(
                        id_a = id_a,
                        id_b = id_b,
                        entity_type_a = entity.Entity._get_entity_type(id_a),
                        entity_type_b = entity.Entity._get_entity_type(id_b),
                        interaction_type = rec.itype,
                        effect = rec.effect,
                    )
                )

    return set.union(*cplexes.values()) if complexes else result



def kegg_medicus_complexes(max_entity_variations = 10):
    """
    Extracts a `dict` of protein complexes from the KEGG MEDICUS database.

    max_entity_variations : int
        In KEGG MEDICUS many molecular entities are protein families or
        families of often large and nested protein complexes. By this option
        you can limit largest number of variants a single entity might yield,
        so you won't end up with one complex yielding hundreds of
        combinatiorial variants.
    """

    cplexes = kegg_medicus_interactions(
        max_entity_variations = max_entity_variations,
        complexes = True,
    )

    cplexes = dict((cplex.__str__(), cplex) for cplex in cplexes)

    return cplexes




def kegg_dbget(entry):
    """
    Retrieves an entry (e.g. compounds, network modules) by the KEGG DBGET
    interface (kegg.jp/dbget-bin/www_bget).
    """

    rexa = re.compile(r'\xa0+')
    stripchars = '\r\n; '
    reffields = {'Authors', 'Title', 'Journal'}

    result = {}

    if isinstance(entry, int):

        entry = 'hsa:%u' % entry

    if entry.isdigit():

        entry = 'hsa:%s' % entry

    url = urls.urls['kegg_pws']['dbget'] % entry
    c = curl.Curl(url, silent = True, large = False)
    soup = bs4.BeautifulSoup(c.result, 'html.parser')

    tbl = soup.find_all('table', limit = 4)

    if not tbl:

        return None

    tbl = tbl[-1]

    collecting_ref = False
    last_ref = {}

    for row in tbl.findChildren('tr', recursive = False):

        key = row.find('th').text.strip()
        td = row.find('td')

        if collecting_ref:

            if key in reffields:

                last_ref[key] = td.text
                continue

            else:

                if 'References' not in result:

                    result['References'] = []

                result['References'].append(last_ref)
                last_ref = {}
                collecting_ref = False

        if key == 'Reference':

            collecting_ref = True
            last_ref['PMID'] = re.findall(r'\d+', td.text)[-1]
            continue

        subtbl = td.find_all('table')

        if subtbl:

            value = {}

            for st in subtbl:

                for subrow in st.find_all('tr'):

                    subtd = subrow.find_all('td')

                    if len(subtd) > 1 and subtd[1].text:

                        value[rexa.sub('', subtd[0].text)] = (
                            subtd[1].text.strip(stripchars)
                        )

                    else:

                        subcontent = rexa.sub(' ', subtd[0].text).split()

                        if len(subcontent) > 1:

                            value[subcontent[0]] = (
                                subcontent[1].strip(stripchars)
                            )

        else:

            value = rexa.sub(' ', td.text).strip(stripchars)

            if '\n' in value:

                value = [
                    lval.strip(stripchars)
                    for lval in re.split(r'\s*[\n\r]+\s*', value)
                ]


        if key == 'Entry':

            value, result['Type'] = next(value.items().__iter__())

        result[key] = value

    return result
