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
from future.utils import iteritems

import os
import pickle
import re
import itertools

import xml.etree.cElementTree as ET

import pypath.share.progress as progress
import pypath.utils.taxonomy as taxonomy
import pypath.internals.intera as intera
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.uniprot_db as uniprot_db
import pypath.inputs.common as inputs_common
import pypath.inputs.homologene as homologene
import pypath.utils.mapping as mapping
import pypath.share.common as common
import pypath.share.session as session

_logger = session.Logger(name = 'phosphosite_input')


def phosphosite_enzyme_substrate(
        raw = True,
        organism = 'human',
        strict = True,
    ):
    """
    Downloads and preprocesses phosphorylation site data from PhosphoSitePlus.
    """

    url = urls.urls['psite_kin']['url']
    c = curl.Curl(
        url,
        silent = False,
        compr = 'gz',
        encoding = 'iso-8859-1',
        large = True,
    )
    orto = {}
    data = c.result
    cols = {
        'kinase': 2,
        'kinase_org': 3,
        'substrate': 6,
        'substrate_org': 8,
        'residue': 9,
        'motif': 11
    }
    data = inputs_common.read_table(
        cols = cols,
        fileObject = data,
        sep = '\t',
        hdr = 4,
    )
    result = []
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(_*)([A-Za-z]+)(_*)')

    for r in data:

        if organism is None or \
            ((r['kinase_org'] == organism or not strict) and \
            r['substrate_org'] == organism):

            if r['kinase_org'] != organism:
                korg = r['kinase_org']
                # attempting to map by orthology:
                if korg in taxonomy.taxa and organism in taxonomy.taxa:

                    ktaxid = taxonomy.taxa[korg]
                    taxid = taxonomy.taxa[organism]

                    if korg not in orto:

                        orto[korg] = homologene.homologene_dict(
                            ktaxid,
                            taxid,
                            'refseqp',
                        )

                    korg_refseq = mapping.map_name(r['kinase'],
                                                    'uniprot',
                                                    'refseqp',
                                                    ktaxid)

                    kin_uniprot = \
                        list(
                            itertools.chain(
                                *map(
                                    lambda ors:
                                        mapping.map_name(ors,
                                                        'refseqp',
                                                        'uniprot',
                                                        taxid),
                                    itertools.chain(
                                        *map(
                                            lambda rs:
                                                orto[korg][rs],
                                            filter(
                                                lambda rs:
                                                    rs in orto[korg],
                                                korg_refseq
                                            )
                                        )
                                    )
                                )
                            )
                        )
            else:
                kin_uniprot = [r['kinase']]

            for kinase in kin_uniprot:

                r['resaa'] = r['residue'][0]
                r['resnum'] = int(non_digit.sub('', r['residue'][1:]))
                mot = motre.match(r['motif'])

                # excluding e.g. Q12809_VAR_014388
                r['substrate'] = r['substrate'].split('_')[0]
                sisoform = 1 if '-' not in r['substrate'] else \
                    int(r['substrate'].split('-')[1])
                r['substrate'] = r['substrate'].split('-')[0]

                kisoform = (
                    1 if '-' not in kinase else int(kinase.split('-')[1])
                )
                kinase = kinase.split('-')[0]

                r['substrate'] = r['substrate'].split('-')[0]

                if mot:
                    r['start'] = r['resnum'] - 7 + len(mot.groups()[0])
                    r['end'] = r['resnum'] + 7 - len(mot.groups()[2])
                    r['instance'] = r['motif'].replace('_', '').upper()
                else:
                    r['start'] = None
                    r['end'] = None
                    r['instance'] = None

                if raw:
                    r['kinase'] = kinase
                    result.append(r)
                else:
                    res = intera.Residue(r['resnum'], r['resaa'],
                                         r['substrate'],
                                         isoform = sisoform)

                    mot = intera.Motif(
                        r['substrate'],
                        r['start'],
                        r['end'],
                        instance = r['instance'],
                        isoform = sisoform)

                    ptm = intera.Ptm(protein = r['substrate'],
                                    residue = res,
                                    motif = mot,
                                    typ = 'phosphorylation',
                                    source = 'PhosphoSite',
                                    isoform = sisoform)

                    dom = intera.Domain(protein = kinase, isoform = kisoform)

                    dommot = intera.DomainMotif(
                        domain = dom, ptm = ptm, sources = ['PhosphoSite'])

                    result.append(dommot)

    return result


def phosphosite_ptm_orthology():
    """
    Returns an orthology translation dict of phosphosites
    based on phosphorylation sites table from PhosphoSitePlus.
    In the result all PTMs represented by a tuple of the following
    6 elements: UniProt ID, isoform (int), residue one letter code,
    residue number (int), NCBI Taxonomy ID (int), modification type.

    :param int source: Source taxon (NCBI Taxonomy).
    :param int target: Target taxon (NCBI Taxonomy).
    """

    result = {}

    nondigit = re.compile(r'[^\d]+')

    unknown_taxa = set([])

    for typ in common.psite_mod_types:

        groups = {}

        url = urls.urls['psite_%s' % typ[0]]['url']
        c = curl.Curl(url, silent = False, large = True)

        data = c.result

        for _ in xrange(4):

            __ = next(data)

        for r in data:

            r = common.decode(r, 'utf-8').split('\t')

            if len(r) < 10:

                continue

            uniprot = r[2]
            isoform = 1 if '-' not in uniprot else int(uniprot.split('-')[1])
            uniprot = uniprot.split('-')[0]
            aa = r[4][0]
            num = int(nondigit.sub('', r[4]))
            if r[6] not in taxonomy.taxa:
                unknown_taxa.add(r[6])
                continue

            tax = taxonomy.taxa[r[6]]
            group = int(r[5])

            this_site = (uniprot, isoform, aa, num, tax, typ[1])

            if group not in groups:
                groups[group] = set([])

            groups[group].add(this_site)

        for group, sites in iteritems(groups):

            for site1 in sites:

                for site2 in sites:

                    if site1[4] == site2[4]:

                        continue

                    if site1 not in result:

                        result[site1] = {}

                    if site2[4] not in result[site1]:

                        result[site1][site2[4]] = set([])

                    result[site1][site2[4]].add(site2)

    if len(unknown_taxa):

        _logger._log(
            'Unknown organisms encountered: %s' %
            ', '.join(sorted(unknown_taxa))
        )

    return result


def phosphosite_ptms(organism = 'human'):
    """
    Downloads the phosphorylation site dataset from PhosphoSitePlus.
    """

    result = []
    url = urls.urls['psite_p']['url']
    nondigit = re.compile(r'[^\d]+')
    remot = re.compile(r'(_*)([A-Za-z]+)(_*)')

    c = curl.Curl(url, silent = False, large = True)
    data = c.result

    for _ in xrange(4):

        _ = next(c.result)

    for r in data:

        r = r.split('\t')

        if len(r) > 9 and (organism is None or r[6] == organism):

            uniprot = r[2]
            isoform = 1 if '-' not in uniprot else int(uniprot.split('-')[1])
            uniprot = uniprot.split('-')[0]
            typ = r[3].lower()

            if len(typ) == 0:

                typ = r[4].split('-')[1] if '-' in r[4] else None

            aa = r[4][0]
            num = int(nondigit.sub('', r[4]))
            motif = remot.match(r[9])

            if motif:

                start = num - 7 + len(motif.groups()[0])
                end = num + 7 - len(motif.groups()[2])
                instance = r[9].replace('_', '').upper()

            else:

                start = None
                end = None
                instance = None

            res = intera.Residue(
                num,
                aa,
                uniprot,
                isoform = isoform,
            )
            mot = intera.Motif(
                uniprot,
                start,
                end,
                instance = instance,
                isoform = isoform,
            )
            ptm = intera.Ptm(
                uniprot,
                typ = typ,
                motif = mot,
                residue = res,
                evidences = 'PhosphoSite',
                isoform = isoform,
            )
            result.append(ptm)

    return result


def phosphosite_regsites():
    """
    Downloads and preprocesses the regulatory sites dataset from
    PhosphoSitePlus. This data provides information about which
    proteins a PTM disrupts or induces the interaction with.
    """

    kwds_pos = {
        'enzymatic activity, induced',
        'activity, induced',
        'protein stabilization',
        'receptor inactivation, inhibited',
        'receptor desensitization, inhibited',
        'receptor internalization, inhibited',
        'receptor recycling, induced'
    }

    kwds_neg = {
        'enzymatic activity, inhibited',
        'activity, inhibited',
        'protein degradation',
        'receptor inactivation, induced',
        'receptor desensitization, induced',
        'receptor internalization, induced',
        'receptor recycling, inhibited'
    }

    url = urls.urls['psite_reg']['url']
    c = curl.Curl(url, silent = False, compr = 'gz',
                  encoding = 'iso-8859-1', large = True)
    data = c.result
    cols = {
        'uniprot': 3,
        'organism': 6,
        'mod': 7,
        'on_function': 11,
        'on_process': 12,
        'on_interact': 13,
        'pmids': 15,
        'comments': 19
    }

    data = inputs_common.read_table(
        cols = cols,
        fileObject = data,
        sep = '\t',
        hdr = 4,
    )
    regsites = {}

    for r in data:

        interact = [[y.replace(')', '').strip() for y in x.split('(')]
                    for x in r['on_interact'].strip().split(';') if len(x) > 0]
        induces = [x[0] for x in interact if x[1] == 'INDUCES']
        disrupts = [x[0] for x in interact if x[1] == 'DISRUPTS']
        mod = r['mod']
        modt = r['mod'].split('-')
        mod = list(modt[0])
        aa = mod.pop(0)
        modt = modt[1]
        res = ''.join(mod)
        isoform = (
            int(r['uniprot'].split('-')[1])
            if '-' in r['uniprot']
            else 1
        )
        uniprot = r['uniprot'].split('-')[0]

        if uniprot not in regsites:
            regsites[uniprot] = []

        function = set(map(lambda f: f.strip(),
                           r['on_function'].split(';')))

        regsites[uniprot].append({
            'aa':       aa,
            'res':      res,
            'modt':     modt,
            'organism': r['organism'],
            'pmids':    set(map(lambda f: f.strip(),
                                r['pmids'].split(';'))),
            'induces':  induces,
            'disrupts': disrupts,
            'isoform':  isoform,
            'function': function,
            'process':  set(map(lambda f: f.strip(),
                                r['on_process'].split(';'))),
            'comments': r['comments'],
            'positive': bool(kwds_pos & function),
            'negative': bool(kwds_neg & function)
        })

    return regsites


def phosphosite_regsites_one_organism(organism = 9606):
    """
    Returns PhosphoSitePlus regulatory sites translated to
    one organism by orthology. Residue numbers will be translated
    where necessary, while gene symbols will be translated to
    UniProt IDs of the given organism.
    This works with human, mouse or rat.

    :param int organism:
        NCBI Taxonomy ID of the target organism. In this
        method possible values are human, mouse or rat, as these species
        provide the vast majority of the data, and are close enough to each
        other that the sites can be safely translated between orthologous
        proteins by sequence alignement.
    """

    def genesymbols2uniprots(genesymbols, tax):
        return (
            set(
                itertools.chain(
                    *map(
                        lambda gs:
                            mapping.map_name(
                                gs,
                                'genesymbol',
                                'uniprot',
                                ncbi_tax_id = tax,
                            ),
                        genesymbols
                    )
                )
            )
        )

    def translate_uniprots(uniprots, homo):
        return (
            set(
                itertools.chain(
                    *map(
                        lambda usrc:
                            homo[usrc] if usrc in homo else [],
                        uniprots
                    )
                )
            )
        )

    result = {}

    organisms = set([9606, 10090, 10116])

    mod_types = dict(common.psite_mod_types2)

    regsites = phosphosite_regsites()

    other_organisms = organisms - set([organism])

    homology = (
        dict(
            map(
                lambda other:
                    (
                        other,
                        homologene.homologene_uniprot_dict(
                            source = other,
                            target = organism,
                        )
                    ),
                other_organisms
            )
        )
    )

    ptm_homology = phosphosite_ptm_orthology()

    proteome = uniprot_db.all_uniprots(
        organism = organism,
        swissprot = 'YES',
    )

    for substrate, regs in iteritems(regsites):

        subs = []

        if substrate in proteome:
            subs = [substrate]
        else:
            for other, homo in iteritems(homology):
                if substrate in homo:
                    subs = homo[substrate]

        for sub in subs:

            if sub not in result:
                result[sub] = {}

            for reg in regs:

                reg_organism = taxonomy.taxa[reg['organism']]

                if reg_organism not in organisms:
                    continue

                if reg['modt'] not in mod_types:

                    _logger._log(
                        'Unknown PhosphoSite modification '
                        'type code: %s' % reg['modt']
                    )
                    continue

                mod_type = mod_types[reg['modt']]
                resnum = int(reg['res'])

                psite_key = (
                    substrate,
                    reg['isoform'],
                    reg['aa'],
                    resnum,
                    reg_organism,
                    mod_type,
                )

                if reg_organism != organism:

                    regs_target = []
                    disrupts    = []
                    induces     = []

                    if psite_key in ptm_homology:

                        if organism in ptm_homology[psite_key]:

                            regs_target = ptm_homology[psite_key][organism]

                    if len(regs_target):

                        disrupts = genesymbols2uniprots(
                            reg['disrupts'],
                            reg_organism,
                        )
                        disrupts = translate_uniprots(
                            disrupts,
                            homology[reg_organism],
                        )
                        induces  = genesymbols2uniprots(
                            reg['induces'],
                            reg_organism,
                        )
                        induces  = translate_uniprots(
                            induces,
                            homology[reg_organism],
                        )

                else:

                    regs_target = [psite_key]

                    disrupts = genesymbols2uniprots(reg['disrupts'], organism)
                    induces  = genesymbols2uniprots(reg['induces'], organism)

                for regt in regs_target:

                    modkey = (regt[2], regt[3], regt[5])

                    if modkey not in result[sub]:

                        result[sub][modkey] = {
                            'induces':  set([]),
                            'disrupts': set([]),
                            'pmids':    set([]),
                            'isoforms': set([]),
                            'process':  set([]),
                            'function': set([]),
                            'positive': False,
                            'negative': False,
                            'comments': []
                        }

                    result[sub][modkey]['induces'].update(induces)
                    result[sub][modkey]['disrupts'].update(disrupts)
                    result[sub][modkey]['process'].update(reg['process'])
                    result[sub][modkey]['function'].update(reg['function'])
                    result[sub][modkey]['isoforms'].update([regt[1]])
                    result[sub][modkey]['pmids'].update(reg['pmids'])
                    result[sub][modkey]['positive'] = \
                        result[sub][modkey]['positive'] or reg['positive']
                    result[sub][modkey]['negative'] = \
                        result[sub][modkey]['negative'] or reg['negative']
                    if len(reg['comments']):
                        result[sub][modkey]['comments'].append(reg['comments'])


    return result


def regsites_tab(regsites, outfile = None):
    """
    Exports PhosphoSite regulatory sites as a tabular file, all
    IDs translated to UniProt.
    """

    header = [
        'uniprot_a', 'isoform_a', 'a_res_aa', 'a_res_num', 'a_mod_type',
        'effect', 'uniprot_b', 'references'
    ]
    result = []
    for uniprot, regsite in iteritems(regsites):
        isoform = '1'
        uniprot = uniprot.split('-')
        if len(uniprot) > 1:
            isoform = uniprot[1]
        uniprot = uniprot[0]
        for r in regsite:
            if r['organism'] == 'human':
                for i in r['induces']:
                    other = mapping.map_name(i, 'genesymbol', 'uniprot')
                    for o in other:
                        if o != 'unmapped':
                            result.append([
                                uniprot, isoform, r['aa'], r['res'], r['modt'],
                                '+', o
                            ])
                for i in r['disrupts']:
                    other = mapping.map_name(i, 'genesymbol', 'uniprot')
                    for o in other:
                        if o != 'unmapped':
                            result.append([
                                uniprot, isoform, r['aa'], r['res'], r['modt'],
                                '-', o, ';'.join(r['pmids'])
                            ])
    if outfile is not None:
        out = '\t'.join(header) + '\n'
        for r in result:
            out += '\t'.join(r) + '\n'
        with open(outfile, 'w') as f:
            f.write(out)
    return result



def phosphosite_interactions(cache = True, ncbi_tax_id = 9606):
    """
    Downloads curated and HTP data from Phosphosite,
    from preprocessed cache file if available.
    Processes BioPAX format.
    Returns list of interactions.
    """

    curated_cache = urls.files['phosphosite']['curated']
    noref_cache = urls.files['phosphosite']['noref']

    if cache and os.path.exists(curated_cache) and os.path.exists(noref_cache):

        return (
            pickle.load(open(curated_cache, 'rb')),
            pickle.load(open(noref_cache, 'rb')),
        )

    result_curated = []
    result_noref = []
    url = urls.urls['psite_bp']['url']
    c = curl.Curl(url, silent = False, large = True)
    bpax = c.gzfile
    xml = ET.parse(bpax)
    xmlroot = xml.getroot()
    bpprefix = '{http://www.biopax.org/release/biopax-level3.owl#}'
    rdfprefix = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    proteins = {}
    for p in xmlroot.iter(bpprefix + 'ProteinReference'):
        psid = p.attrib[rdfprefix + 'ID']
        db = p.find(bpprefix + 'xref').find(bpprefix + 'UnificationXref').find(
            bpprefix + 'db').text
        up = p.find(bpprefix + 'xref').find(bpprefix + 'UnificationXref').find(
            bpprefix + 'id').text
        tax = ''
        if p.find(bpprefix + 'organism') is not None:
            tmp = p.find(bpprefix + 'organism')
            if rdfprefix + 'resource' in tmp.attrib:
                tax = tmp.attrib[rdfprefix + 'resource'].split('_')[1]
        if db == 'UniProtKB':
            up = up[0:6]
        proteins[psid] = {'id': up, 'db': db, 'species': tax, 'psid': psid}
    evidences = {}
    for p in xmlroot.iter(bpprefix + 'EvidenceCodeVocabulary'):
        evid = p.attrib[rdfprefix + 'ID'].split('_')[1]
        evname = p.find(bpprefix + 'term').text
        evidences[evid] = evname
    ev_short = {'0113': 'WB', '0427': 'MS', '0074': 'MA', '0421': 'AB'}
    nosrc = []
    notgt = []
    norefs = []
    noev = []
    noth = []
    edges = []

    for c in xmlroot.findall(bpprefix + 'Catalysis'):
        if rdfprefix + 'resource' in c.find(bpprefix + 'controller').attrib:
            src = 'po_' + \
                c.find(
                    bpprefix + 'controller').attrib[rdfprefix + 'resource'].split('_')[1]
        else:
            srcProt = c.find(bpprefix + 'controller').find(bpprefix +
                                                           'Protein')
            if srcProt is not None:
                src = 'po_' + srcProt.attrib[rdfprefix + 'ID'].split('_')[1]
            else:
                nosrc.append(c)
        tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                       'ProteinReference')
        tgt = next(tgtProt, None)
        if tgt is not None:
            tgt = tgt.attrib[rdfprefix + 'ID']
        else:
            tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                           'entityReference')
            tgt = next(tgtProt, None)
            if tgt is not None:
                if rdfprefix + 'resource' in tgt.attrib:
                    tgt = tgt.attrib[rdfprefix + 'resource'][1:]
            else:
                tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                               'left')
                tgt = next(tgtProt, None)
                if tgt is not None:
                    if rdfprefix + 'resource' in tgt.attrib:
                        tgt = 'po_' + \
                            tgt.attrib[rdfprefix + 'resource'].split('_')[1]
                else:
                    notgt.append(c)
        refs = c.iter(bpprefix + 'PublicationXref')
        pmids = []
        for r in refs:
            pm = r.attrib[rdfprefix + 'ID'].split('_')
            if pm[0] == 'pmid':
                pmids.append(pm[1])
        refs = c.iter(bpprefix + 'evidence')
        for r in refs:
            rrefs = r.iter(bpprefix + 'xref')
            for rr in rrefs:
                if rdfprefix + 'resource' in rr.attrib:
                    pm = rr.attrib[rdfprefix + 'resource'].split('_')
                    if pm[0] == 'pubmed':
                        pmids.append(pm[1])
        evs = []
        for e in c.iter(bpprefix + 'evidenceCode'):
            if rdfprefix + 'resource' in e.attrib:
                evs.append(ev_short[e.attrib[rdfprefix + 'resource'].split('_')
                                    [1]])
            else:
                ev = e.find(bpprefix + 'EvidenceCodeVocabulary')
                evs.append(ev_short[ev.attrib[rdfprefix + 'ID'].split('_')[1]])
        for e in c.iter(bpprefix + 'evidence'):
            if rdfprefix + 'resource' in e.attrib:
                ev = e.attrib[rdfprefix + 'resource'].split('_')
                if len(ev) == 4:
                    if len(ev[3]) == 4:
                        evs.append(ev_short[ev[3]])
        if (src is not None and tgt is not None and src in proteins and
                tgt in proteins and proteins[src]['id'] is not None and
                proteins[tgt]['id'] is not None):
            edges.append({
                'src': proteins[src],
                'tgt': proteins[tgt],
                'pmids': list(set(pmids)),
                'evs': list(set(evs))
            })
            if len(evs) == 0:
                noev.append(c)
            if len(pmids) == 0:
                norefs.append(c)
            if len(evs) == 0 and len(pmids) == 0:
                noth.append(c)

    if ncbi_tax_id:

        all_uniprots = uniprot_db.all_uniprots(organism = ncbi_tax_id)

    for e in edges:

        if (
            ncbi_tax_id and (
                e['src']['id'] not in all_uniprots or
                e['tgt']['id'] not in all_uniprots
            )
        ):

            continue

        this_iaction = [
            e['src']['id'], e['tgt']['id'], e['src']['species'],
            e['tgt']['species'], ';'.join(e['evs']), ';'.join(e['pmids'])
        ]

        if len(this_iaction[-1]) > 0:

            result_curated.append(this_iaction)

        else:

            result_noref.append(this_iaction)

    pickle.dump(result_curated, open(curated_cache, 'wb'))
    pickle.dump(result_noref, open(noref_cache, 'wb'))
    return result_curated, result_noref


def phosphosite_interactions_new(cache = True):
    """
    Downloads curated and HTP data from Phosphosite,
    from preprocessed cache file if available.
    Processes BioPAX format.
    Returns list of interactions.
    """

    curated_cache = urls.files['phosphosite']['curated']
    noref_cache = urls.files['phosphosite']['noref']

    if (
        cache and
        os.path.exists(curated_cache) and
        os.path.exists(noref_cache)
    ):

        with open(curated_cache, 'rb') as fp:

            data_curated = pickle.load(fp)

        with open(noref_cache, 'rb') as fp:

            data_noref = pickle.load(fp)

        return data_curated, data_noref


    def collect_items(tagname, process_method):

        result = {}

        for p in xmlroot.iter(tagname):

            key, value = process_method(p)
            result[key] = value

        return result


    def process_protein(protein):

        protein_id = protein.attrib['%sID' % rdfprefix]
        database = (
            protein.find(
                '%sxref' % bpprefix
            ).find(
                '%sUnificationXref' % bpprefix
            ).find(
                '%sdb' % bpprefix
            ).text
        )
        identifier = (
            protein.find(
                '%sxref' % bpprefix
            ).find(
                '%sUnificationXref' % bpprefix
            ).find(
                '%sid' % bpprefix
            ).text
        )

        organism = None
        e_organism = protein.find('%sorganism' % bpprefix)
        if (
            e_organism is not None and
            '%sresource' % rdfprefix in e_organism.attrib
        ):
            organism = (
                e_organism.attrib['%sresource' % rdfprefix].split('_')[1]
            )

        return protein_id, (databas, identifier, organism)

    def process_site(site):

        site_id = site.attrib['%sID' % rdfprefix]
        site_offset = site.find('%ssequencePosition').text

        return site_id, site_offset


    def process_modification(seqmodvoc):

        mod_id = seqmodvoc.attrib['%sID' % rdfprefix]
        residue, mod = mod_id.split('_').split('-')

        return mod_id, (residue, mod)


    def get_resource(elem, resource_tag):

        res_attr = '%sresource' % rdfprefix

        if res_attr in elem.attrib:

            return elem.attrib[res_attr][1:]

        else:

            return elem.find(resource_tag).attrib['%sID' % rdfprefix]


    def process_feature(feature):

        feature_id = feature.attrib['%sID' % rdfprefix]
        site = get_resource(
            feature.find('%sfeatureLocation' % bpprefix),
            '%sSequenceSite' % bpprefix,
        )
        modification = get_resource(
            feature.find('%smodificationType' % bpprefix),
            '%sSequenceModificationVocabulary' % bpprefix,
        )

        return feature_id, (site, modification)


    result_curated = []
    result_noref = []
    bpprefix = '{http://www.biopax.org/release/biopax-level3.owl#}'
    rdfprefix = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'

    url = urls.urls['psite_bp']['url']
    c = curl.Curl(url, silent = False, large = True)
    bpax = c.gzfile
    xml = ET.parse(bpax)
    xmlroot = xml.getroot()


    proteins = collect_items(
        '%sProtein' % bpprefix,
        process_method = process_protein,
    )
    sites = collect_items(
        '%sSequenceSite' % bpprefix,
        process_method = process_site,
    )
    modifications = collect_items(
        '%sSequenceModificationVocabulary' % bpprefix,
        process_method = process_modification,
    )
    features = collect_items(
        '%sModificationFeature' % bpprefix,
        process_method = process_feature,
    )

    evidences = {}
    for p in xmlroot.iter(bpprefix + 'EvidenceCodeVocabulary'):
        evid = p.attrib[rdfprefix + 'ID'].split('_')[1]
        evname = p.find(bpprefix + 'term').text
        evidences[evid] = evname
    ev_short = {'0113': 'WB', '0427': 'MS', '0074': 'MA', '0421': 'AB'}
    nosrc = []
    notgt = []
    norefs = []
    noev = []
    noth = []
    edges = []

    for c in xmlroot.findall(bpprefix + 'Catalysis'):
        if rdfprefix + 'resource' in c.find(bpprefix + 'controller').attrib:
            src = 'po_' + \
                c.find(
                    bpprefix + 'controller').attrib[rdfprefix + 'resource'].split('_')[1]
        else:
            srcProt = c.find(bpprefix + 'controller').find(bpprefix +
                                                           'Protein')
            if srcProt is not None:
                src = 'po_' + srcProt.attrib[rdfprefix + 'ID'].split('_')[1]
            else:
                nosrc.append(c)
        tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                       'ProteinReference')
        tgt = next(tgtProt, None)
        if tgt is not None:
            tgt = tgt.attrib[rdfprefix + 'ID']
        else:
            tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                           'entityReference')
            tgt = next(tgtProt, None)
            if tgt is not None:
                if rdfprefix + 'resource' in tgt.attrib:
                    tgt = tgt.attrib[rdfprefix + 'resource'][1:]
            else:
                tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                               'left')
                tgt = next(tgtProt, None)
                if tgt is not None:
                    if rdfprefix + 'resource' in tgt.attrib:
                        tgt = 'po_' + \
                            tgt.attrib[rdfprefix + 'resource'].split('_')[1]
                else:
                    notgt.append(c)
        refs = c.iter(bpprefix + 'PublicationXref')
        pmids = []
        for r in refs:
            pm = r.attrib[rdfprefix + 'ID'].split('_')
            if pm[0] == 'pmid':
                pmids.append(pm[1])
        refs = c.iter(bpprefix + 'evidence')
        for r in refs:
            rrefs = r.iter(bpprefix + 'xref')
            for rr in rrefs:
                if rdfprefix + 'resource' in rr.attrib:
                    pm = rr.attrib[rdfprefix + 'resource'].split('_')
                    if pm[0] == 'pubmed':
                        pmids.append(pm[1])
        evs = []
        for e in c.iter(bpprefix + 'evidenceCode'):
            if rdfprefix + 'resource' in e.attrib:
                evs.append(ev_short[e.attrib[rdfprefix + 'resource'].split('_')
                                    [1]])
            else:
                ev = e.find(bpprefix + 'EvidenceCodeVocabulary')
                evs.append(ev_short[ev.attrib[rdfprefix + 'ID'].split('_')[1]])
        for e in c.iter(bpprefix + 'evidence'):
            if rdfprefix + 'resource' in e.attrib:
                ev = e.attrib[rdfprefix + 'resource'].split('_')
                if len(ev) == 4:
                    if len(ev[3]) == 4:
                        evs.append(ev_short[ev[3]])
        if (src is not None and tgt is not None and src in proteins and
                tgt in proteins and proteins[src]['id'] is not None and
                proteins[tgt]['id'] is not None):
            edges.append({
                'src': proteins[src],
                'tgt': proteins[tgt],
                'pmids': list(set(pmids)),
                'evs': list(set(evs))
            })
            if len(evs) == 0:
                noev.append(c)
            if len(pmids) == 0:
                norefs.append(c)
            if len(evs) == 0 and len(pmids) == 0:
                noth.append(c)
    for e in edges:
        this_iaction = [
            e['src']['id'], e['tgt']['id'], e['src']['species'],
            e['tgt']['species'], ';'.join(e['evs']), ';'.join(e['pmids'])
        ]
        if len(this_iaction[-1]) > 0:
            result_curated.append(this_iaction)
        else:
            result_noref.append(this_iaction)
    pickle.dump(result_curated, open(curated_cache, 'wb'))
    pickle.dump(result_noref, open(noref_cache, 'wb'))
    return result_curated, result_noref


def _phosphosite_filter_organism(psite_data, ncbi_tax_id = 9606):

    all_uniprots = uniprot_db.all_uniprots(organism = ncbi_tax_id)

    return [
        rec
        for rec in psite_data
        if rec[0] in all_uniprots and rec[1] in all_uniprots
    ]


def phosphosite_interactions_curated(ncbi_tax_id = 9606):
    """
    Loads literature curated PhosphoSite data,
    from preprocessed cache file if available.
    Returns list of interactions.
    """

    curated_cache = urls.files['phosphosite']['curated']
    if not os.path.exists(curated_cache):
        curated, noref = phosphosite_interactions(ncbi_tax_id = ncbi_tax_id)
        result = curated
    else:
        result = pickle.load(open(curated_cache, 'rb'))

    return _phosphosite_filter_organism(result, ncbi_tax_id)


def phosphosite_interactions_noref(ncbi_tax_id = 9606):
    """
    Loads HTP PhosphoSite data,
    from preprocessed cache file if available.
    Returns list of interactions.
    """

    noref_cache = urls.files['phosphosite']['noref']
    if not os.path.exists(noref_cache):
        curated, noref = phosphosite_interactions(ncbi_tax_id = ncbi_tax_id)
        result = noref
    else:
        result = pickle.load(open(noref_cache, 'rb'))

    return _phosphosite_filter_organism(result, ncbi_tax_id)


def phosphosite_directions(organism = 'human'):
    """
    From curated and HTP PhosphoSite data generates a
    list of directions.
    """

    curated, noref = phosphosite_interactions()

    return [
        i[:2] for i in curated + noref if i[2] == organism and i[3] == organism
    ]


def phosphosite_interactions_all():

    return phosphosite_interactions_curated() + phosphosite_interactions_noref()
