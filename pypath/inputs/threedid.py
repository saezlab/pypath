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

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.inputs.pfam as pfam_input
import pypath.inputs.pdb as pdb_input
import pypath.inputs.uniprot_db as uniprot_db
import pypath.utils.pdb as pdb_utils
import pypath.share.cache as cache
import pypath.internals.intera as intera
import pypath.share.progress as progress


def get_3did_ddi(residues = False, ddi_flat = None, organism = 9606):

    if ddi_flat is None:
        c = curl.Curl(urls.urls['3did_ddi']['url'], silent = False)
        data = c.result
        tmpfile = '3did_flat_tmp'

        if data is None:
            return None

        with open(tmpfile, 'w') as f:
            f.write(data)

        lnum = data.count('\n')
        del data

    else:
        tmpfile = ddi_flat

    u_pfam, pfam_u = pfam_input.pfam_uniprot(organism = organism)
    u_pdb, pdb_u = pdb_input.pdb_chains()

    if pfam_u is None or pdb_u is None:

        return None

    ddi = {}
    interfaces = {}
    pdblist = {}
    ddi_collect = False
    con_collect = False

    with open(tmpfile, 'r') as f:

        prg = progress.Progress(lnum, 'Reading data', 33)

        for l in f:
            prg.step()

            if l.startswith('#=') and con_collect:

                interfaces[(uniprot1, uniprot2, pdb)].append(this_interface)
                con_collect = False

            if l.startswith('#=ID'):
                # new domain pair: attach previous to results:
                if ddi_collect:
                    for u1 in uniprots1:
                        for u2 in uniprots2:
                            if u1 != u2 and len(pdblist) > 0:
                                if (u1, u2) not in ddi:
                                    ddi[(u1, u2)] = {}

                                if (pfam1, pfam2) not in ddi[(u1, u2)]:
                                    ddi[(u1, u2)][(pfam1, pfam2)] = {
                                        'pdbs': pdblist
                                    }

                    ddi_collect = False

                pdblist = {}
                l = l.split('\t')
                pfam1 = l[3].split('(')[1].split('.')[0]
                pfam2 = l[4].split('.')[0]
                uniprots1 = [] if pfam1 not in pfam_u else pfam_u[pfam1]
                uniprots2 = [] if pfam2 not in pfam_u else pfam_u[pfam2]

                if len(set(uniprots1) | set(uniprots2)) > 1:

                    ddi_collect = True

            elif l.startswith('#=3D'):

                l = l.split('\t')
                pdb = l[1]
                chain1 = l[2].split(':')[0]
                chain2 = l[3].split(':')[0]

                if (
                    pdb in pdb_u and
                    chain1 in pdb_u[pdb] and
                    chain2 in pdb_u[pdb]
                ):

                    uniprot1 = pdb_u[pdb][chain1]['uniprot']
                    uniprot2 = pdb_u[pdb][chain2]['uniprot']

                    if uniprot1 != uniprot2:

                        if pdb not in pdblist:

                            pdblist[pdb] = []

                        pdblist[pdb] = common.add_to_list(
                            pdblist[pdb],
                            (uniprot1, uniprot2),
                        )

                    if residues:

                        if chain1 != chain2:

                            if (
                                pdb_u[pdb][chain1]['offset'] is not None and
                                pdb_u[pdb][chain2]['offset'] is not None and
                                pdb_u[pdb][chain1]['uniprot'] !=
                                pdb_u[pdb][chain2]['uniprot']
                            ):

                                con_collect = True
                                offset1 = pdb_u[pdb][chain1]['offset']
                                offset2 = pdb_u[pdb][chain2]['offset']

                                this_interface = intera.Interface(
                                    uniprot1,
                                    uniprot2,
                                    source = '3DID',
                                    pdb = pdb,
                                )

                                key = (uniprot1, uniprot2, pdb)

                                if key not in interfaces:

                                    interfaces[key] = []

                            else:

                                con_collect = False

            elif not residues or not con_collect:

                continue

            else:

                l = l.split('\t')

                if len(l) > 3:

                    rnum1 = int(common.non_digit.sub('', l[2])) + offset1
                    rnum2 = int(common.non_digit.sub('', l[3])) + offset2

                    this_interface.add_residues(
                        (rnum1, l[0], uniprot1),
                        (rnum2, l[1], uniprot2),
                    )

        prg.terminate()
        prg = progress.Progress(len(ddi), 'Processing interfaces', 99)

        if residues:

            for u, v1 in iteritems(ddi):

                prg.step()

                for d, v2 in iteritems(v1):
                    for p in v2['pdbs'].keys():

                        key = (u[0], u[1], p)

                        if key in interfaces:

                            ddi[u][d]['interfaces'] = interfaces[key]

        prg.terminate()

    if ddi_flat is None:

        os.remove(tmpfile)

    if residues:

        return ddi, interfaces

    else:

        return ddi


def get_3did(ddi_flat = None, res = True, organism = 9606, pickl = True):

    resultfile = os.path.join(cache.get_cachedir(), '3did_ddi.pickle')

    if pickl and os.path.exists(resultfile):

        result = pickle.load(open(resultfile, 'rb'))

        if len(result) == 1:
            return result

        else:
            return result[0], result[1]

    if ddi_flat is None:

        c = curl.Curl(urls.urls['3did_ddi']['url'], silent = False)
        data = c.result
        tmpfile = '3did_flat_tmp'

        if data is None:

            return None

        with open(tmpfile, 'w') as f:

            f.write(data)

        lnum = data.count('\n')
        del data

    elif os.path.exists(ddi_flat):

        tmpfile = ddi_flat

    else:

        return None

    u_pdb, pdb_u = pdb_input.pdb_chains()
    all_unip = set(uniprot_db.all_uniprots(organism = organism))

    if all_unip is None or pdb_u is None:

        return None

    ddi = []
    interfaces = []
    pdb = pdb_prev = intf = None
    skip = True
    rmap = pdb_utils.ResidueMapper()

    with open(tmpfile, 'r') as f:

        prg = progress.Progress(
            lnum,
            'Processing 3DID domain-domain interactions',
            33,
        )

        for l in f:

            prg.step()
            l = l.split('\t')

            if l[0].startswith('#=ID'):

                pfam1 = l[3].split('.')[0][2:]
                pfam2 = l[4].split('.')[0]

            elif l[0].startswith('#=3D'):

                pdb_prev = pdb
                skip = True
                pdb = l[1]
                chain1 = l[2][0]
                chain2 = l[3][0]
                uniprot1 = uniprot2 = None

                if pdb != pdb_prev:

                    rmap.clean()

                if pdb in pdb_u:

                    if chain1 in pdb_u[pdb]:

                        uniprot1 = pdb_u[pdb][chain1]['uniprot']

                    if chain2 in pdb_u[pdb]:

                        uniprot2 = pdb_u[pdb][chain2]['uniprot']

                if (
                    uniprot1 is not None and
                    uniprot2 is not None and
                    uniprot1 in all_unip and
                    uniprot2 in all_unip and
                    uniprot1 != uniprot2
                ):

                    skip = False

                    if intf is not None:

                        interfaces.append(intf)

                    intf = intera.Interface(uniprot1, uniprot2, '3DID', pdb)
                    u1start = u1end = u2start = u2end = {}

                    if l[2].count('-') == 1:

                        start1 = int(
                            common.non_digit.sub('', l[2][2:].split('-')[0])
                        )
                        end1 = int(
                            common.non_digit.sub('', l[2][2:].split('-')[1])
                        )
                        u1start = rmap.get_residue(
                            pdb,
                            start1,
                            chain = chain1,
                        )
                        u1end = rmap.get_residue(pdb, end1, chain = chain1)

                    if l[3].count('-') == 1:

                        start2 = int(
                            common.non_digit.sub('', l[3][2:].split('-')[0])
                        )
                        end2 = int(
                            common.non_digit.sub('', l[3][2:].split('-')[1])
                        )
                        u2start = rmap.get_residue(
                            pdb,
                            start2,
                            chain = chain2,
                        )
                        u2end = rmap.get_residue(pdb, end2, chain = chain2)

                    u1start = u1start.resnum if u1start else None
                    u1end = u1end.resnum if u1end else None
                    u2start = u2start.resnum if u2start else None
                    u2end = u2end.resnum if u2end else None

                    dom1 = intera.Domain(
                        uniprot1,
                        domain = pfam1,
                        start = u1start,
                        end = u1end,
                        isoform = 1,
                    )

                    dom2 = intera.Domain(
                        uniprot2,
                        domain = pfam2,
                        start = u2start,
                        end = u2end,
                        isoform = 1,
                    )

                    dd = intera.DomainDomain(dom1, dom2, [pdb], '3DID')
                    ddi.append(dd)

            elif not skip and res and not l[0].startswith('//'):

                conv1 = rmap.get_residue(
                    pdb,
                    int(common.non_digit.sub('', l[2])),
                    chain = chain1,
                )
                conv2 = rmap.get_residue(
                    pdb,
                    int(common.non_digit.sub('', l[3])),
                    chain = chain2,
                )

                if conv1 and conv2:

                    intf.add_residues(
                        (conv1.resnum, l[0], uniprot1),
                        (conv2.resnum, l[1], uniprot2),
                    )

        interfaces.append(intf)
        prg.terminate()

    if ddi_flat is None:

        os.remove(tmpfile)

    if res:

        pickle.dump([ddi, interfaces], open(resultfile, 'wb'))

        return ddi, interfaces

    else:

        pickle.dump([ddi], open(resultfile, 'wb'))

        return ddi


def get_3did_dmi(dmi_flat = None):

    resultfile = os.path.join(cache.get_cachedir(), '3did_dmi.pickle')

    if os.path.exists(resultfile):

        return pickle.load(open(resultfile, 'rb'))

    if dmi_flat is None:

        c = curl.Curl(urls.urls['3did_dmi']['url'], silent = False)
        data = c.result
        tmpfile = '3did_dmi_flat_tmp'

        if data is None:
            return None

        with open(tmpfile, 'w') as f:
            f.write(data)

        lnum = data.count('\n')
        del data

    elif os.path.exists(dmi_flat):

        tmpfile = dmi_flat

    else:

        return None

    u_pdb, pdb_u = pdb_input.pdb_chains()

    if pdb_u is None:

        return None

    dmi = {}
    rmap = pdb_utils.ResidueMapper()

    with open(tmpfile, 'r') as f:

        prg = progress.Progress(
            lnum,
            'Processing 3DID domain-motif interactions',
            1,
        )

        for l in f:

            prg.step()

            l = l.strip().split()

            if l[0].startswith('#=ID'):

                domain = l[3]

            if l[0].startswith('#=PT'):

                regex = l[1]

            if l[0].startswith('#=3D'):

                pdb = l[1]
                chain1 = l[2].split(':')[0]
                chain2 = l[3].split(':')[0]

                if l[2].count('-') == 1 and l[3].count('-') == 1:

                    pdb_region1 = [
                        int(common.non_digit.sub('', x))
                        for x in l[2].split(':')[1].split('-')
                    ]
                    pdb_region2 = [
                        int(common.non_digit.sub('', x))
                        for x in l[3].split(':')[1].split('-')
                    ]
                    u1start = rmap.get_residue(
                        pdb,
                        pdb_region1[0],
                        chain = chain1,
                    )
                    u1end = rmap.get_residue(
                        pdb,
                        pdb_region1[1],
                        chain = chain1,
                    )
                    u2start = rmap.get_residue(
                        pdb,
                        pdb_region2[0],
                        chain = chain2,
                    )
                    u2end = rmap.get_residue(
                        pdb,
                        pdb_region2[1],
                        chain = chain2,
                    )

                    if u1start and u2start and u1end and u2end:

                        uniprot_key = (
                            u1start[chain1]['uniprot'],
                            u2start[chain2]['uniprot'],
                        )

                        residue_key = (
                            u1start[chain1]['resnum'],
                            u1end[chain1]['resnum'],
                            u2start[chain2]['resnum'],
                            u2end[chain2]['resnum'],
                        )

                        if uniprot_key not in dmi:

                            dmi[uniprot_key] = {}

                        if residue_key not in dmi[uniprot_key]:

                            dmi[uniprot_key][residue_key] = []

                        dmi[uniprot_key][residue_key].append({
                            'pdb': pdb,
                            'regex': regex,
                            'instance': l[4],
                            'domain': domain,
                            'contacts': int(non_digit.sub('', l[5])),
                            'topology': int(non_digit.sub('', l[6])),
                        })

        prg.terminate()

    if dmi_flat is None:

        os.remove(tmpfile)

    pickle.dump(dmi, open(resultfile, 'wb'))

    return dmi


def process_3did_dmi():

    dmi = get_3did_dmi()

    if dmi is None:

        return None

    dname_pfam, pfam_dname = pfam_input.pfam_names()
    dname_re = re.compile(r'(.*)(_[A-Z]{3}_)(.*)')
    dmi2 = {}
    prg = progress.Progress(len(dmi), 'Processing data', 11)

    for uniprots, dmis in iteritems(dmi):
        prg.step()

        if uniprots not in dmi2:
            dmi2[uniprots] = []

        for regions, dmi_list in iteritems(dmis):
            new = True

            for dm in dmi_list:
                if new:
                    pfam = None
                    dname = None
                    mname = None
                    name_match = dname_re.match(dm['domain'])

                    if name_match:
                        dname = name_match.groups(0)[0]
                        mname = ''.join(name_match.groups(0)[1:])[1:]

                    if dname in dname_pfam:
                        pfam = dname_pfam[dname][0]

                    domain = pfam if pfam is not None else dname
                    domain_name = 'pfam' if pfam is not None else 'domain_name'
                    dom = intera.Domain(
                        uniprots[0],
                        domain = domain,
                        domain_id_type = domain_name,
                        start = regions[0],
                        end = regions[1])
                    mot = intera.Motif(
                        uniprots[1],
                        regions[2],
                        regions[3],
                        instance = dm['instance'],
                        regex = dm['regex'],
                        motif_name = mname)
                    ptm = intera.Ptm(uniprots[1], motif = mot, source = '3DID')
                    dommot = intera.DomainMotif(dom, ptm, sources = '3DID')
                    new = False

                dommot.add_pdbs(dm['pdb'])

            dmi2[uniprots].append(dommot)

    prg.terminate()

    return dmi2
