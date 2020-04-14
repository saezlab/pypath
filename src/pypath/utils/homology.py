#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems
from past.builtins import xrange, range

import os
import sys
import itertools
import importlib as imp
import re
import time
import datetime
import json
import pickle

import timeloop

import pypath.utils.mapping as mapping
import pypath.share.common as common
import pypath.internals.intera as intera
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.uniprot as uniprot_input
import pypath.utils.seq as _se
import pypath.share.session as session_mod
import pypath.utils.taxonomy as taxonomy
import pypath.share.cache as cache_mod

timeloop.app.logging.disable(level = 9999)
_homology_cleanup_timeloop = timeloop.Timeloop()

_logger = session_mod.Logger(name = 'homology')
_log = _logger._log


class HomologyManager(session_mod.Logger):


    def __init__(self, cleanup_period = 10, lifetime = 300):

        session_mod.Logger.__init__(self, name = 'homology')


        @_homology_cleanup_timeloop.job(
            interval = datetime.timedelta(
                seconds = cleanup_period
            )
        )
        def _cleanup():

            self._remove_expired()


        _homology_cleanup_timeloop.start(block = False)

        self.lifetime = lifetime
        self.tables = {}
        self.expiry = {}
        self.cachedir = cache_mod.get_cachedir()

        self._log('HomologyManager has been created.')


    def which_table(self, target, source = 9606, only_swissprot = True):

        key = (source, target, only_swissprot)

        self.expiry[key] = time.time()

        if key not in self.tables:

            self.load(key)

        if key in self.tables:

            return self.tables[key]


    def load(self, key):

        cachefile = common.md5(json.dumps(key))
        cachefile = os.path.join(self.cachedir, cachefile)

        if os.path.exists(cachefile):

            self.tables[key] = pickle.load(open(cachefile, 'rb'))

            self._log(
                'Homology table from taxon %u to %u (only SwissProt: %s) '
                'has been loaded from `%s`.' % (key + (cachefile,))
            )

        else:

            self.tables[key] = self._load(key)
            pickle.dump(self.tables[key], open(cachefile, 'wb'))
            self._log(
                'Homology table from taxon %u to %u (only SwissProt: %s) '
                'has been saved to `%s`.' % (key + (cachefile,))
            )


    def _load(self, key):

        return ProteinHomology(
            target = key[1],
            source = key[0],
            only_swissprot = key[2],
        )


    def translate(
            self,
            source_id,
            target,
            source = 9606,
            only_swissprot = True,
        ):

        table = self.which_table(
            target = target,
            source = source,
            only_swissprot = only_swissprot,
        )

        return table.translate(protein = source_id, source = source)


    def _remove_expired(self):

        for key, last_used in list(self.expiry.items()):

            if time.time() - last_used > self.lifetime and key in self.tables:

                self._log(
                    'Removing homology table from taxon %u to %u '
                    '(only SwissProt: %s)' % key
                )

                del self.tables[key]
                del self.expiry[key]


    def __del__(self):

        if hasattr(_homology_cleanup_timeloop, 'stop'):

            _homology_cleanup_timeloop.stop()


def get_homologene():
    """
    Downloads the recent release of the NCBI HomoloGene database.
    Returns file pointer.
    """

    url = urls.urls['homologene']['url']

    c = curl.Curl(
        url = url,
        silent = False,
        large = True,
        timeout = 1800,
        ignore_content_length = True,
    )

    return c.result


def homologene_dict(source, target, id_type):
    """
    Returns orthology translation table as dict, obtained
    from NCBI HomoloGene data.

    :param int source: NCBI Taxonomy ID of the source species (keys).
    :param int target: NCBI Taxonomy ID of the target species (values).
    :param str id_type: ID type to be used in the dict. Possible values:
        'RefSeq', 'Entrez', 'GI', 'GeneSymbol'.
    """
    ids = {
        'refseq': 5,
        'refseqp': 5,
        'genesymbol': 3,
        'gi': 4,
        'entrez': 2
    }

    try:
        id_col = ids[id_type.lower()]
    except KeyError:
        _log(
            'Unknown ID type: `%s`. Please use RefSeq, '
            'Entrez, GI or GeneSymbol.' % id_type
        )
        raise

    hg = get_homologene()
    hgroup = None
    result = {}

    for l in hg:

        l = l.strip().split('\t')
        this_hgroup = l[0].strip()

        if this_hgroup != hgroup:
            this_source = None
            this_target = None
            hgroup = this_hgroup

        this_taxon = int(l[1].strip())
        if this_taxon == source:
            this_source = l[id_col]
        elif this_taxon == target:
            this_target = l[id_col]

        if this_source is not None and this_target is not None \
            and len(this_source) and len(this_target):
            if this_source not in result:
                result[this_source] = set([])
            result[this_source].add(this_target)

    return result


def homologene_uniprot_dict(source, target, only_swissprot = True):
    """
    Returns orthology translation table as dict from UniProt to Uniprot,
    obtained from NCBI HomoloGene data. Uses RefSeq and Entrez IDs for
    translation.

    :param int source: NCBI Taxonomy ID of the source species (keys).
    :param int target: NCBI Taxonomy ID of the target species (values).
    :param bool only_swissprot: Translate only SwissProt IDs.
    """
    result = {}

    hge = homologene_dict(source, target, 'entrez')
    hgr = homologene_dict(source, target, 'refseq')

    all_source = set(uniprot_input.all_uniprots(
                                                organism = source,
                                                swissprot = 'YES'
                                            ))

    if not only_swissprot:
        all_source_trembl = uniprot_input.all_uniprots(
                                                       organism = source,
                                                       swissprot = 'NO'
                                                   )
        all_source.update(set(all_source_trembl))

    for u in all_source:

        source_e = mapping.map_name(u, 'uniprot', 'entrez', source)
        source_r = mapping.map_name(u, 'uniprot', 'refseqp', source)
        target_u = set([])
        target_r = set([])
        target_e = set([])

        for e in source_e:
            if e in hge:
                target_e.update(hge[e])

        for r in source_r:
            if r in hgr:
                target_r.update(hgr[r])

        for e in target_e:
            target_u.update(
                mapping.map_name(e, 'entrez', 'uniprot', target)
            )

        for r in target_r:
            target_u.update(
                mapping.map_name(e, 'refseqp', 'uniprot', target)
            )


        target_u = \
            itertools.chain(
                *map(
                    lambda tu:
                        mapping.map_name(tu, 'uniprot', 'uniprot', target),
                    target_u
                )
            )

        result[u] = sorted(list(target_u))

    return result


class SequenceContainer(session_mod.Logger):

    def __init__(self, preload_seq = [], isoforms = True):
        """
        This is an object to store sequences of multiple
        organisms and select the appropriate one.
        """

        session_mod.Logger.__init__(self, name = 'homology')

        self.seq_isoforms = isoforms

        for taxon in preload_seq:

            self.load_seq(taxon)

    def load_seq(self, taxon):

        if not hasattr(self, 'seq'):
            self.seq = {}

        taxon = taxon or self.ncbi_tax_id

        if taxon not in self.seq:

            self.seq[taxon] = _se.swissprot_seq(
                organism = taxon,
                isoforms = self.seq_isoforms
            )

    def get_seq(self, protein, taxon = None):

        if taxon is not None:

            if taxon not in self.seq:

                self.load_seq(taxon)

            if protein in self.seq[taxon]:

                return self.seq[taxon][protein]

        else:

            for taxon, seq in iteritems(self.seq):

                if protein in seq:

                    return seq[protein]


class Proteomes(object):

    def __init__(self, preload_prot = [], swissprot_only = True):

        if not hasattr(self, '_taxonomy'):
            self._taxonomy = {}
            self._proteomes = {}

        for taxon in preload_prot:

            self.load_proteome(taxon, swissprot_only)

    def load_proteome(self, taxon, swissprot_only = True):

        key = (taxon, swissprot_only)

        if key not in self._proteomes:

            self._proteomes[key] = (
                set(uniprot_input.all_uniprots(*key))
            )

            for protein in self._proteomes[key]:

                self._taxonomy[protein] = key

            if not swissprot_only:

                self.load_proteome(taxon, True)

    def get_taxon(self, protein, swissprot_only = True):

        if not swissprot_only or self.is_swissprot(protein):

            return self._taxonomy[protein][0]

    def get_taxon_trembl(self, protein):

        if self.has_protein(protein):

            return self._taxonomy[protein][0]

    def has_protein(self, protein):

        return protein in self._taxonomy

    def is_swissprot(self, protein):

        return self.has_protein(protein) and self._taxonomy[protein][1]


class ProteinHomology(Proteomes):

    def __init__(
            self,
            target,
            source = None,
            only_swissprot = True,
        ):
        """
        This class translates between homologous UniProt IDs of
        2 organisms based on NCBI HomoloGene data.
        Uses RefSeq and Entrez IDs for translation.


        :param int target: NCBI Taxonomy ID of the organism
                           to be translated to.
        :param int source: NCBI Taxonomy ID of the default organism
                           to be translated from.
        :param bool only_swissprot: Whether only SwissProt or Trembl IDs
                                    should be used.
        :mapper pypath.mapping.Mapper mapper: A Mapper object.
        """

        self.homo = {}
        self.only_swissprot = only_swissprot
        self.target = target
        self.source = source
        self.set_default_source(source)

        Proteomes.__init__(self)
        self.load_proteome(self.target, self.only_swissprot)

        if source is not None:
            self.homologene_uniprot_dict(source)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def set_default_source(self, source = None):

        self.source = source or self.source


    def get_source(self, source = None):

        source = source or self.source

        if source is None:
            raise ValueError('No source NCBI Taxonomy ID provided.')
        else:
            return source


    def translate(self, protein, source = None):
        """
        For one UniProt ID of the source organism returns all orthologues
        from the target organism.
        """

        if isinstance(protein, common.list_like):

            return list(set(
                itertools.chain(*(
                    self.translate(p, source = source)
                    for p in protein
                ))
            ))

        if self.get_taxon(protein) == self.target:
            return [protein]

        source = self.get_source(source)

        if source not in self.homo:

            self.homologene_uniprot_dict(source)

        if protein in self.homo[source]:

            return self.homo[source][protein]

        else:

            return []


    def homologene_uniprot_dict(self, source):
        """
        Builds orthology translation table as dict from UniProt to Uniprot,
        obtained from NCBI HomoloGene data. Uses RefSeq and Entrez IDs for
        translation.
        """

        source = self.get_source(source)

        self.homo[source] = {}

        hge = homologene_dict(source, self.target, 'entrez')
        hgr = homologene_dict(source, self.target, 'refseq')

        self.load_proteome(source, self.only_swissprot)

        for u in self._proteomes[(source, self.only_swissprot)]:

            source_e = mapping.map_name(
                u, 'uniprot', 'entrez', source)
            source_r = mapping.map_name(
                u, 'uniprot', 'refseqp', source)
            target_u = set([])
            target_r = set([])
            target_e = set([])

            for e in source_e:
                if e in hge:
                    target_e.update(hge[e])

            for r in source_r:
                if r in hgr:
                    target_r.update(hgr[r])

            for e in target_e:
                target_u.update(set(mapping.map_name(
                    e, 'entrez', 'uniprot', self.target)))

            for r in target_r:
                target_u.update(set(mapping.map_name(
                    e, 'refseqp', 'uniprot', self.target)))

            target_u = \
                itertools.chain(
                    *map(
                        lambda tu:
                            mapping.map_name(
                                tu, 'uniprot', 'uniprot', self.target),
                        target_u
                    )
                )

            self.homo[source][u] = sorted(list(target_u))


class PtmHomology(ProteinHomology, SequenceContainer):


    def __init__(self, target, source = None, only_swissprot = True,
             strict = True):

        ProteinHomology.__init__(
            self,
            target,
            source,
            only_swissprot,
        )

        SequenceContainer.__init__(self)
        self.load_seq(taxon = self.target)

        self.reptm = re.compile(r'([A-Z\d]{6,10})_([A-Z])(\d*)')

        self.strict = strict

        self.ptm_orthology()


    def translate_site(
            self,
            protein,
            res,
            offset,
            isoform = 1,
            typ = 'phosphorylation',
            source_taxon = None,
        ):
        """
        Translates one PTM site.
        """

        result = set()

        self.set_default_source(source_taxon)

        source = self.get_source(source_taxon)

        sourceptm = (protein, isoform, res, offset, source, typ)

        if self.get_taxon(protein.identifier) == self.target:
            result.add(sourceptm)
            return result

        if sourceptm in self.ptmhomo:

            if self.target in self.ptmhomo[sourceptm]:

                result = self.ptmhomo[sourceptm]

        if not result and not self.strict:

            tsubs = ProteinHomology.translate(
                self,
                protein.identifier,
                source = source,
            )

            for tsub in tsubs:

                se = self.get_seq(tsub)

                if se is None:
                    continue

                for toffset in xrange(offset, offset + 3):

                    for i in se.isoforms():

                        tres = se.get(toffset, isoform = i)

                        if tres == res:

                            result.add((
                                tsub,
                                i,
                                tres,
                                toffset,
                                self.target,
                                typ
                            ))

                    if result:
                        break

        return result


    def translate_domain(self, domain):

        return (
            list(
                map(
                    lambda x:
                        intera.Domain(
                            protein = x,
                            ncbi_tax_id = self.target,
                        ),
                    ProteinHomology.translate(
                        self,
                        domain.protein.identifier,
                        source = self.get_source()
                    )
                )
            )
        )


    def translate_ptm(self, ptm):

        tptms = self.translate_site(
            ptm.protein,
            ptm.residue.name,
            ptm.residue.number,
            ptm.residue.isoform,
            ptm.typ,
        )

        result = []

        for x in tptms:

            se = self.get_seq(x[0])

            if (se is None or x[1] not in se.isof) and self.strict:
                continue

            res = intera.Residue(
                number = x[3],
                name = x[2],
                protein = x[0],
                isoform = x[1],
                ncbi_tax_id = self.target,
            )
            start, end, region = (
                se.get_region(x[3], isoform = x[1])
                if se is not None and x[1] in se.isof
                else (None, None, None)
            )
            mot = intera.Motif(
                protein = x[0],
                start = start,
                end = end,
                instance = region,
                isoform = x[1],
                ncbi_tax_id = self.target,
            )

            ptm = intera.Ptm(
                protein = x[0],
                motif = mot,
                residue = res,
                typ = x[5],
                isoform = x[1],
                evidences = ptm.evidences,
                ncbi_tax_id = self.target,
            )

            result.append(ptm)

        return result


    def translate_domain_motif(self, dmotif):

        ds = self.translate_domain(dmotif.domain)
        ps = self.translate_ptm(dmotif.ptm)

        return (
            list(
                map(
                    lambda x:
                        intera.DomainMotif(
                            x[0],
                            x[1],
                            evidences = dmotif.evidences,
                        ),
                        itertools.product(ds, ps)
                )
            )
        )


    def translate_residue(self, residue):

        return (
            list(
                map(
                    lambda r:
                        intera.Residue(r[3], r[2], r[0], isoform = r[1]),
                    self.translate_site(
                        residue.protein,
                        residue.name,
                        residue.number,
                        residue.isoform
                    )
                )
            )
        )


    def translate(self, x, return_strings = False, **kwargs):
        """
        Translates anything:

        - one UniProt ID
        - one PTM provided as tuple of (UniProt, amino acid, offest)
        - one PTM provided as string (e.g. `P00533_S231`)
        - instance from pypath.intera: DomainMotif, Domain or Ptm

        Additional arguments can be isoform and typ (modification type).

        """

        result = []

        if type(x) is tuple:

            result = self.translate_site(*x, **kwargs)

        elif type(x) in common.char_types:

            ptm = self.reptm.match(x)

            if ptm is not None:

                result = self.translate_site(ptm[0], ptm[1],
                                             int(ptm[2]), **kwargs)

        if return_strings:

            result = list(map(lambda r:
                              '%s_%s%u' % (r[0], r[2], r[3]),
                              result))

        elif isinstance(x, intera.Ptm):

            result = self.translate_ptm(x)

        elif isinstance(x, intera.Domain):

            result = self.translate_domain(x)

        elif isinstance(x, intera.DomainMotif):

            result = self.translate_domain_motif(x)

        return result


    def ptm_orthology(self):
        """
        Creates an orthology translation dict of phosphosites
        based on phosphorylation sites table from PhosphoSitePlus.
        In the result all PTMs represented by a tuple of the following
        6 elements: UniProt ID, isoform (int), residue one letter code,
        residue number (int), NCBI Taxonomy ID (int), modification type.

        """

        self.ptmhomo = {}

        nondigit = re.compile(r'[^\d]+')

        unknown_taxa = set([])

        for typ in common.psite_mod_types:

            groups = {}

            url = urls.urls['psite_%s' % typ[0]]['url']
            c = curl.Curl(url, silent=False, large=True)

            data = c.result

            for _ in xrange(4):
                null = next(data)

            for r in data:

                r = r.split('\t')

                if len(r) < 10:

                    continue

                uniprot = r[2]
                isoform = (
                    1
                        if '-' not in uniprot else
                    int(uniprot.split('-')[1])
                )
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

                        if site1 not in self.ptmhomo:

                            self.ptmhomo[site1] = {}

                        if site2[4] not in self.ptmhomo[site1]:

                            self.ptmhomo[site1][site2[4]] = set([])

                        self.ptmhomo[site1][site2[4]].add(site2)

        if len(unknown_taxa):
            self._log(
                'Unknown taxa encountered: %s' % (
                    ', '.join(sorted(unknown_taxa))
                )
            )


def init():

    globals()['manager'] = HomologyManager()


def get_manager():

    if 'manager' not in globals():

        init()

    return globals()['manager']


def translate(source_id, target, source = 9606, only_swissprot = True):

    manager = get_manager()

    return manager.translate(
        source_id  = source_id,
        target = target,
        source = source,
        only_swissprot = only_swissprot,
    )
