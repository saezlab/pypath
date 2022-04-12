#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
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
import collections
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
import pypath.inputs.homologene as homologene_input
import pypath.inputs.biomart as biomart
import pypath.utils.seq as _se
import pypath.share.session as session
import pypath.share.settings as settings
import pypath.utils.taxonomy as taxonomy
import pypath.share.cache as cache_mod

_homology_cleanup_timeloop = timeloop.Timeloop()
_homology_cleanup_timeloop.logger.setLevel(9999)

_logger = session.Logger(name = 'homology')
_log = _logger._log


class Ortholog(
        collections.namedtuple(
            'OrthologBase',
            (
                'uniprot',
                'resource',
                'ensembl_hc',
                'ensembl_type',
            ),
        )
    ):

    def __new__(
            cls,
            uniprot,
            resource,
            ensembl_hc = None,
            ensembl_type = None
        ):

        return super(Ortholog, cls).__new__(
            cls,
            uniprot,
            resource,
            ensembl_hc = ensembl_hc,
            ensembl_type = ensembl_type,
        )


    def __str__(self):

        return self.uniprot


    def __repr__(self):

        return '<Ortholog %s (%s)>' % (self.uniprot, self.resource)


    def __eq__(self, other):

        return self.__str__() == other.__str__()


    def __hash__(self):

        return super().__hash__()


class HomologyManager(session.Logger):


    def __init__(
            self,
            cleanup_period = 10,
            lifetime = 300,
            homologene = None,
            ensembl = None,
            ensembl_hc = None,
            ensembl_types = None,
        ):

        session.Logger.__init__(self, name = 'homology')


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

        for param in ('homologene', 'ensembl', 'ensembl_hc', 'ensembl_types'):

            setattr(
                self,
                param,
                settings.get('homology_%s' % param)
                    if locals()[param] is None else
                locals()[param]
            )

        self.ensembl_types = common.to_set(self.ensembl_types)
        self.ensembl_types = {
            x if x.startswith('ortholog_') else ('ortholog_%s' % x)
            for x in self.ensembl_types
        }

        self._log('HomologyManager has been created.')


    def which_table(self, target, source = 9606, only_swissprot = True):

        key = (source, target, only_swissprot)

        self.expiry[key] = time.time()

        if key not in self.tables:

            self.load(key)

        if key in self.tables:

            return self.tables[key]


    def load(self, key):

        cachefile = '%s-homology-%u-%u.pickle' % (
            common.md5(json.dumps(key)),
            key[0],
            key[1],
        )
        cachefile = os.path.join(self.cachedir, cachefile)

        if os.path.exists(cachefile):

            self.tables[key] = pickle.load(open(cachefile, 'rb'))
            self.tables[key].homologene = self.homologene
            self.tables[key].ensembl = self.ensembl
            self.tables[key].ensembl_hc = self.ensembl_hc
            self.tables[key].ensembl_types = self.ensembl_types

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
            homologene = self.homologene,
            ensembl = self.ensembl,
            ensembl_hc = self.ensembl_hc,
            ensembl_types = self.ensembl_types,
        )


    def translate(
            self,
            source_id,
            target,
            source = 9606,
            only_swissprot = True,
            homologene = None,
            ensembl = None,
            ensembl_hc = None,
            ensembl_types = None,
        ):
        """
        For one or more UniProt ID of the source organism returns all
        orthologs from the target organism.

        Args:
            source_id (str,list): UniProt ID of one or more protein in the
                source organism.
            target (int,str): The target organism.
            source (int,str): The source organism.
            homologene (bool): Use NCBI HomoloGene data for ortholog lookup.
            ensembl (bool): Use Ensembl data for ortholog lookup.
            ensembl_hc (bool): Use only high confidence orthology relations
                from Ensembl. By default it is True. You can also set it
                by the `ensembl_hc` attribute.
            ensembl_types (list): The Ensembl orthology relationship types
                to use. Possible values are `one2one`, `one2many` and
                `many2many`. By default only `one2one` is used. You can
                also set this parameter by the `ensembl_types` attribute.

        Returns:
            Set of UniProt IDs of homologous proteins in the target taxon.
        """

        table = self.which_table(
            target = target,
            source = source,
            only_swissprot = only_swissprot,
        )

        homologene = self.homologene if homologene is None else homologene
        ensembl = self.ensembl if ensembl is None else ensembl
        ensembl_hc = self.ensembl_hc if ensembl_hc is None else ensembl_hc
        ensembl_types = (
            self.ensembl_types if ensembl_types is None else ensembl_types
        )
        ensembl_types = common.to_set(ensembl_types)
        ensembl_types = {
            x if x.startswith('ortholog_') else ('ortholog_%s' % x)
            for x in ensembl_types
        }

        return table.translate(
            protein = source_id,
            source = source,
            homologene = homologene,
            ensembl = ensembl,
            ensembl_hc = ensembl_hc,
            ensembl_types = ensembl_types,
        )


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


class SequenceContainer(session.Logger):


    def __init__(self, preload_seq = [], isoforms = True):
        """
        This is an object to store sequences of multiple
        organisms and select the appropriate one.
        """

        if not hasattr(self, '_logger'):

            session.Logger.__init__(self, name = 'homology')

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
            homologene = True,
            ensembl = True,
            ensembl_hc = True,
            ensembl_types = None,
        ):
        """
        This class translates between homologous UniProt IDs of two organisms
        based on NCBI HomoloGene and Ensembl data. In case of HomoloGene,
        the UniProt-UniProt translation table is created by translating the
        source organism UniProts to RefSeq and Entrez IDs, finding the
        homologues (orthologues) for these IDs, and then translating them
        to the target organism UniProt IDs. In case of Ensembl, we obtain
        data with Ensembl protein identifiers and translate those to UniProt.

        Args:
            target (int): NCBI Taxonomy ID of the target organism.
            source (int): NCBI Taxonomy ID of the default source organism.
                Multiple source organisms can be used on the same instance.
            only_swissprot (bool): Use only SwissProt IDs.
            homologene (bool): Use homology information from NCBI HomoloGene.
            ensembl (bool): Use homology information from Ensembl.
            ensembl_hc (bool): Use only the high confidence
                orthology relations from Ensembl.
            ensembl_types (list): Ensembl orthology relation types to use.
                Possible values are `one2one`, `one2many` and `many2many`.
                By default only `one2one` is used.
        """

        self.orthologs = {}
        self.only_swissprot = only_swissprot
        self.target = taxonomy.ensure_ncbi_tax_id(target)
        self.source = source
        self._default_source = 9606
        self.set_default_source(source)
        for param in ('homologene', 'ensembl', 'ensembl_hc', 'ensembl_types'):

            setattr(
                self,
                param,
                settings.get('homology_%s' % param)
                    if locals()[param] is None else
                locals()[param]
            )

        self.ensembl_types = common.to_set(self.ensembl_types)
        self.ensembl_types = {
            x if x.startswith('ortholog_') else ('ortholog_%s' % x)
            for x in self.ensembl_types
        }

        Proteomes.__init__(self)
        self.load_proteome(self.target, self.only_swissprot)

        if source is not None:

            self.load(source)


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def load(self, source = None):

        if self.homologene:

            self.load_homologene(source)

        if self.ensembl:

            self.load_ensembl(source)


    def ensure_source_taxon(self, source):

        if source not in self.orthologs:

            self.load(source = source)


    def set_default_source(self, source = None):

        self.source = self.get_source(source)


    def get_source(self, source = None):

        source = source or self.source or self._default_source
        _source = taxonomy.ensure_ncbi_tax_id(source)

        if _source is None:

            msg = 'Unknown organism: `%s`.' % str(source)
            _log(msg)
            raise ValueError(msg)

        else:

            return _source


    def translate(
            self,
            protein,
            source = None,
            homologene = None,
            ensembl = None,
            ensembl_hc = None,
            ensembl_types = None,
        ):
        """
        For one UniProt ID of the source organism returns all orthologues
        from the target organism.

        Args:
            protein (str,list): UniProt ID of one or more protein in the
                source organism.
            source (int,str): The source organism.
            homologene (bool): Use NCBI HomoloGene data for ortholog lookup.
            ensembl (bool): Use Ensembl data for ortholog lookup.
            ensembl_hc (bool): Use only high confidence orthology relations
                from Ensembl. By default it is True. You can also set it
                by the `ensembl_hc` attribute.
            ensembl_types (list): The Ensembl orthology relationship types
                to use. Possible values are `one2one`, `one2many` and
                `many2many`. By default only `one2one` is used. You can
                also set this parameter by the `ensembl_types` attribute.

        Returns:
            Set of UniProt IDs of homologous proteins in the target taxon.
        """

        protein = (
            (protein,)
                if hasattr(protein, 'components') else
            common.to_list(protein)
        )

        source = self.get_source(source)

        if source is None:

            msg = (
                'Can not translate without knowing the organism of the input.'
            )
            _log(msg)
            raise ValueError(msg)

        homologene = self.homologene if homologene is None else homologene
        ensembl = self.ensembl if ensembl is None else ensembl
        ensembl_hc = self.ensembl_hc if ensembl_hc is None else ensembl_hc
        ensembl_types = (
            self.ensembl_types if ensembl_types is None else ensembl_types
        )
        ensembl_types = common.to_set(ensembl_types)
        ensembl_types = {
            x if x.startswith('ortholog_') else ('ortholog_%s' % x)
            for x in ensembl_types
        }

        result = {
            p for p in protein
            if self.get_taxon(p) == self.target
        }

        self.ensure_source_taxon(source)

        for p in protein:

            result.update(
                {
                    o.uniprot
                    for o in self.orthologs[source][p]
                    if (
                        (
                            homologene and
                            o.resource == 'HomoloGene'
                        ) or (
                            ensembl and
                            o.resource == 'Ensembl' and (
                                not ensembl_hc or
                                o.ensembl_hc
                            ) and
                            o.ensembl_type in ensembl_types
                        )
                    )
                }
            )

        return result


    def load_homologene(self, source):
        """
        Builds orthology translation table as dict from UniProt to Uniprot,
        obtained from NCBI HomoloGene data. Uses RefSeq and Entrez IDs for
        translation.
        """

        source = self.get_source(source)

        if source not in self.orthologs:

            self.orthologs[source] = collections.defaultdict(set)

        hge = homologene_input.homologene_dict(source, self.target, 'entrez')
        hgr = homologene_input.homologene_dict(source, self.target, 'refseq')

        self.load_proteome(source, self.only_swissprot)

        _log(
            'Loading homology data from NCBI HomoloGene '
            'between organisms `%u` and `%u`' % (source, self.target)
        )

        for u in self._proteomes[(source, self.only_swissprot)]:

            source_e = mapping.map_name(u, 'uniprot', 'entrez', source)
            source_r = mapping.map_name(u, 'uniprot', 'refseqp', source)
            target_u = set()
            target_r = set()
            target_e = set()

            for e in source_e:

                if e in hge:

                    target_e.update(hge[e])

            for r in source_r:

                if r in hgr:

                    target_r.update(hgr[r])

            for e in target_e:

                target_u.update(
                    set(
                        mapping.map_name(e, 'entrez', 'uniprot', self.target)
                    )
                )

            for r in target_r:

                target_u.update(
                    set(
                        mapping.map_name(e, 'refseqp', 'uniprot', self.target)
                    )
                )

            target_u = (
                itertools.chain(
                    *map(
                        lambda tu:
                            mapping.map_name(
                                tu, 'uniprot', 'uniprot', self.target
                            ),
                        target_u
                    )
                )
            )

            self.orthologs[source][u].update(
                {
                    Ortholog(u, 'HomoloGene')
                    for u in target_u
                }
            )


    def load_ensembl(self, source):

        source = self.get_source(source)
        target_organism = taxonomy.ensure_ensembl_name(self.target)
        source_organism = taxonomy.ensure_ensembl_name(source)

        _log(
            'Loading homology data from Ensembl '
            'between organisms `%u` and `%u`' % (source, self.target)
        )

        if not target_organism or not source_organism:

            self.ensembl[source] = {}
            _log(
                'No Ensembl homology data available between '
                'organisms `%s` and `%s`.' % (source, self.target)
            )
            return

        if source not in self.orthologs:

            self.orthologs[source] = collections.defaultdict(set)

        attr_target_ensp = '%s_homolog_ensembl_peptide' % target_organism
        attr_conf = '%s_homolog_orthology_confidence' % target_organism
        attr_type = '%s_homolog_orthology_type' % target_organism

        ensembl_data = biomart.biomart_homology(
            source_organism = source,
            target_organism = self.target,
        )

        for r in ensembl_data:

            source_uniprots = mapping.map_name(
                r.ensembl_peptide_id,
                'ensp',
                'uniprot',
                source,
            )

            target_uniprots = mapping.map_name(
                getattr(r, attr_target_ensp),
                'ensp',
                'uniprot',
                self.target,
            )

            for u in source_uniprots:

                self.orthologs[source][u].update(
                    {
                        Ortholog(
                            u,
                            'Ensembl',
                            ensembl_hc = getattr(r, attr_conf) == '1',
                            ensembl_type = getattr(r, attr_type),
                        )
                        for u in target_uniprots
                    }
                )


class PtmHomology(ProteinHomology, SequenceContainer):


    def __init__(
        self,
        target,
        source = None,
        only_swissprot = True,
        strict = True
    ):

        if not hasattr(self, '_logger'):

            session.Logger.__init__(self, name = 'homology')

        ProteinHomology.__init__(
            self,
            target = target,
            source = source,
            only_swissprot = only_swissprot,
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
    """
    Creates an instance of the homology manager. Stores it in the module
    namespace.
    """

    globals()['manager'] = HomologyManager()


def get_manager():
    """
    Returns the homology manager, an object which loads and unloads the
    homology lookup tables as necessary, and provides the interface for
    querying the homology data. Normally an instance of the manager
    belongs to the module, and if it does not exist yet, will be created
    automatically.
    """

    if 'manager' not in globals():

        init()

    return globals()['manager']


def translate(
        source_id,
        target,
        source = 9606,
        only_swissprot = True,
        homologene = None,
        ensembl = None,
        ensembl_hc = None,
        ensembl_types = None,
    ):
    """
    Homology translation. For a UniProt ID, finds the corresponding
    homologous (orthologous) genes in another organism.

    Args:
        source_id (str,list): UniProt ID of one or more protein in the
                source organism.
            target (int,str): The target organism.
            source (int,str): The source organism.
        only_swissprot (bool): Use only SwissProt IDs. For human and some
            popular model organisms this is advisible, as almost all proteins
            have reviewed record in UniProt.
        homologene (bool): Use NCBI HomoloGene data for ortholog lookup.
        ensembl (bool): Use Ensembl data for ortholog lookup.
        ensembl_hc (bool): Use only high confidence orthology relations
            from Ensembl. By default it is True. You can also set it
            by the `ensembl_hc` attribute.
        ensembl_types (list): The Ensembl orthology relationship types
            to use. Possible values are `one2one`, `one2many` and
            `many2many`. By default only `one2one` is used. You can
            also set this parameter by the `ensembl_types` attribute.

    Returns:
        Set of UniProt IDs of orthologous gene products in the target
        organism.
    """

    manager = get_manager()

    return manager.translate(
        source_id  = source_id,
        target = target,
        source = source,
        only_swissprot = only_swissprot,
        homologene = homologene,
        ensembl = ensembl,
        ensembl_hc = ensembl_hc,
        ensembl_types = ensembl_types,
    )
