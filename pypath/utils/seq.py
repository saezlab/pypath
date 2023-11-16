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

"""
Utilities for working with protein sequences.
"""

from __future__ import annotations

from future.utils import iteritems

import os
import sys
import re
import collections

import pypath.share.common as common
import pypath.inputs.pfam as pfam_input
import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.uniprot_db as uniprot_db
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.taxonomy as taxonomy


def swissprot_seq(
        organism: str | int = 9606,
        reviewed: bool = True,
        isoforms: bool = False,
    ) -> dict[str, Seq]:
    """
    All UniProt sequences for an organism.

    Loads all sequences for an organism, optionally
    for all isoforms, by default only first isoform.

    Args:
        organism:
            Name or NCBI Taxonomy ID of the organism.
        reviewed:
            Load only reviewed (SwissProt) sequences.
        isoforms:
            Load all isoforms, not only the first one.

    Returns:
        A dict with UniProt IDs as keys and `Seq` objects
        as values.
    """

    data = uniprot_input.uniprot_data(
        fields = 'sequence',
        organism = organism,
        reviewed = reviewed or None,
    )
    result = {u: Seq(u, s) for u, s in data.items()}

    if isoforms:

        data = get_isoforms(organism = organism)

        for unip, isoforms in iteritems(data):

            for isof, seq in iteritems(isoforms):

                if unip in result:

                    result[unip].add_seq(seq, isof)

    return result


def get_isoforms(organism = 9606):
    """
    Loads UniProt sequences for all isoforms.
    """

    if organism in taxonomy.phosphoelm_taxids:
        organism = taxonomy.phosphoelm_taxids[organism]

    reorg = re.compile(r'OS=([A-Z][a-z]+\s[a-z]+)')
    result = {}
    url = urls.urls['unip_iso']['url']
    c = curl.Curl(url, silent=False)
    data = c.result
    data = read_fasta(data)

    for header, seq in iteritems(data):

        org = reorg.findall(header)

        if len(org) > 0 and org[0] == organism:

            prot = header.split('|')[1].split('-')
            unip = prot[0]
            isof = int(prot[1])

            if unip not in result:
                result[unip] = {}

            result[unip][isof] = seq

    return result


def read_fasta(fasta):
    """
    Parses a fasta file.
    Returns dict with headers as keys and sequences as values.
    """

    result = {}
    fasta = re.split(r'\n>', fasta)

    for section in fasta:

        section = section.strip().split('\n')
        label = section.pop(0)
        seq = ''.join(section)
        result[label] = seq

    return result


class Resource(object):


    def __init__(self, loader, name = None):
        """
        Represents a resource of sequence features,
        e.g. domains, motifs or post-translational
        modification sites.
        """

        self.loader = loader
        self.db = {}
        self.name = name


    def load(self, ncbi_tax_id = 9606):
        """
        Loads the data from the resource for a given organism.
        """

        if ncbi_tax_id not in self.db:

            self.db[ncbi_tax_id] = list(
                self.processor(self.loader(ncbi_tax_id = ncbi_tax_id))
            )


    def unload(self, ncbi_tax_id = None):
        """
        Removes data in order to free up memory.
        """

        if ncbi_tax_id in self.db:

            del self.db[ncbi_tax_id]

        elif ncbi_tax_id is None:

            self.db = {}

    def processor(self, raw):
        """
        Preprocesses the features loaded from a resource.
        """

        for feature in raw:

            yield feature

    def iterprotein(self, uniprot, ncbi_tax_id = 9606):
        """
        Iterates over the features of one protein.
        """

        self.load(ncbi_tax_id)

        if uniprot in self.db[ncbi_tax_id]:

            for feature in self.db[ncbi_tax_id][uniprot]:

                yield feature

    def iterdb(self, ncbi_tax_id = 9606):
        """
        Iterates over all proteins and features of one organism.
        """

        self.load(ncbi_tax_id)

        for uniprot in self.db[ncbi_tax_id]:

            for feature in self.iterprotein(uniprot, ncbi_tax_id):

                yield (uniprot, feature)


class Pfam(Resource):

    def __init__(self):
        """
        Provides Pfam domains as sequence features.
        """

        Resource.__init__(self, loader = None, name = 'Pfam')

        def loader(ncbi_tax_id = 9606):

            all_up = uniprot_db.all_uniprots(organism = ncbi_tax_id)

            return (
                pfam_input.pfam_regions(
                    uniprots = all_up,
                    value = 'uniprot',
                    keepfile = True,
                )
            )

        self.loader = loader


class Seq(object):


    def __init__(self, protein, sequence, isoform=1):
        """
        This class is to look up or match
        residues and regions in sequences of
        proteins, by default in the canonical
        sequence, and optionally in other isoforms.
        """

        self.isof = {}
        self.protein = protein
        self.canonical = isoform
        self.add_seq(sequence, isoform)


    def add_seq(self, sequence, isoform):
        self.isof[isoform] = sequence


    def match(self, pattern, start, end=None, isoform=None):

        instance = self.get(start, end, isoform)
        pattern = pattern.upper()

        return instance == pattern


    def get(self, start, end=None, isoform=None):
        isoform = isoform if isoform is not None else self.canonical
        end = end if end is not None else start
        return None if len(self.isof[isoform]) < max(start, end) or min(start, end) < 1 \
            else self.isof[isoform][start - 1:end]


    def isoforms(self):
        return list(self.isof.keys())


    def has_isoform(self, isoform):
        return isoform in self.isof


    def get_region(self,
                   residue=None,
                   start=None,
                   end=None,
                   flanking=7,
                   isoform=None):
        isoform = self.canonical if isoform is None else isoform
        if residue is None and start is None and end is None:
            return (1, len(self.isof[isoform]), self.isof[isoform])
        if residue is not None and residue > len(self.isof[isoform]):
            return (None, None, None)
        start = start if start is not None else residue - flanking
        end = end if end is not None else residue + flanking
        start = max(start, 1)
        end = min(end, len(self.isof[isoform]))
        return (start, end, self.isof[isoform][start - 1:end])


    def findall(self, fragment):
        """
        Looks up a sequence fragment in the sequences.
        Yields tuples of isoform and offset.
        """

        SeqLookup = collections.namedtuple(
            'SeqLookup',
            ['isoform', 'offset'],
        )

        for iso, se in iteritems(self.isof):

            offset = 0

            while True:

                offset = se.find(fragment, offset)

                if offset == -1:

                    break

                yield SeqLookup(iso, offset)

                offset += 1

    def get_biopython(self, isoform = 1):

        isoform = int(isoform)

        if isoform not in self.isof:

            raise ValueError('No isoform %u available for protein `%s`.' % (
                isoform, self.protein))

        try:
            import Bio.Seq
            import Bio.SeqRecord

            srec = Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(self.isof[isoform],
                            Bio.Alphabet.ProteinAlphabet()),
                id = self.protein
            )

            srec.annotations['isoform'] = isoform

            return srec

        except ImportError:
            sys.stdout.write('\t:: Module `Bio` (biopython)'\
                'could not be imported.\n')
            sys.stdout.flush()


    def export_fasta(self, fname = None, sequences = None):

        sequences = sequences or [self]
        fname = fname or '%s.fasta' % self.protein

        try:
            import Bio.SeqIO

            Bio.SeqIO.write([s.get_biopython() for s in sequences],
                            fname, 'fasta')

        except ImportError:
            sys.stdout.write('\t:: Module `Bio` (biopython)'\
                'could not be imported.\n')
            sys.stdout.flush()


    def multiple_alignment(self, sequences, outfile = None,
                           method = 'ClustalW', param = {}):

        try:
            import Bio.Align.Applications

        except ImportError:

            sys.stdout.write('\t:: Module `Bio` (biopython)'\
                'could not be imported.\n')
            sys.stdout.flush()
            return

        session = common.random_string()

        method = method.capitalize()

        infile  = os.path.join('cache', '_align.%s.tmp.fasta' % session)
        keep_outfile = outfile is not None
        outfile = outfile or os.path.join('cache',
                                          '_align.%s.tmp.aln' % session)
        self.export_fasta(infile, sequences)

        if method == 'Muscle':
            if 'clw' not in param:
                param['clw'] = True
            param['input'] = infile
            param['out']   = outfile

        if method == 'Clustalomega':
            method = 'ClustalOmega'
            param['infile'] = infile
            param['outfile'] = outfile

        if method == 'Clustalw':
            param['cmd'] = 'clustalw2'
            param['infile'] = infile
            param['outfile'] = outfile

        app = getattr(Bio.Align.Applications, '%sCommandline' % method)

        cmd = app(**param)
        cmd()

        os.remove(infile)

        import Bio.AlignIO

        if method.lower() == 'clustalw' or param['clw']:
            aln = Bio.AlignIO.read(outfile, 'clustal')
        else:
            aln = None

        if not keep_outfile:
            os.remove(outfile)

        return aln
