#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# cogent and sqlalchemy modules need to be installed:
# pip2 install cogent
# pip2 install sqlalchemy
from cogent.db.ensembl import HostAccount, Species, Genome
from pypath import mapping

Release = 78
account = HostAccount('ensembldb.ensembl.org', 'anonymous', '')

human = Genome(Species='human', Release=Release, account=account)

# UniProt, seq offset, residue, isoform
positions = [('P00533', 40, 'Q', 1), ('P60520', 30, 'P', 1)]

m = mapping.Mapper()
m.load_uniprot_mappings(['ensg'], bi=True)

positions_ens = []
for p in positions:
    ensgs = m.map_name(p[0], 'uniprot', 'ensg')
    for ensg in ensgs:
        genes = human.getGenesMatching(StableId=ensg)
        for gene in genes:
            positions_ens.append(
                tuple([ensg, gene.Location, gene.CanonicalTranscript.Exons] +
                      list(p)))

# another attempts with biopython --
# (it works, if you can map all proteins to RefSeq Gene IDs)

from Bio.SeqUtils.Mapper import CoordinateMapper
from Bio import SeqIO
from Bio import Entrez
from pypath import mapping, dataio

Entrez.email = 'denes@ebi.ac.uk'

refsg = dataio.curl('ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/'
                    'RefSeqGene/LRG_RefSeqGene')

prot2gen = {}
refsg = [x.split('\t') for x in refsg.split('\n') if len(x) > 0][1:]
for r in refsg:
    prot2gen[r[7].split('.')[0]] = r[3].split('.')[0]


def get_first_CDS(parser):
    for rec in parser:
        for feat in rec.features:
            if feat.type == "CDS" and len(feat.location.parts) > 1:
                return feat

# UniProt, seq offset, residue, isoform
positions = [('P00533', 40, 'Q', 1), ('P60520', 30, 'P', 1)]

m = mapping.Mapper()
m.load_uniprot_mappings(['refseqp'], bi=True)

positions_g = []
for p in positions:
    refseq = m.map_name(p[0], 'uniprot', 'refseqp')
    for rsp in refseq:
        rsg = prot2gen[rsp]
        record = Entrez.esearch(db='nucleotide', term=rsg)
        record = Entrez.read(record)
        if len(record['IdList']) > 0:
            gi = record['IdList'][0]
            handle = Entrez.efetch(
                db="nucleotide", id=gi, rettype="gb", retmode="text")
            seq = SeqIO.parse(handle, 'genbank')
            exons = get_first_CDS(seq)

exons
