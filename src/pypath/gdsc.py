#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This module (pypath.gdsc) is not available for public,
#  please use only within EBI.
#
#  Website: http://www.ebi.ac.uk/~denes
#

# external modules:
import os
import sys
import re
import numpy
import pandas

# from this package:
import intera
import seq
import progress
from common import wcl
from dataio import curl, RemoteFile

# constants
EBI_HOST = 'ebi-005.ebi.ac.uk'

class GDSC(object):
    
    def __init__(self, datadir = None):
        self.datadir = datadir if datadir is not None \
            else '/home/denes/Dokumentumok/pw/data'
    
    def read_mutations(self, sample_col = 'COSMIC_ID', 
            attributes = [], datadir = None, infile = None):
        '''
        The file processed by Luz
        Original data:
        /nfs/research2/saezrodriguez/jsr-gdsc/RAW/Genomic/
        All_variants_cell-lines_22102014.xlsx
        Original location: 
        /nfs/research2/saezrodriguez/jsr-gdsc/PROCESSED/Genomic/
        All_variants_cell-lines_22102014_indelMod_ANNOVfun.txt
        '''
        result = {}
        datadir = datadir if datadir is not None \
            else '/home/denes/Dokumentumok/pw/data'
        infile = infile if infile is not None \
            else 'All_variants_cell-lines_22102014_indelMod_ANNOVfun.txt'
        infile = os.path.join(self.datadir, infile)
        size = os.path.getsize(infile)
        remut = re.compile(r'([A-Z])([0-9]*)([A-Z])')
        prg = progress.Progress(size, 'Reading mutations', 7)
        with open(infile, 'r') as f:
            hdr = f.readline().split('\t')
            cols = {}
            for i, col in enumerate(hdr):
                cols[col] = i
            for l in f:
                prg.step(len(l))
                l = l.split('\t')
                if not l[5].startswith('syn'):
                    uniprot = l[11]
                    mutation = remut.findall(l[12])
                    if len(mutation) != 0:
                        mutation = mutation[0]
                        this_attrs = {}
                        for attr in attributes:
                            if attr in cols:
                                this_attrs[attr] = l[cols[attr]].strip()
                        smpl = 0 if sample_col not in cols else \
                            int(l[cols[sample_col]].strip())
                        ori = intera.Residue(mutation[1], mutation[0], uniprot)
                        mtd = intera.Residue(mutation[1], mutation[2], 
                            uniprot, mutated = True)
                        mut = intera.Mutation(ori, mtd, smpl, this_attrs)
                        if uniprot not in result:
                            result[uniprot] = []
                        result[uniprot].append(mut)
        prg.terminate()
        return result

    def read_transcriptomics(self, infile = 'normalized'):
        '''
        This fun reads the raw transcriptomics data
        original data:
        /nfs/research2/saezrodriguez/jsr-gdsc/RAW/Transcriptomic/01_RMAproc_basal_exp.csv
        returns a dict of dicts
        '''
        datadir = datadir if datadir is not None else '/home/denes/Dokumentumok/pw/data'
        files = {
            'raw': '01_RMAproc_basal_exp.csv',
            'normalized': '02_d_norm_basal_exp.csv',
            'tissue': '03_d_norm_tissue_centered_basal_exp.csv'
        }
        infile = infile if infile not in files else files[infile]
        infile = os.path.join(self.datadir, infile)
        result = {}
        cols = {}
        size = os.path.getsize(infile)
        sys.stdout.write('\t:: Reading transcriptomics data from\n'\
            '\t\t%s\n'%infile)
        sys.stdout.flush()
        prg = progress.Progress(size, 'Reading transcriptomics data', 7)
        with open(infile, 'r') as f:
            # cosmic ids as integers:
            try:
                hdr = [int(x) for x in f.readline().strip().split(',')[1:]]
            except:
                f.seek(0)
                print f.readline().split(',')
            for l in f:
                prg.step(len(l))
                l = l.split(',')
                result[l[0]] = dict(zip(hdr, [float(x) for x in l[1:]]))
        prg.terminate()
        return result

    def cell_lines(self, infile = None):
        self.cosmic = {}
        self.cellines = {}
        infile = 'cell_line_list.tab' if infile is None else infile
        infile = os.path.join(self.datadir, infile)
        with open(infile, 'r') as f:
            headers = [hdr.lower().strip().replace(' ', '_') \
                for hdr in f.readline().strip().split('\t')]
            for r in f:
                r = [x.strip() for x in r.split('\t')]
                self.cellines[r[0]] = {}
                for index, hdr in enumerate(headers):
                    self.cellines[r[0]][hdr] = r[index]
                if r[2].isdigit():
                    self.cosmic[r[2]] = r[0]
    
    def cell_line_lookup(self, field, value):
        result = []
        value = set(value) if type(value) is list else set([value])
        for celll, data in self.cellines.iteritems():
            if field in data and data[field] in value and \
                data['cosmic_id'] in self.cosmic:
                result.append(data['cosmic_id'])
        return result
    
    def remote_open(self, filename, user, passwd, host = None, port = 22):
        host = host if host is not None else EBI_HOST
        return RemoteFile(filename, user, host, passwd, port)
    
    def get_microarray(self, f, cosmix = None, uniprots = None, header = None, sep = '\t'):
        '''
        Reads the expression vector of one sample from processed microarray data.
        
        @f : generator or file like object
            File on the server can be opened with remote_open(), 
            local file with open(). `f` migh be any iterable object
            yielding one line per loop.
        @cosmix : list
            List of sample IDs (column header in file), here we use COSMIC IDs.
        @uniprots : list
            If None, expression of all genes will be returned, otherwise
            only of those in the list.
        @header : list
            If None, the first element of the iterable/file `f` will be used as header.
        @sep : str
            Field separator of file.
        '''
        if type(f) is file:
            rownames = []
            if header is None:
                l = f.readline()
            for l in f:
                l = l.strip().split(sep)
                rownames.append(l[0])
        elif type(f) is list:
            rownames = [l[0] for l in f]
        elif type(f) is RemoteFile:
            rownames = f.rowns()
        if uniprots is not None:
            uniprots = uniprots if type(uniprots) is set else set(uniprots)
            rownames = list(set(rownames) & uniprots)
        colnums = None
        tbl = f if type(f) is not RemoteFile else f.open()
        prg = progress.Progress(len(rownames), 'Reading transcriptomic data', 1)
        for l in tbl:
            prg.step()
            l = l if type(l) is list else l.rstrip().split(sep)
            if colnums is None and header is not None:
                colnums = dict([(cosmic, header.index(cosmic) + 1) \
                    for cosmic in cosmix if len(cosmic) > 0 and cosmic in header])
                # result = pandas.DataFrame(index = rownames, columns = header)
                result = numpy.array([[0.0] * len(header)] * len(rownames))
            if header is None:
                header = l[1:]
                if cosmix is None:
                    cosmix = header
                continue
            if uniprots is None or l[0] in uniprots:
                # result.loc[l[0]] = [float(l[colnums[cosmic]]) for cosmic in header]
                result[rownames.index(l[0])] = \
                    [float(l[colnums[cosmic]]) for cosmic in header]
        result = pandas.DataFrame(result, index = rownames, columns = header)
        prg.terminate()
        return result