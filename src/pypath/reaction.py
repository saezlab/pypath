#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import imp
import re
import gzip
import tarfile
import zipfile
import struct
from lxml import etree

import pypath.mapping as mapping
import pypath.common as common
import pypath.progress as progress
import pypath.curl as curl

class BioPaxReader(object):
    
    def __init__(self, biopax, source, cleanup_period = 800):
        self.biopax = biopax
        self.source = source
        
        self.cleanup_period = cleanup_period
    
        # string constants
        self.bppref = '{http://www.biopax.org/release/biopax-level3.owl#}'
        self.rdfpref = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
        self.rdfid = '%sID' % self.rdfpref
        self.rdfab = '%sabout' % self.rdfpref
        self.rdfres = '%sresource' % self.rdfpref
        self.bpprot = '%sProtein' % self.bppref
        self.bpcplx = '%sComplex' % self.bppref
        self.bpprre = '%sProteinReference' % self.bppref
        self.bpreac = '%sBiochemicalReaction' % self.bppref
        self.bpcata = '%sCatalysis' % self.bppref
        self.bpctrl = '%sControl' % self.bppref
        self.bpcoma = '%sComplexAssembly' % self.bppref
        self.bppstp = '%sPathwayStep' % self.bppref
        self.bpuxrf = '%sUnificationXref' % self.bppref
        self.bpstoi = '%sStoichiometry' % self.bppref
        self.bppubr = '%sPublicationXref' % self.bppref
        self.bppath = '%sPathway' % self.bppref
        self.bpfrfe = '%sFragmentFeature' % self.bppref
        self.bpseqi = '%sSequenceInterval' % self.bppref
        self.bpseqs = '%sSequenceSite' % self.bppref
        self.bpmodf = '%sModificationFeature' % self.bppref
        self.bpmodv = '%sSequenceModificationVocabulary' % self.bppref
        self.bpmphe = '%smemberPhysicalEntity' % self.bppref
        self.bperef = '%sentityReference' % self.bppref
        self.bpxref = '%sxref' % self.bppref
        self.bprelr = '%sRelationshipXref' % self.bppref
        self.bpcsto = '%scomponentStoichiometry' % self.bppref
        self.bpstoc = '%sstoichiometricCoefficient' % self.bppref
        self.bpphye = '%sphysicalEntity' % self.bppref
        self.bpcted = '%scontrolled' % self.bppref
        self.bpcter = '%scontroller' % self.bppref
        self.bpctyp = '%scontrolType' % self.bppref
        self.bpleft = '%sleft' % self.bppref
        self.bprgth = '%sright' % self.bppref
        self.bpsprc = '%sstepProcess' % self.bppref
        self.bpfeat = '%sfeature' % self.bppref
        self.bpfelo = '%sfeatureLocation' % self.bppref
        self.bpibeg = '%ssequenceIntervalBegin' % self.bppref
        self.bpiend = '%ssequenceIntervalEnd' % self.bppref
        self.bpseqp = '%ssequencePosition' % self.bppref
        self.bpmoty = '%smodificationType' % self.bppref
        self.bppcom = '%spathwayComponent' % self.bppref
        self.bpterm = '%sterm' % self.bppref
        self.bpdb = '%sdb' % self.bppref
        self.bpid = '%sid' % self.bppref
        self.upStr = 'UniProt'
        
        self.proteins = {}
        self.pfamilies = {}
        self.complexes = {}
        self.cvariations = {}
        self.prefs = {}
        self.ids = {}
        self.reactions = {}
        self.cassemblies = {}
        self.stoichiometries = {}
        self.catalyses = {}
        self.controls = {}
        self.pwsteps = {}
        self.pubrefs = {}
        self.fragfeas = {}
        self.seqints = {}
        self.seqsites = {}
        self.modfeas = {}
        self.seqmodvocs = {}
        self.pathways = {}
        
        self.methods = {
            self.bpprot: 'protein',
            self.bpcplx: 'cplex',
            self.bpprre: 'pref',
            self.bpuxrf: 'uxref',
            self.bprelr: 'uxref',
            self.bpstoi: 'stoichiometry',
            self.bpreac: 'reaction',
            self.bpcoma: 'cassembly',
            self.bpcata: 'catalysis',
            self.bpctrl: 'control',
            self.bppstp: 'pwstep',
            self.bppubr: 'pubref',
            self.bpfrfe: 'fragfea',
            self.bpseqi: 'seqint',
            self.bpseqs: 'seqsite',
            self.bpmodf: 'modfea',
            self.bpmodv: 'seqmodvoc',
            self.bppath: 'pathway'
        }
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def process(self):
        self.biopax_size()
        self.set_progress()
        self.init_etree()
        self.iterate()
        self.close_biopax()
    
    def biopax_size(self):
        self.bp_filesize = 0
        if type(self.biopax) is tarfile.ExFileObject:
            self.bp_filesize = biopax.size
        elif type(self.biopax) is gzip.GzipFile:
            f = open(self.biopax.name, 'rb')
            f.seek(-4, 2)
            self.bp_filesize = struct.unpack('<I', f.read())[0]
            f.close()
        elif hasattr(self.biopax, 'name') and os.path.exists(self.biopax.name):
            self.bp_filesize = os.path.getsize(self.biopax.name)
    
    def init_etree(self):
        self.bp = etree.iterparse(self.biopax, events = ('end',))
        self.used_elements = []
    
    def set_progress(self):
        self.prg = progress.Progress(self.bp_filesize,
            'Processing %s from BioPAX XML' % self.source, 1)
    
    def iterate(self):
        self.fpos = self.biopax.tell()
        try:
            for ev, elem in self.bp:
                # step the progressbar:
                new_fpos = self.biopax.tell()
                self.prg.step(new_fpos - self.fpos)
                self.fpos = new_fpos
                self.next_elem = elem
                self.next_id = self.next_elem.get(self.rdfid) \
                    if self.rdfid in elem.attrib \
                    else self.next_elem.get(self.rdfab)
                if self.next_elem.tag in self.methods:
                    method = getattr(self, self.methods[self.next_elem.tag])
                    method()
                self.used_elements.append(self.next_elem)
                self.cleanup_hook()
        except etree.XMLSyntaxError as e:
            self.prg.terminate(status = 'failed')
            sys.stdout.write('\n\t:: Syntax error in BioPAX:\n\t\t%s\n' % str(e))
            sys.stdout.flush()
        self.prg.terminate()
    
    def cleanup_hook(self):
        if len(self.used_elements) > self.cleanup_period:
            for _ in xrange(int(self.cleanup_period / 2)):
                e = self.used_elements.pop()
                e.clear()
    
    def close_biopax(self):
        del self.bp
        self.biopax.close()
    
    def protein(self):
        entref = self.next_elem.find(self.bperef)
        if entref is not None:
            protein = self.get_none(entref.get(self.rdfres))
            self.proteins[self.next_id] = {
                'protein': protein,
                'seqfeatures': self._bp_collect_resources(self.bpfeat),
                'modfeatures': self._bp_collect_resources(self.bpfeat)
            }
        else:
            self.protein_families[self.next_id] = \
                self._bp_collect_resources(self.bpmphe)
    
    def pref(self):
        self.prefs[self.next_id] = \
            self._bp_collect_resources(self.bpxref)
    
    def uxref(self):
        db = self.next_elem.find(self.bpdb)
        if db is not None:
            id_type = db.text.lower()
            i = self.next_elem.find(self.bpid)
            if i is not None:
                self.ids[self.next_id] = (id_type, i.text)
    
    def cplex(self):
        if self.next_elem.find(self.bpcsto) is not None:
            self.complexes[self.next_id] = \
                self._bp_collect_resources(self.bpcsto)
        else:
            self.cvariations[self.next_id] = \
                self._bp_collect_resources(self.bpmphe)
    
    def stoichiometry(self):
        self.stoichiometries[self.next_id] = (
            self.next_elem.find(self.bpphye).get(self.rdfres).replace('#', ''),
            int(float(self.next_elem.find(self.bpstoc).text))
        )
    
    def reaction(self):
        self.reactions[self.next_id] = {
            'refs': self._bp_collect_resources(self.bpxref),
            'left': self._bp_collect_resources(self.bpleft),
            'right': self._bp_collect_resources(self.bprgth)
        }
    
    def cassembly(self):
        self.cassemblies[self.next_id] = {
            'refs': self._bp_collect_resources(self.bpxref),
            'left': self._bp_collect_resources(self.bpleft),
            'right': self._bp_collect_resources(self.bprgth)
        }
    
    def catalysis(self):
        cter = self.next_elem.find(self.bpcter)
        cted = self.next_elem.find(self.bpcted)
        if cter is not None and cted is not None:
            typ = self.next_elem.find(self.bpctyp)
            self.catalyses[self.next_id] = {
                'controller': self.get_none(cter.get(self.rdfres)),
                'controlled': self.get_none(cted.get(self.rdfres)),
                'type': '' if typ is None else typ.text
            }
    
    def control(self):
        cter = self.next_elem.find(self.bpcter)
        cted = self.next_elem.find(self.bpcted)
        if cter is not None and cted is not None:
            typ = self.next_elem.find(self.bpctyp)
            self.controls[self.next_id] = {
                'refs': _bp_collect_resources(elem, bpxref),
                'controller': cter.get(rdfres).replace('#', ''),
                'controlled': cted.get(rdfres).replace('#', ''),
                'type': '' if typ is None else typ.text
            }
    
    def pwstep(self):
        self.pwsteps[_id] = self._bp_collect_resources(self.bppstp)
    
    def pubref(self):
        pmid = self.next_elem.find(self.bpid)
        if pmid is not None:
            self.pubrefs[self.next_id] = pmid.text
    
    def fragfea(self):
        self.fragfeas[self.next_id] = \
            self.next_elem.find(self.bpfelo).get(self.rdfres).replace('#', '')
    
    def seqint(self):
        beg = self.next_elem.find(self.bpibeg)
        end = self.next_elem.find(self.bpiend)
        self.seqints[self.next_id] = (
            beg.get(self.rdfres).replace('#', '') if beg is not None else None,
            end.get(self.rdfres).replace('#', '') if end is not None else None
        )
    
    def seqsite(self):
        eqp = self.next_elem.find(self.bpseqp)
        if seqp is not None:
            self.seqsites[self.next_id] = int(seqp.text)
    
    def modfea(self):
        felo = self.next_elem.find(self.bpfelo)
        moty = self.next_elem.find(self.bpmoty)
        if felo is not None and moty is not None:
            self.modfeas[self.next_id] = (
                self.next_elem.find(self.bpfelo).\
                    get(self.rdfres).replace('#', ''),
                self.next_elem.find(self.bpmoty).\
                    get(self.rdfres).replace('#', '')
            )
    
    def seqmodvoc(self):
        term = self.next_elem.find(self.bpterm)
        if term is not None:
            self.seqmodvocs[self.next_id] = term.text
    
    def pathway(self):
        try:
            self.pathways[self.next_id] = {
                'reactions': self._bp_collect_resources(self.bppcom),
                'pathways': self._bp_collect_resources(self.bppcom)
            }
        except TypeError:
            sys.stdout.write('Wrong type at element:\n')
            sys.stdout.write(etree.tostring(self.next_elem))
            sys.stdout.flush()
    
    def get_none(self, something):
        if something is not None:
            something.replace('#', '')
        return something
    
    def _bp_collect_resources(self, tag, restype = None):
        return \
        map(
            lambda e:
                e.get(self.rdfres).replace('#', ''),
            filter(
                lambda e:
                    self.rdfrs in e.attrib and (\
                        restype is None or \
                        e.get(self.rdfres).replace('#', '').startswith(restype)
                    ),
                self.next_elem.iterfind(tag)
            )
        )

class RePath(object):
    
    def __init__(self, mapper = None, ncbi_tax_id = 9606):
        self.ncbi_tax_id = ncbi_tax_id
        self.inputs = {}
        self.sources = set([])
        self.mapper = mapping.Mapper(ncbi_tax_id) if mapper is None else mapper
        self.species = {}
        self.reactions = {}
        self.references = {}
    
    def add_source(self, source, infile):
        self.sources.add(source)
        self.inputs[source] = infile
    
    def read_biopax(self, source):
        pass

class Reaction(object):
    
    def __init__(self):
        pass

class Species(object):
    
    def __init__(self):
        pass
