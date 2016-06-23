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
    
    def __init__(self, biopax, source, cleanup_period = 800,
        file_from_archive = None):
        self.biopax = biopax
        self.source = source
        self.file_from_archive = file_from_archive
        self.cleanup_period = cleanup_period
        self.biopax_tmp_file = None
        self.cachedir = 'cache'
    
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
        self.open_biopax()
        self.biopax_size()
        self.extract()
        self.set_progress()
        self.init_etree()
        self.iterate()
        self.close_biopax()
    
    def open_biopax(self):
        if type(self.biopax) is curl.FileOpener:
            self.opener = biopax
        else:
            self.opener = curl.FileOpener(self.biopax)
        if type(self.opener.result) is dict:
            if self.file_from_archive is None or \
                self.file_from_archive not in self.opener.result:
                self.file_from_archive = \
                    sorted(list(self.opener.result.keys()))[0]
            self._biopax = self.opener.result[self.file_from_archive]
        else:
            self._biopax = self.opener.fileobj
    
    def biopax_size(self):
        self.bp_filesize = self.opener.sizes[self.file_from_archive] \
            if hasattr(self.opener, 'sizes') else self.opener.size
    
    def extract(self):
        if self.opener.type != 'plain':
            self.biopax_tmp_file = os.path.join(
                self.cachedir, 'biopax.processing.tmp.owl')
            prg = progress.Progress(self.bp_filesize,
                                    'Extracting %s from %s compressed file' % \
                                        (self.file_from_archive,
                                         self.opener.type),
                                    1000000)
            with open(self.biopax_tmp_file, 'wb') as tmpf:
                while True:
                    chunk = self._biopax.read(100000)
                    prg.step(len(chunk))
                    if not len(chunk):
                        break
                    tmpf.write(chunk)
            self._biopax = open(self.biopax_tmp_file, 'rb')
            prg.terminate()
    
    def init_etree(self):
        self.bp = etree.iterparse(self._biopax, events = ('start', 'end'))
        _, self.root = next(self.bp)
        self.used_elements = []
    
    def set_progress(self):
        self.prg = progress.Progress(self.bp_filesize,
            'Processing %s from BioPAX XML' % self.source, 33)
    
    def iterate(self):
        self.fpos = self._biopax.tell()
        try:
            for ev, elem in self.bp:
                # step the progressbar:
                new_fpos = self._biopax.tell()
                self.prg.step(new_fpos - self.fpos)
                self.fpos = new_fpos
                self.next_elem = elem
                self.next_event = ev
                self.next_id = self.next_elem.get(self.rdfid) \
                    if self.rdfid in elem.attrib \
                    else self.next_elem.get(self.rdfab)
                if ev == 'end' and self.next_elem.tag in self.methods:
                    method = getattr(self, self.methods[self.next_elem.tag])
                    method()
                    self.root.clear()
                #self.used_elements.append(self.next_elem)
                #self.cleanup_hook()
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
        self._biopax.close()
    
    def protein(self):
        entref = self.next_elem.find(self.bperef)
        if entref is not None:
            protein = self.get_none(entref.get(self.rdfres)).replace('#', '')
            self.proteins[self.next_id] = {
                'protein': protein,
                'seqfeatures': self._bp_collect_resources(self.bpfeat),
                'modfeatures': self._bp_collect_resources(self.bpfeat)
            }
        else:
            self.pfamilies[self.next_id] = \
                self._bp_collect_resources(self.bpmphe)
    
    def pref(self):
        self.prefs[self.next_id] = \
            self._bp_collect_resources(self.bpxref)
    
    def uxref(self):
        db = self.next_elem.find(self.bpdb)
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
        snum = self.next_elem.find(self.bpstoc).text
        self.stoichiometries[self.next_id] = (
            self.next_elem.find(self.bpphye).get(self.rdfres).replace('#', ''),
            int(float(snum))
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
                'refs': self._bp_collect_resources(self.bpxref),
                'controller': cter.get(self.rdfres).replace('#', ''),
                'controlled': cted.get(self.rdfres).replace('#', ''),
                'type': '' if typ is None else typ.text
            }
    
    def pwstep(self):
        self.pwsteps[self.next_id] = self._bp_collect_resources(self.bppstp)
    
    def pubref(self):
        pmid = self.next_elem.find(self.bpid)
        if pmid is not None:
            self.pubrefs[self.next_id] = pmid.text
    
    def fragfea(self):
        felo = self.next_elem.find(self.bpfelo).get(self.rdfres)
        self.fragfeas[self.next_id] = felo.replace('#', '')
    
    def seqint(self):
        beg = self.next_elem.find(self.bpibeg)
        end = self.next_elem.find(self.bpiend)
        self.seqints[self.next_id] = (
            beg.get(self.rdfres).replace('#', '') if beg is not None else None,
            end.get(self.rdfres).replace('#', '') if end is not None else None
        )
    
    def seqsite(self):
        seqp = self.next_elem.find(self.bpseqp)
        if seqp is not None and seqp.text is not None:
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
        list(
            map(
                lambda e:
                    e.get(self.rdfres).replace('#', ''),
                filter(
                    lambda e:
                        self.rdfres in e.attrib and (\
                            restype is None or \
                            e.get(self.rdfres).replace('#', '').startswith(restype)
                        ),
                    self.next_elem.iterfind(tag)
                )
            )
        )

class MolecularEntity(object):
    
    def __init__(self, identifier, id_type):
        self.id = identifier
        self.id_type = id_type
        self.sources = set([])
        self.attrs = {}
    
    def __str__(self):
        return '%s (%s)' % (self.id, self.id_type)
    
    def __hash__(self):
        return hash(self.__str__())
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __repr__(self):
        return '%s (%s)' % (self.id, self.id_type)
    
    def add_source(self, source):
        self.sources.add(source)
        if source not in self.attrs:
            self.attrs[source] = {}
    
    def __iadd__(self, other):
        self.sources = self.sources | other.sources
        return self

class Protein(MolecularEntity):
    
    def __init__(self, uniprot, sources = []):
        super(Protein, self).__init__(uniprot, 'uniprot')
        for source in sources:
            self.add_source(source)

class EntitySet(object):
    
    def __init__(self, members, sources = set([])):
        self.members = sorted(common.uniqList(members))
        self.sources = sources
        self.originals = {}
        self.type = None
    
    def __str__(self):
        return ';'.join(
                    sorted(
                        map(
                            lambda x:
                                str(x),
                            list(self.members)
                        )
                    )
                )
    
    def __repr__(self):
        return '%s: %s' (self.__class__.__name__, self.__str__())
    
    def __hash__(self):
        return hash(self.__str__())
    
    def __iter__(self):
        for m in self.members:
            yield m
    
    def __iadd__(self, other):
        self.members.extend(other.members)
        self.sources = self.sources | other.sources
        return globals()[self.__class__.__name](self.members)
    
    def add_source(self, source):
        self.sources.add(source)

class Complex(EntitySet):
    
    def __init__(self, members, source):
        sources = set([source]) if type(source) is not set else source
        super(Complex, self).__init__(members, sources)
        self.type = 'complex'

class ProteinFamily(EntitySet):
    
    def __init__(self, members, source):
        sources = set([source]) if type(source) is not set else source
        super(ProteinFamily, self).__init__(members, sources)
        self.type = 'pfamily'

class RePath(object):
    
    def __init__(self, mapper = None, ncbi_tax_id = 9606,
                default_id_types = {}):
        self.ncbi_tax_id = ncbi_tax_id
        self.parsers = {}
        self.sources = set([])
        self.mapper = mapping.Mapper(ncbi_tax_id) if mapper is None else mapper
        self.species = {}
        self.proteins = {}
        self.reactions = {}
        self.references = {}
        self.id_types = {}
        self.default_id_types = {
            'protein': 'uniprot'
        }
        self.default_id_types.update(default_id_types)
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def add_dataset(self, parser, id_types = {},
                    process_id = lambda x: {'id': x}):
        self.id_types.update(id_types)
        self.source = parser.source
        self.parser = parser
        self.id_processor = process_id
        self.merge()
    
    def merge(self):
        self.sources.add(self.source)
        self.parsers[self.source] = self.parser
        self.merge_proteins()
    
    def merge_proteins(self):
        def get_protein_ids(pref):
            
            pids = []
            if pref in self.parser.prefs:
                uxrefs = self.parser.prefs[pref]
                pids = \
                    common.uniqList(
                        map(
                            lambda uxref:
                                self.parser.ids[uxref],
                            filter(
                                lambda uxref:
                                    uxref in self.parser.ids,
                                uxrefs
                            )
                        )
                    )
            return pids
        
        def map_protein_ids(ids):
            print('mapping ids: %s' % ids)
            target_ids = []
            id_attrs = {}
            for id_type, _id in ids:
                if id_type in self.id_types:
                    std_id_type = self.id_types[id_type]
                else:
                    std_id_type = id_type
                id_a = self.id_processor(_id)
                id_attrs[id_a['id']] = id_a
                print('mapping from: %s, to: %s' % (std_id_type, self.default_id_types['protein']))
                target_ids.extend(
                    self.mapper.map_name(id_a['id'], std_id_type,
                                        self.default_id_types['protein'])
                )
            return target_ids, id_attrs
        
        for pid, p in iteritems(self.parser.proteins):
            ids = get_protein_ids(p['protein'])
            target_ids, id_attrs = map_protein_ids(ids)
            print('got ids: %s' % target_ids)
            for target_id in target_ids:
                print('target id: %s' % target_id)
                protein = Protein(target_id, sources = set([self.source]))
                if target_id in self.proteins:
                    self.proteins[target_id] += protein
                else:
                    self.proteins[target_id] = protein

class Reaction(object):
    
    def __init__(self):
        pass

class Species(object):
    
    def __init__(self):
        pass

# Stoichiometry15
# UnificationXref330
