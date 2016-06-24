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
import pypath.intera as intera
import pypath.uniprot_input as uniprot_input

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

class AttributeHandler(object):
    
    def __init__(self):
        pass
    
    def add_source(self, source):
        self.sources.add(source)
        if source not in self.attrs:
            self.attrs[source] = {}
    
    def merge_attrs(self, attrs):
        self.attrs = common.merge_dicts(self.attrs, attrs)
    
    def update_attr(self, attr):
        self.attrs = common.dict_set_path(self.attrs, attr)
    
    def __iadd__(self, other):
        '''
        Members or ids of entities should never change, as
        these are their unique, hashable and comparable
        attributes.
        
        __iadd__ operator is for merging entities with
        identical members or ids.
        '''
        self.sources = self.sources | other.sources
        self.merge_attrs(other.attrs)
        return self

class MolecularEntity(AttributeHandler):
    
    def __init__(self, identifier, id_type, sources = [], attrs = None):
        super(MolecularEntity, self).__init__()
        self.id = identifier
        self.id_type = id_type
        self.sources = set([])
        self.attrs = {}
        for source in sources:
            self.add_source(source)
        if attrs is not None:
            self.merge_attrs(attrs)
    
    def __str__(self):
        return '%s (%s)' % (self.id, self.id_type)
    
    def __hash__(self):
        return hash(self.__str__())
    
    def __eq__(self, other):
        return self.__str__() == other.__str__()
    
    def __gt__(self, other):
        return self.__str__() > other.__str__()
    
    def __lt__(self, other):
        return self.__str__() < other.__str__()
    
    def __repr__(self):
        return '%s (%s)' % (self.id, self.id_type)

class Protein(MolecularEntity):
    
    def __init__(self, protein_id, id_type = 'uniprot',
                 sources = [], attrs = None):
        super(Protein, self).__init__(protein_id, id_type,
                                      sources = sources, attrs = attrs)

class EntitySet(AttributeHandler):
    
    def __init__(self, members, sources = []):
        super(EntitySet, self).__init__()
        self.members = sorted(common.uniqList(members))
        self.sources = set([])
        self.attrs = {}
        for source in sources:
            self.add_source(source)
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
        return '%s: %s' % (self.__class__.__name__, self.__str__())
    
    def __hash__(self):
        return hash(self.__str__())
    
    def __eq__(self, other):
        return self.__str__() == other.__str__()
    
    def __gt__(self, other):
        return self.__str__() > other.__str__()
    
    def __lt__(self, other):
        return self.__str__() < other.__str__()
    
    def __iter__(self):
        for m in self.members:
            yield m
    
    def add_source(self, source):
        self.sources.add(source)
        if source not in self.attrs:
            self.attrs[source] = {}

class Complex(EntitySet):
    
    def __init__(self, members, source):
        sources = set([source]) if type(source) is not set else source
        super(Complex, self).__init__(members, sources)
        self.type = 'complex'
    
    def get_stoichiometries(self, source, cid, with_pids = False):
        if source in self.sources and cid in self.attrs[source]:
            return \
                list(
                    map(
                        lambda memb:
                            (memb[0], memb['stoi'], memb['pid']) \
                                if with_pids else (memb[0], memb['stoi']),
                        iteritems(self.attrs[source][cid])
                    )
                )

class ProteinFamily(EntitySet):
    
    def __init__(self, members, source):
        sources = set([source]) if type(source) is not set else source
        super(ProteinFamily, self).__init__(members, sources)
        self.type = 'pfamily'

class RePath(object):
    
    def __init__(self, mapper = None, ncbi_tax_id = 9606,
                default_id_types = {}, modifications = True,
                seq = None):
        
        self.ncbi_tax_id = ncbi_tax_id
        self.modifications = modifications
        self.parsers = {}
        self.sources = set([])
        self.seq = seq
        self.mapper = mapping.Mapper(ncbi_tax_id) if mapper is None else mapper
        
        self.species = {}
        self.proteins = {}
        self.pfamilies = {}
        self.complexes = {}
        self.reactions = {}
        self.mods = {}
        self.frags = {}
        
        self.rproteins = {}
        self.rpfamilies = {}
        self.rcomplexes = {}
        self.rreactions = {}
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
        self.merge_pfamilies()
        if self.modifications:
            self.merge_modifications()
        self.merge_complexes()
        self.remove_defaults()
    
    def remove_defaults(self):
        self.parser = None
        self.source = None
    
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
            target_ids = []
            id_attrs = {}
            for id_type, _id in ids:
                if id_type in self.id_types:
                    std_id_type = self.id_types[id_type]
                else:
                    std_id_type = id_type
                id_a = self.id_processor(_id)
                id_attrs[id_a['id']] = id_a
                target_ids.extend(
                    self.mapper.map_name(id_a['id'], std_id_type,
                                        self.default_id_types['protein'])
                )
                target_ids = common.uniqList(target_ids)
                if len(target_ids) > 1:
                    sys.stdout.write('\t:: Ambiguous ID translation: from %s to %s\n' % (ids, target_ids))
                elif len(target_ids) == 0:
                    target_ids = None
                else:
                    target_ids = target_ids[0]
            return target_ids, id_attrs
        
        self.rproteins[self.source] = {}
        for pid, p in iteritems(self.parser.proteins):
            ids = get_protein_ids(p['protein'])
            target_id, id_attrs = map_protein_ids(ids)
            if target_id is None:
                continue
            attrs = {
                self.source: {
                    'prefs': set([p['protein']]),
                    'pids': {
                        pid: {}
                    },
                    'originals': set([])
                }
            }
            for original_id, id_a in iteritems(id_attrs):
                attrs[self.source]['originals'].add(original_id)
                for k, v in iteritems(id_a):
                    if k != 'id':
                        attrs[self.source]['pids'][pid][k] = set([])
                        attrs[self.source]['pids'][pid][k].add(v)
            protein = Protein(target_id,
                                sources = set([self.source]), attrs = attrs)
            if target_id in self.proteins:
                self.proteins[target_id] += protein
            else:
                self.proteins[target_id] = protein
            self.rproteins[self.source][pid] = target_id
    
    def preprocess_seqmodvoc(self):
        kws = common.mod_keywords[self.source]
        aas = common.aanames
        self.seqmod_dict = {}
        for modkey, modname in iteritems(self.parser.seqmodvocs):
            for mod_std_name, kwlist in kws:
                if all(map(lambda kw: kw in modname, kwlist)):
                    this_aa = None
                    for aan, aa in iteritems(aas):
                        if aan in modname:
                            this_aa = aa
                            break
                    self.seqmod_dict[modkey] = (mod_std_name, this_aa)
    
    def merge_modifications(self):
        
        self.load_sequences()
        self.preprocess_seqmodvoc()
        
        def get_protein(pid):
            proteins = []
            if pid in self.rproteins[self.source]:
                _id = self.rproteins[self.source][pid]
                if 'isoform' in \
                    self.proteins[_id].attrs[self.source]['pids'][pid]:
                    for isof in self.proteins[_id].attrs[self.source]\
                        ['pids'][pid]['isoform']:
                        proteins.append((_id, isof))
                else:
                    proteins.append((_id, None))
            return proteins
        
        def get_seqsite(seqsite):
            if seqsite in self.parser.seqsites:
                return int(self.parser.seqsites[seqsite])
        
        def get_residue(protein, isof, resnum, resname):
            if protein in self.seq:
                if isof is not None and isof in self.seq[protein].isof:
                    sresname = self.seq[protein].get(resnum, isoform = isof)
                    if sresname == resname or resname is None:
                        return sresname, isof
                for isof in self.seq[protein].isoforms():
                    sresname = self.seq[protein].get(resnum, isoform = isof)
                    if sresname == resname or resname is None:
                        return sresname, isof
            return resname, isof
        
        for pid, p in iteritems(self.parser.proteins):
            proteins = get_protein(pid)
            
            for modfea in p['modfeatures']:
                
                if modfea in self.parser.fragfeas:
                    seqint = self.parser.fragfeas[modfea]
                    if seqint in self.parser.seqints:
                        start = get_seqsite(self.parser.seqints[seqint][0])
                        end = get_seqsite(self.parser.seqints[seqint][1])
                    if start is not None and end is not None:
                        for protein, isof in proteins:
                            if protein in self.seq:
                                if self.seq[protein].has_isoform(isof):
                                    frag = (protein, isof, start, end)
                                    if frag in self.frags:
                                        self.frags[frag].add_source(self.source)
                                    else:
                                        instance = \
                                            self.seq[protein].get(start, end, isof)
                                        mot = intera.Motif(protein, start, end,
                                                        isoform = isof,
                                                        instance = instance,
                                                        source = self.source)
                                        self.frags[frag] = mot
                                    self.proteins[protein].update_attr(
                                        [self.source, 'pids', pid, 'frags', frag]
                                    )
                
                if modfea in self.parser.modfeas:
                    resnum = get_seqsite(self.parser.modfeas[modfea][0])
                    seqmodvoc = self.parser.modfeas[modfea][1]
                    if seqmodvoc in self.seqmod_dict:
                        typ, resname = self.seqmod_dict[seqmodvoc]
                        for protein, isof in proteins:
                            if protein in self.seq:
                                resname, isof = \
                                    get_residue(protein, isof, resnum, resname)
                                mod = (protein, isof, resnum, resname)
                                if mod in self.mods:
                                    self.mods[mod].add_source(self.source)
                                else:
                                    res = intera.Residue(resnum, resname,
                                                            protein, isoform = isof)
                                    start, end, instance = self.seq[protein].get_region(
                                                            resnum, isoform = isof)
                                    mot = intera.Motif(protein, start, end,
                                                        isoform = isof,
                                                        instance = instance)
                                    ptm = intera.Ptm(protein, motif = mot,
                                                        residue = res,
                                                        source = self.source,
                                                        isoform = isof,
                                                        typ = typ)
                                    self.mods[mod] = ptm
                                self.proteins[protein].update_attr(
                                    [self.source, 'pids', pid, 'mods', mod]
                                )
    
    def merge_pfamilies(self):
        self.rpfamilies[self.source] = {}
        this_round = set(list(self.parser.pfamilies.keys()))
        next_round = []
        prev_round = -1
        while len(this_round) - prev_round != 0:
            prev_round = len(this_round)
            for pfid in this_round:
                pids = self.parser.pfamilies[pfid]
                subpf_unproc = \
                    any(map(lambda pid: pid in self.parser.pfamilies, pids))
                if subpf_unproc:
                    next_round.append(pfid)
                    continue
                proteins = \
                    list(
                        map(
                            lambda pid:
                                (self.rproteins[self.source][pid], pid),
                            filter(
                                lambda pid:
                                    pid in self.rproteins[self.source],
                                pids
                            )
                        )
                    )
                subpfs = \
                    list(
                        map(
                            lambda pid:
                                (self.rpfamilies[self.source][pid], pid),
                            filter(
                                lambda pid:
                                    pid in self.rpfamilies[self.source],
                                pids
                            )
                        )
                    )
                for spf, spfid in subpfs:
                    spfmembs = \
                        list(
                            map(
                                lambda p:
                                    (p[0], p[1]['pid']),
                                iteritems(self.pfamilies[spf].attrs[self.source][spfid])
                            )
                        )
                    proteins.extend(spfmembs)
                
                members = sorted(common.uniqList(map(lambda p: p[0], proteins)))
                pf = ProteinFamily(members, source = self.source)
                members = tuple(members)
                pf.attrs[self.source][pfid] = {}
                for protein, pid in proteins:
                    pf.attrs[self.source][pfid][protein] = {}
                    pf.attrs[self.source][pfid][protein]['pid'] = pid
                if members not in self.pfamilies:
                    self.pfamilies[members] = pf
                else:
                    self.pfamilies[members] += pf
                self.rpfamilies[self.source][pfid] = members
                
            this_round = next_round
            next_round = []
    
    def merge_complexes(self):
        self.rcomplexes[self.source] = {}
        this_round = set(list(self.parser.complexes.keys()))
        next_round = []
        prev_round = -1
        while len(this_round) - prev_round != 0:
            prev_round = len(this_round)
            for cid in this_round:
                stois = self.parser.complexes[cid]
                pids = list(map(lambda stoi:
                                    self.parser.stoichiometries[stoi], stois))
                subc_unproc = \
                    any(
                        map(
                            lambda pid:
                                pid[0] in self.parser.complexes,
                            pids
                        )
                    )
                if subc_unproc:
                    next_round.append(cid)
                    continue
                proteins = \
                    list(
                        map(
                            lambda pid:
                                (self.rproteins[self.source][pid[0]],
                                    pid[1], pid[0]),
                            filter(
                                lambda pid:
                                    pid[0] in self.rproteins[self.source],
                                pids
                            )
                        )
                    )
                subcplexs = \
                    list(
                        map(
                            lambda pid:
                                (self.rcomplexes[self.source][pid[0]],
                                    pid[1], pid[0]),
                            filter(
                                lambda pid:
                                    pid[0] in self.rcomplexes[self.source],
                                pids
                            )
                        )
                    )
                for sc, scstoi, scid in subcplexs:
                    scmembs = sc.get_stoichiometries(self.source,
                                                     scid, with_pids = True)
                    scmembs = list(map(lambda p:
                                           (p[0], p[1] * scstoi, p[2]), scmembs))
                    proteins.extend(scmembs)
                
                members = sorted(common.uniqList(map(lambda p: p[0], proteins)))
                cplex = Complex(members, source = self.source)
                members = tuple(members)
                cplex.attrs[self.source][cid] = {}
                for protein, stoi, pid in proteins:
                    cplex.attrs[self.source][cid][protein] = {}
                    cplex.attrs[self.source][cid][protein]['pid'] = pid
                    cplex.attrs[self.source][cid][protein]['stoi'] = stoi
                if members not in self.complexes:
                    self.complexes[members] = cplex
                else:
                    self.complexes[members] += cplex
                self.rcomplexes[self.source][cid] = members
                
            this_round = next_round
            next_round = []
    
    def merge_reactions(self):
        def get_side(ids):
            members = []
            for _id in ids:
                for typ in ('proteins', 'pfamilies', 'complexes'):
                    r = getattr(self, 'r%s'%typ)[self.source]
                    if _id in r:
                        members.append(getattr(self, typ)[r[_id]])
            return ReactionSide(members, source = self.source)
        
        for rid, reac in iteritems(self.parser.reactions):
            left = get_side(reac['left'])
            right = get_side(reac['right'])
            
            reaction = Reaction(left, right, )
            
            self.reactions
    
    def iterate_reactions(self):
        pass
    
    def load_sequences(self):
        if self.seq is None:
            self.seq = uniprot_input.swissprot_seq(self.ncbi_tax_id,
                                                   isoforms = True)

class ReactionSide(AttributeHandler):
    
    def __init__(self, members, sources = []):
        super(ReactionSide, self).__init__()
        
        self.members = sorted(members)
        self.source = set([])
        self.attrs = {}
        for source in sources:
            self.add_source(source)
    
    def __hash__(self):
        return hash(self.__str__())
    
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        return 'ReactionSide: (%s)' % \
            ('; '.join(map(lambda m: m.__str__(), self.members)))
    
    def __eq__(self, other):
        return self.__str__() == other.__str__()

class Reaction(AttributeHandler):
    
    def __init__(self, left, right, source = []):
        super(Reaction, self).__init__()
        
        self.left = ReactionSide(left, source)
        self.right = ReactionSide(left, source)
    
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        return 'Reaction: LEFT(%s) --> RIGHT(%s)' % \
            (self.left.__str__(), self.right.__str__())
    
    def __hash__(self):
        return hash(self.__str__())
    
    def __eq__(self, other):
        return self.__str__() == other.__str__()
