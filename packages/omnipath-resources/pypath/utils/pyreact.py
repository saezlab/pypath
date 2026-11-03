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

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import importlib as imp
import re
import gzip
import tarfile
import zipfile
import struct
import itertools
import bs4
import time
from lxml import etree
try:
    import cPickle as pickle
except:
    import pickle

import pypath.utils.mapping as mapping
import pypath.utils.reflists as reflists
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.share.progress as progress
import pypath.share.curl as curl
import pypath.internals.intera as intera
import pypath.internals.refs as refs
import pypath.inputs.uniprot as uniprot_input
import pypath.resources.urls as urls
import pypath.utils.seq as seq
import pypath.share.session as session_mod
import pypath.share.cache as cache


class BioPaxReader(session_mod.Logger):
    """
    This class parses a BioPAX file and exposes its content easily accessible
    for further processing. First it opens the file, if necessary it extracts
    from the archive. Then an `lxml.etree.iterparse` object is created, so the
    iteration is efficient and memory requirements are minimal. The iterparse
    object is iterated then, and for each tag included in the
    `BioPaxReader.methods` dict, the appropriate method is called. These me-
    thods extract information from the BioPAX entity, and store it in arbit-
    rary data structures: strings, lists or dicts. These are stored in dicts
    where keys are the original IDs of the tags, prefixed with the unique ID
    of the parser object. This is necessary to give a way to merge later the
    result of parsing more BioPAX files. For example, `id42` may identify
    EGFR in one file, but AKT1 in the other. Then, the parser of the first
    file has a unique ID of a 5 letter random string, the second parser a
    different one, and the molecules with the same ID can be distinguished
    at merging, e.g. EGFR will be `ffjh2@id42` and AKT1 will be `tr9gy@id42`.
    The methods and the resulted dicts are named after the BioPAX elements,
    sometimes abbreviated. For example, `BioPaxReader.protein()` processes
    the `<bp:Protein>` elements, and stores the results in
    `BioPaxReader.proteins`.

    In its current state, this class does not parse every information and
    all BioPax entities. For example, nucleic acid related entities and
    interactions are omitted. But these easily can be added with minor mo-
    difications.
    """

    def __init__(
            self,
            biopax,
            source,
            cleanup_period=800,
            file_from_archive=None,
            silent=False,
        ):
        """
        :param str,FileOpener biopax: either a filename, or a FileOpener
        object; if string is supplied, the FileOpener will be created in-
        ternally

        :param str source: the name of the data source, e.g. *Reactome*

        :param int cleanup_period: the number of last elements stored
        during the iteration of lxml.etree.iterparse; lower number
        results lower memory usage, but might risk that an element is
        deleted before it has been processed. Default is 800, which is
        a safe option.

        :param str file_from_archive: in case of processing an archive
        which may contain multiple files (tar.gz or zip), the path of
        the file to be processed needs to be supplied.
        E.g. *BioPax/Homo_sapiens.owl*.

        :param bool silent: whether print status messages and progress
        bars during processing. If you process large number of small
        files, better to set False, in case of one large file, True.
        The default is *False*.
        """

        session_mod.Logger.__init__(self, name = 'biopax')

        self.biopax = biopax
        self.source = source
        self.file_from_archive = file_from_archive
        self.cleanup_period = cleanup_period
        self.biopax_tmp_file = None
        self.cachedir = cache.get_cachedir()
        self.parser_id = common.random_string()
        self.silent = silent

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
        self.bpmodu = '%sModulation' % self.bppref
        self.bpcoma = '%sComplexAssembly' % self.bppref
        self.bptran = '%sTransport' % self.bppref
        self.bptrbr = '%sTransportWithBiochemicalReaction' % self.bppref
        self.bpmoli = '%sMolecularInteraction' % self.bppref
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
        self.bpmerf = '%smemberEntityReference' % self.bppref
        self.bperef = '%sentityReference' % self.bppref
        self.bpxref = '%sxref' % self.bppref
        self.bpdinm = '%sdisplayName' % self.bppref
        self.bprelr = '%sRelationshipXref' % self.bppref
        self.bpcsto = '%scomponentStoichiometry' % self.bppref
        self.bpcomp = '%scomponent' % self.bppref
        self.bpstoc = '%sstoichiometricCoefficient' % self.bppref
        self.bpphye = '%sphysicalEntity' % self.bppref
        self.bpcted = '%scontrolled' % self.bppref
        self.bpcter = '%scontroller' % self.bppref
        self.bpctyp = '%scontrolType' % self.bppref
        self.bppart = '%sparticipant' % self.bppref
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
        self.interactions = {}
        self.transports = {}
        self.transwreas = {}
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
        self.interactions_not2 = {}

        self.methods = {
            self.bpprot: 'protein',
            self.bpcplx: 'cplex',
            self.bpprre: 'pref',
            self.bpuxrf: 'uxref',
            self.bprelr: 'rxref',
            self.bpstoi: 'stoichiometry',
            self.bpreac: 'reaction',
            self.bpcoma: 'cassembly',
            self.bptrbr: 'reaction',
            self.bptran: 'reaction',
            self.bpmoli: 'interaction',
            self.bpcata: 'catalysis',
            self.bpmodu: 'catalysis',
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
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def process(self, silent=False):
        """
        This method executes the total workflow of BioPax processing.

        :param bool silent: whether to print status messages and progress bars.
        """
        self.silent = silent
        self.open_biopax()
        self.biopax_size()
        self.extract()
        self.set_progress()
        self.init_etree()
        if self.bp is not None:

            self.iterate()
            self.close_biopax()

            if len(self.interactions_not2):

                self._log(
                    '%u interactions have not exactly 2 '
                    'participants (%s).' % (
                        len(self.interactions_not2),
                        self.source,
                    )
                )

        else:

            self._log(
                'XML syntax error or empty file encountered. '
                'Skipping to next file or resource.'
            )

    def open_biopax(self):
        """
        Opens the BioPax file. This method should not be called directly,
        ``BioPaxReader.process()`` calls it.
        """
        opener_args = (
            {
                'default_mode': 'rb',
            }
            if self.file_from_archive is None else
            {
                'files_needed': [self.file_from_archive],
            }
        )
        if type(self.biopax) is curl.FileOpener:
            self.opener = self.biopax
        else:
            self.opener = curl.FileOpener(self.biopax, **opener_args)
        if type(self.opener.result) is dict:
            if self.file_from_archive is None or \
                    self.file_from_archive not in self.opener.result:
                self.file_from_archive = \
                    sorted(list(self.opener.result.keys()))[0]
            self._biopax = self.opener.result[self.file_from_archive]
        elif self.opener.type == 'gz':
            self._biopax = self.opener.gzfile
        else:
            self._biopax = self.opener.fileobj

    def biopax_size(self):
        """
        Gets the uncompressed size of the BioPax XML. This is needed in
        order to have a progress bar. This method should not be called
        directly, ``BioPaxReader.process()`` calls it.
        """
        self.bp_filesize = self.opener.sizes[self.file_from_archive] \
            if hasattr(self.opener, 'sizes') else self.opener.size

    def extract(self):
        """
        Extracts the BioPax file from compressed archive. Creates a
        temporary file. This is needed to trace the progress of
        processing, which is useful in case of large files.
        This method should not be called directly,
        ``BioPaxReader.process()`` calls it.
        """

        if self.opener.type != 'plain':

            self.biopax_tmp_file = os.path.join(
                self.cachedir,
                'biopax.processing.tmp.owl',
            )

            self._log(
                'Extracting %s from %s compressed file.' % (
                    self.file_from_archive,
                    self.opener.type
                )
            )

            if not self.silent:
                prg = progress.Progress(
                    self.bp_filesize, 'Extracting %s from %s compressed file' %
                    (self.file_from_archive, self.opener.type), 1000000)

            with open(self.biopax_tmp_file, 'wb') as tmpf:

                while True:

                    chunk = self._biopax.read(100000)
                    if not self.silent:
                        prg.step(len(chunk))
                    if not len(chunk):
                        break
                    if hasattr(chunk, 'encode'):
                        chunk = chunk.encode('utf8')
                    tmpf.write(chunk)

            self._biopax = open(self.biopax_tmp_file, 'rb')

            if not self.silent:
                prg.terminate()

    def init_etree(self):
        """
        Creates the ``lxml.etree.iterparse`` object.
        This method should not be called directly,
        ``BioPaxReader.process()`` calls it.
        """
        try:

            self.bp = etree.iterparse(self._biopax, events=('start', 'end'))
            _, self.root = next(self.bp)

        except etree.XMLSyntaxError:

            self.bp = None

        self.used_elements = []

    def set_progress(self):
        """
        Initializes a progress bar.
        This method should not be called directly,
        ``BioPaxReader.process()`` calls it.
        """
        if not self.silent:
            self.prg = progress.Progress(self.bp_filesize,
                                         'Processing %s from BioPAX XML' %
                                         self.source, 33)

    def iterate(self):
        """
        Iterates the BioPax XML and calls the appropriate methods
        for each element.
        This method should not be called directly,
        ``BioPaxReader.process()`` calls it.
        """
        self.fpos = self._biopax.tell()
        try:
            for ev, elem in self.bp:
                # step the progressbar:
                new_fpos = self._biopax.tell()
                if not self.silent:
                    self.prg.step(new_fpos - self.fpos)
                self.fpos = new_fpos
                self.next_elem = elem
                self.next_event = ev
                self.next_id = '%s@%s' % \
                    (self.parser_id, self.next_elem.get(self.rdfid)) \
                    if self.rdfid in elem.attrib \
                    else '%s@%s' % \
                    (self.parser_id, self.next_elem.get(self.rdfab))
                if ev == 'end' and self.next_elem.tag in self.methods:
                    method = getattr(self, self.methods[self.next_elem.tag])
                    method()
                    self.root.clear()
                # self.used_elements.append(self.next_elem)
                # self.cleanup_hook()
        except etree.XMLSyntaxError as e:
            if not self.silent:
                self.prg.terminate(status='failed')
            sys.stdout.write('\n\t:: Syntax error in BioPAX:\n\t\t%s\n' %
                             str(e))
            sys.stdout.flush()
        if not self.silent:
            self.prg.terminate()

    def cleanup_hook(self):
        """
        Removes the used elements to free up memory.
        This method should not be called directly,
        ``BioPaxReader.iterate()`` calls it.
        """
        if len(self.used_elements) > self.cleanup_period:
            for _ in xrange(int(self.cleanup_period / 2)):
                e = self.used_elements.pop()
                e.clear()

    def close_biopax(self):
        """
        Deletes the iterator and closes the file object.
        This method should not be called directly,
        ``BioPaxReader.process()`` calls it.
        """
        del self.bp
        self._biopax.close()

    def protein(self):
        entref = self.next_elem.find(self.bperef)
        if entref is not None:
            protein = self.get_none(entref.get(self.rdfres)).replace('#', '')
            self.proteins[self.next_id] = {
                'protein': '%s@%s' % (self.parser_id, protein),
                'seqfeatures': self._bp_collect_resources(self.bpfeat),
                'modfeatures': self._bp_collect_resources(self.bpfeat)
            }
        else:
            self.pfamilies[self.next_id] = \
                self._bp_collect_resources(self.bpmphe)

    def pref(self):
        self.prefs[self.next_id] = {}
        self.prefs[self.next_id]['uxrefs'] = \
            self._bp_collect_resources(self.bpxref)
        self.prefs[self.next_id]['prefs'] = \
            self._bp_collect_resources(self.bpmerf)

    def uxref(self):
        db = self.next_elem.find(self.bpdb)
        id_type = db.text.lower()
        i = self.next_elem.find(self.bpid)
        if i is not None:
            self.ids[self.next_id] = (id_type, i.text)

    def rxref(self):
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
        elif self.next_elem.find(self.bpmphe) is not None:
            self.cvariations[self.next_id] = \
                self._bp_collect_resources(self.bpmphe)
        elif self.next_elem.find(self.bpcomp) is not None:
            self.complexes[self.next_id] = \
                self._bp_collect_resources(self.bpcomp)

    def stoichiometry(self):
        snum = self.next_elem.find(self.bpstoc).text
        self.stoichiometries[self.next_id] = ('%s@%s' % (
            self.parser_id,
            self.next_elem.find(self.bpphye).get(self.rdfres).replace(
                '#', '')), int(float(snum)))

    def interaction(self):
        part = self._bp_collect_resources(self.bpleft)
        if len(part) == 2:
            self.interactions[self.next_id] = {
                'refs': self._bp_collect_resources(self.bpxref),
                'left': [part[0]],
                'right': [part[1]],
                'type': self.next_elem.tag.split('}')[-1]
            }
        else:
            self.interactions_not2[self.next_id] = len(part)

    def reaction(self):
        self.reactions[self.next_id] = {
            'refs': self._bp_collect_resources(self.bpxref),
            'left': self._bp_collect_resources(self.bpleft),
            'right': self._bp_collect_resources(self.bprgth),
            'type': self.next_elem.tag.split('}')[-1]
        }

    def cassembly(self):
        self.cassemblies[self.next_id] = {
            'refs': self._bp_collect_resources(self.bpxref),
            'left': self._bp_collect_resources(self.bpleft),
            'right': self._bp_collect_resources(self.bprgth),
            'type': self.next_elem.tag.split('}')[-1]
        }

    def catalysis(self):
        cter = self.next_elem.find(self.bpcter)
        cted = self.next_elem.find(self.bpcted)
        if cter is not None and cted is not None:
            typ = self.next_elem.find(self.bpctyp)
            self.catalyses[self.next_id] = {
                'controller': '%s@%s' %
                (self.parser_id, self.get_none(cter.get(self.rdfres))),
                'controlled': '%s@%s' %
                (self.parser_id, self.get_none(cted.get(self.rdfres))),
                'type': '' if typ is None else typ.text
            }

    def control(self):
        cter = self.next_elem.find(self.bpcter)
        cted = self.next_elem.find(self.bpcted)
        if cter is not None and cted is not None:
            typ = self.next_elem.find(self.bpctyp)
            self.controls[self.next_id] = {
                'refs': self._bp_collect_resources(self.bpxref),
                'controller': '%s@%s' %
                (self.parser_id, cter.get(self.rdfres).replace('#', '')),
                'controlled': '%s@%s' %
                (self.parser_id, cted.get(self.rdfres).replace('#', '')),
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
        self.fragfeas[self.next_id] = '%s@%s' % \
            (self.parser_id, felo.replace('#', ''))

    def seqint(self):
        beg = self.next_elem.find(self.bpibeg)
        end = self.next_elem.find(self.bpiend)
        self.seqints[self.next_id] = (
            '%s@%s' % (self.parser_id, beg.get(self.rdfres).replace('#', ''))
            if beg is not None else None,
            '%s@%s' % (self.parser_id, end.get(self.rdfres).replace('#', ''))
            if end is not None else None)

    def seqsite(self):
        seqp = self.next_elem.find(self.bpseqp)
        if seqp is not None and seqp.text is not None:
            self.seqsites[self.next_id] = int(seqp.text)

    def modfea(self):
        felo = self.next_elem.find(self.bpfelo)
        moty = self.next_elem.find(self.bpmoty)
        if felo is not None and moty is not None:
            self.modfeas[self.next_id] = (
                '%s@%s' % (
                    self.parser_id,
                    self.next_elem.find(self.bpfelo).get(self.rdfres).replace(
                        '#', '')), '%s@%s' %
                (self.parser_id,
                 self.next_elem.find(self.bpmoty).get(self.rdfres).replace(
                     '#', '')))

    def seqmodvoc(self):
        term = self.next_elem.find(self.bpterm)
        if term is not None:
            self.seqmodvocs[self.next_id] = term.text

    def pathway(self):
        name = self.next_elem.find(self.bpdinm)
        if name is not None:
            name = name.text
        try:
            self.pathways[self.next_id] = {
                'components': self._bp_collect_resources(self.bppcom),
                'name': name
            }
        except TypeError:
            sys.stdout.write('Wrong type at element:\n')
            sys.stdout.write('%s %s' %
                             (str(etree.tostring(self.next_elem))[:76], '...'))
            sys.stdout.flush()

    def get_none(self, something):
        if something is not None:
            return something.replace('#', '')
        return something

    def _bp_collect_resources(self, tag, restype=None):
        return \
            list(
                map(
                    lambda e:
                    '%s@%s' %
                    (self.parser_id, e.get(self.rdfres).replace('#', '')),
                    filter(
                        lambda e:
                        self.rdfres in e.attrib and (
                            restype is None or
                            e.get(self.rdfres).replace('#', '')
                            .startswith(restype)
                        ),
                        self.next_elem.iterfind(tag)
                    )
                )
            )


class AttributeHandler(object):
    def __init__(self):
        pass

    def add_source(self, source):
        if type(source) in _const.CHAR_TYPES:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
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


class Entity(AttributeHandler):
    def __init__(self, identifier, id_type, sources=[], attrs=None):
        super(Entity, self).__init__()
        self.id = identifier
        self.id_type = id_type
        self.sources = set([])
        self.attrs = {}
        self.add_source(sources)
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

    def __iter__(self):
        """
        With this method it is possible to iterate ``Entity`` objects
        just like ``EntitySet`` objects.

        Yields the object.
        """
        for i in [self]:
            yield i

    def expand(self):
        """
        With this method it is possible to iterate ``Entity`` objects
        just like ``EntitySet`` objects.

        Yields string.
        """
        for i in [self.id]:
            yield i

    def itermembers(self):
        for m in self.__iter__():
            yield (m, {})


class Protein(Entity):
    def __init__(self, protein_id, id_type='uniprot', sources=[], attrs=None):
        super(Protein, self).__init__(
            protein_id, id_type, sources=sources, attrs=attrs)
        self.type = 'protein'

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def key(self):
        return self.id

    def proteins(self):
        return set([self.id])

    def __ror__(self, other):
        return self.proteins() | other

    def __rand__(self, other):
        return self.proteins() & other

    def __rsub__(self, other):
        return self.proteins() - other


class Reference(Entity):
    def __init__(self, ref_id, id_type='pubmed', sources=[]):
        super(Reference, self).__init__(ref_id, id_type, sources=sources)

    def get_ref(self):
        return refs.Reference(self.ref_id)


class EntitySet(AttributeHandler):
    def __init__(self, members, sources=[], sep=';', parent=None):
        super(EntitySet, self).__init__()
        self.parent = parent
        self.members = sorted(common.unique_list(members))
        self.set = set(self.members)
        self.sources = set([])
        self.attrs = {}
        self.add_source(sources)
        self.originals = {}
        self.type = None
        self.sep = sep

    def __str__(self):
        return self.sep.join(sorted(map(lambda x: str(x), list(self.members))))

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


class Intersecting(object):
    def __init__(self):
        pass

    def __eq__(self, other):
        return self.type == other.type and \
            not self.set.isdisjoint(other.set)


class Complex(EntitySet):
    def __init__(self, members, source, parent=None):
        super(Complex, self).__init__(members, source, parent=parent)
        self.type = 'complex'

    def __str__(self):
        return '<%s>' % EntitySet.__str__(self)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def get_stoichiometries(self, source, cid, with_pids=False):
        if source in self.sources and cid in self.attrs[source]:
            return \
                list(
                    map(
                        lambda memb:
                            (memb[0], memb[1]['stoi'], memb[1]['pid'])
                        if with_pids else (memb[0], memb[1]['stoi']),
                        filter(
                            lambda memb:
                                type(memb[1]) is dict
                                and 'stoi' in memb[1]
                                and 'pid' in memb[1],
                            iteritems(self.attrs[source][cid])
                        )
                    )
                )

    def expand(self):
        for i, m1 in enumerate(self.members):
            for m2 in self.members[i + 1:]:
                yield m1, m2

    def itermembers(self):
        for m in self.members:
            return self.parent.proteins[m]

    def key(self):
        return tuple(self.members)

    def proteins(self):
        return set(self.members)

    def __ror__(self, other):
        return self.proteins() | other

    def __rand__(self, other):
        return self.proteins() & other

    def __rsub__(self, other):
        return self.proteins() - other


class ProteinFamily(Intersecting, EntitySet):
    def __init__(self, members, source, parent=None):
        EntitySet.__init__(self, members, source, parent=parent)
        Intersecting.__init__(self)
        self.type = 'pfamily'

    def __str__(self):
        return '|%s|' % EntitySet.__str__(self)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def _itermembers1(self):
        """
        Iterates protein family, yields proteins.
        """
        for m in self.members:
            if m in self.parent.proteins:
                yield self.parent.proteins[m]
            else:
                this_attrs = dict(map(lambda s:
                                      (s, {'pids': {},
                                           'prefs': {},
                                           'originals': {}}),
                                      self.source))
                for s, a in iteritems(self.attrs):
                    for pfid, aa in iteritems(a):
                        this_attrs[s]['pids'][pfid] = {}
                yield Protein(m, sources=self.sources, attrs=this_attrs)

    def expand(self):
        for pkey in self.__iter__():
            yield self.parent.proteins[pkey]

    def itermembers(self):
        for m in self.expand():
            attrs = \
                dict(
                    map(
                        lambda s:
                            (
                                s,
                                dict(
                                    map(
                                        lambda pf:
                                            (pf[0], pf[1][m.id]['pid']),
                                        iteritems(self.attrs[s])
                                    )
                                )
                            ),
                        self.sources
                    )
                )
            yield (m, attrs)

    def key(self):
        return tuple(self.members)

    def proteins(self):
        return set(self.members)

    def __ror__(self, other):
        return self.proteins() | other

    def __rand__(self, other):
        return self.proteins() & other

    def __rsub__(self, other):
        return self.proteins() - other


class ComplexVariations(Intersecting, EntitySet):
    def __init__(self, members, source, parent=None):
        EntitySet.__init__(self, members, source, sep='|', parent=parent)
        Intersecting.__init__(self)
        self.type = 'cvariations'

    def __str__(self):
        return '<%s>' % EntitySet.__str__(self)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def expand(self):
        for m in self.members:
            for m1, m2 in m.expand():
                yield m1, m2

    def itermembers(self):
        """
        This is a convenient iterator for the expand methods
        of higher classes like ``ReactionSide`` or ``Control``.
        """
        for m in self.__iter__():
            attrs = dict(
                map(
                    lambda s:
                        (
                            s,
                            reduce(
                                lambda c1, c2:
                                    c1 | c2,
                                map(
                                    lambda cid:
                                        m.attrs[s][cid]['children']
                                    if
                                        'children' in m.attrs[s][cid]
                                    else
                                        set([cid]),
                                    filter(
                                        lambda cid:
                                            cid in m.attrs[s],
                                        self.attrs[s]['cids']
                                    )
                                )
                            )
                        ),
                    self.sources
                )
            )
            yield (m, attrs)

    def key(self):
        return self.__str__()

    def proteins(self):
        return \
            reduce(
                lambda m1, m2:
                    m1 | m2,
                self.members,
                set([])
            )

    def __ror__(self, other):
        return self.proteins() | other

    def __rand__(self, other):
        return self.proteins() & other

    def __rsub__(self, other):
        return self.proteins() - other


class PyReact(session_mod.Logger):

    def __init__(
        self,
        ncbi_tax_id=9606,
        default_id_types={},
        modifications=True,
        seq=None,
        silent=False,
        max_complex_combinations=100,
        max_reaction_combinations=100,
    ):


        self.cachedir = cache.get_cachedir()

        self.ncbi_tax_id = ncbi_tax_id
        self.modifications = modifications
        self.parsers = {}
        self.sources = set([])
        self.seq = seq
        self.mapper = mapping
        self.silent = silent
        self.max_complex_combinations = max_complex_combinations
        self.max_reaction_combinations = max_reaction_combinations
        self.slow_complexes = {}
        self.huge_complexes = {}
        self.huge_reactions = {}

        self.refs = {}
        self.species = {}
        self.proteins = {}
        self.pfamilies = {}
        self.complexes = {}
        self.cvariations = {}
        self.reactions = {}
        self.controls = {}
        self.mods = {}
        self.frags = {}

        self.rrefs = {}
        self.rproteins = {}
        self.rpfamilies = {}
        self.rcomplexes = {}
        self.rcvariations = {}
        self.rreactions = {}
        self.rcontrols = {}
        self.references = {}

        self.id_types = {}
        self.default_id_types = {'protein': 'uniprot'}
        self.default_id_types.update(default_id_types)

        session_mod.Logger.__init__(self, name = 'pyreact')

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def load_reactome(self):

        self._log('Loading Reactome.')

        def reactome_id_proc(_id):
            _id = _id.split('-')
            return {
                'id': _id[0],
                'isoform': int(_id[1]) if len(_id) > 1 else 1
            }

        biopax = curl.Curl(
            urls.urls['reactome']['biopax_l3'],
            large=True,
            silent=self.silent,
        )
        parser = BioPaxReader(
            biopax.outfile, 'Reactome',
            file_from_archive='Homo_sapiens.owl'
        )
        parser.process()

        self.add_dataset(
            parser,
            id_types={'uniprot isoform': 'uniprot'},
            process_id=reactome_id_proc)

    def load_acsn(self):

        self._log('Loading ACSN.')

        biopax = curl.Curl(
            urls.urls['acsn']['biopax_l3'],
            large=True,
            silent=self.silent,
            default_mode = 'rb',
        )

        self.parser = BioPaxReader(biopax.outfile, 'ACSN')
        self.parser.process()

        self.add_dataset(self.parser, id_types={'hgnc': 'genesymbol'})
        del self.parser

    def load_kegg(self):

        self._log('Loading KEGG.')

        biopax = curl.Curl(
            urls.urls['kegg_pws']['biopax_l3'],
            large=True,
            silent=self.silent,
        )

        parser = BioPaxReader(biopax.outfile, 'KEGG')
        parser.process()

        self.add_dataset(parser, {'uniprot knowledgebase': 'uniprot'})

    def load_pid(self):

        self._log('Loading NCI-PID.')

        biopax = curl.Curl(
            urls.urls['nci-pid']['biopax_l3'],
            large=True,
            silent=self.silent,
        )

        parser = BioPaxReader(biopax.outfile, 'NCI-PID')
        parser.process()

        self.add_dataset(parser)


    def load_wikipathways(self):

        self._log('Loading WikiPathways.')

        biopaxes = curl.Curl(
            urls.urls['wikipw']['biopax_l3'],
            large=True,
            silent=self.silent
        )

        if not self.silent:
            prg = progress.Progress(
                len(biopaxes.result),
                'Processing multiple BioPAX files',
                1,
                percent=False,
            )

        silent_default = self.silent
        self.silent = True

        for fname in biopaxes.files_multipart.keys():

            if not silent_default:
                prg.step()
            parser = BioPaxReader(
                biopaxes.outfile, 'WikiPathways', file_from_archive=fname)
            parser.process(silent=True)

            self.add_dataset(
                parser,
                id_types={
                    'ensembl': 'ensg',
                    'entrez gene': 'entrez',
                    'hgnc': 'genesymbol'
                })

        if not silent_default:
            prg.terminate()
        self.silent = silent_default


    def load_panther(self):

        self._log('Loading PANTHER.')

        biopaxes = curl.Curl(
            urls.urls['panther']['biopax_l3'], large=True, silent=self.silent)
        if not self.silent:
            prg = progress.Progress(
                len(biopaxes.files_multipart),
                'Processing multiple BioPAX files',
                1,
                percent=False)

        silent_default = self.silent
        self.silent = True

        for fname in biopaxes.files_multipart.keys():

            if not silent_default:
                prg.step()
            parser = BioPaxReader(
                biopaxes.outfile, 'PANTHER', file_from_archive=fname)
            parser.process(silent=True)

            self.add_dataset(
                parser,
                id_types={
                    'ensembl': 'ensg',
                    'entrez gene': 'entrez',
                    'hgnc': 'genesymbol'
                })

        if not silent_default:
            prg.terminate()
        self.silent = silent_default


    def load_netpath(self):

        self._log('Loading NetPath.')

        names = self.netpath_names()

        if not self.silent:
            prg = progress.Progress(
                len(names),
                'Processing multiple BioPAX files',
                1,
                percent=False)

        silent_default = self.silent
        self.silent = True

        for pwnum in names.keys():

            if not silent_default:
                prg.step()
            biopax = curl.Curl(
                urls.urls['netpath_bp']['biopax_l3'] % int(pwnum), silent=True)
            parser = BioPaxReader(biopax.outfile, 'NetPath')
            parser.process(silent=True)
            self.add_dataset(
                parser,
                id_types={
                    'ensembl': 'ensg',
                    'entrez gene': 'entrez',
                    'hgnc': 'genesymbol'
                })

        if not silent_default:
            prg.terminate()
        self.silent = silent_default


    def netpath_names(self):

        repwnum = re.compile(r'_([0-9]+)$')
        result = {}
        url = urls.urls['netpath_names']['url']
        c = curl.Curl(url)
        html = c.result
        soup = bs4.BeautifulSoup(html, 'html.parser')
        for a in soup.find_all('a'):
            if a.attrs['href'].startswith('pathways'):
                num = repwnum.findall(a.attrs['href'])[0]
                name = a.text
                result[num] = name
        return result


    def load_all(self):

        self._log('Loading all databases.')

        self.load_wikipathways()
        self.load_netpath()
        self.load_panther()
        self.load_acsn()
        self.load_pid()
        self.load_reactome()


    def add_dataset(self, parser, id_types={}, process_id=lambda x: {'id': x}):

        self.id_types.update(id_types)
        self.source = parser.source
        self.parser = parser
        self.id_processor = process_id
        self.merge()


    def merge(self):

        if self.source not in self.huge_complexes:
            self.huge_complexes[self.source] = {}
        if self.source not in self.slow_complexes:
            self.slow_complexes[self.source] = {}
        if self.source not in self.huge_reactions:
            self.huge_reactions[self.source] = {}

        if not self.silent:
            self.prg = progress.Progress(
                12, 'Processing %s' % self.source, 1, percent=False)
        self.sources.add(self.source)
        if self.source not in self.parsers:
            self.parsers[self.source] = {}
        self.parsers[self.source][self.parser.parser_id] = self.parser
        self.set_corrections()

        if not self.silent:
            self.prg.step(status='processing references')
        self.merge_refs()

        if not self.silent:
            self.prg.step(status='processing proteins')
        self.merge_proteins()

        if not self.silent:
            self.prg.step(status='processing protein families')
        self.merge_pfamilies()

        if not self.silent:
            self.prg.step(status='processing protein modifications')
        if self.modifications:
            self.merge_modifications()

        if not self.silent:
            self.prg.step(status='processing complexes')
        remaining = self.merge_complexes()

        if not self.silent:
            self.prg.step(status='processing complex variations')
        self.merge_cvariations()

        if not self.silent:
            self.prg.step(status='processing complexes 2')
        remaining = self.merge_complexes(this_round=remaining)

        if not self.silent:
            self.prg.step(status='generating complex variations')
        self.gen_cvariations()

        if not self.silent:
            self.prg.step(status='processing reactions')
        self.merge_reactions()

        if not self.silent:
            self.prg.step(status='processing complex asseblies')
        self.merge_cassemblies()

        if not self.silent:
            self.prg.step(status='processing controls')
        self.merge_controls()

        if not self.silent:
            self.prg.step(status='processing catalyses')
        self.merge_catalyses()

        if not self.silent:
            self.prg.terminate()

            sys.stdout.write('\t:: %u proteins, %u complexes'
                             ', %u reactions and %u'
                             ' controls have been added.\n' %
                             (self.proteins_added, self.complexes_added,
                              self.reactions_added + self.cassemblies_added,
                              self.controls_added + self.catalyses_added))
            sys.stdout.write('\t:: %u complexes have not been expanded '
                             'because number of combinations larger than %u.\n'
                             % (len(self.huge_complexes[self.source]),
                                self.max_complex_combinations))
            sys.stdout.write('\t:: %u complexes took longer '
                             'to process than 5 seconds.\n' %
                             (len(self.slow_complexes[self.source])))
            sys.stdout.write('\t:: %u reactions have not been expanded '
                             'because number of combinations larger than %u.\n'
                             % (len(self.huge_reactions[self.source]),
                                self.max_reaction_combinations))
            sys.stdout.write('\t:: Access them in `huge_complexes`, '
                             '`slow_complexes` and `huge_reactions`.\n')
            sys.stdout.flush()

        self.remove_defaults()

    def remove_defaults(self):
        self.parser = None
        self.source = None

    def set_corrections(self):
        # because not all can follow the standards...
        if self.source == 'ACSN':
            self.pref_correction = lambda l: filter(lambda e: e[6:10] == 'HUGO', l)
        else:
            self.pref_correction = lambda l: l
        if self.source == 'ACSN':
            self.pref_refs = lambda l: filter(lambda e: e[6:12] == 'PubMed', l)
        else:
            self.pref_refs = lambda l: []
        if self.source in [
                'ACSN', 'WikiPathways', 'NetPath', 'PANTHER', 'NCI-PID', 'KEGG'
        ]:
            self.ambiguous_ids_permitted = True
        else:
            self.ambiguous_ids_permitted = False

    def merge_refs(self):

        self.refs_added = 0

        if self.source not in self.rrefs:
            self.rrefs[self.source] = {}

        for refid, pubmed in iteritems(self.parser.pubrefs):

            ref = Reference(pubmed, sources=self.source)
            ref.attrs[self.source]['refid'] = set([])
            ref.attrs[self.source]['refid'].add(refid)

            self.refs_added += 1

            if pubmed in self.refs:
                self.refs[pubmed] += ref
            else:
                self.refs[pubmed] = ref

            self.rrefs[self.source][refid] = pubmed

    def merge_proteins(self):
        def get_protein_ids(pref):

            pids = []
            refs = []
            if pref in self.parser.prefs:
                uxrefs = self.pref_correction(self.parser.prefs[pref][
                    'uxrefs'])
                pids = \
                    common.unique_list(
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
                refids = self.pref_refs(self.parser.prefs[pref]['uxrefs'])
                refs = \
                    common.unique_list(
                        map(
                            lambda refid:
                                self.parser.ids[refid],
                            filter(
                                lambda refid:
                                    refid in self.parser.ids,
                                refids
                            )
                        )
                    )
                # in panther proteinreferences are children of
                # other proteinreferences...
                for subpref in self.parser.prefs[pref]['prefs']:
                    subpids, subrefs = get_protein_ids(subpref)
                    pids.extend(subpids)
                    refs.extend(subrefs)
            return pids, refs

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
                if id_a['id'] is not None:
                    target_ids.extend(
                        self.mapper.map_name(id_a['id'], std_id_type,
                                             self.default_id_types['protein']))
            target_ids = filter(
                lambda p: reflists.check(
                    p,
                    self.default_id_types['protein'],
                    self.ncbi_tax_id
                ),
                target_ids)
            target_ids = common.unique_list(target_ids)
            if len(target_ids) > 1:
                if not self.ambiguous_ids_permitted:
                    sys.stdout.write('\t:: Ambiguous ID '
                                     'translation: from %s to %s\n' %
                                     (ids, target_ids))
            elif len(target_ids) == 0:
                target_ids = None
            else:
                target_ids = list(target_ids)[0]
            return target_ids, id_attrs

        self.rproteins[self.source] = {}
        self.proteins_added = 0
        self.pfamilies_added = 0

        for pid, p in iteritems(self.parser.proteins):
            ids, pubmeds = get_protein_ids(p['protein'])
            target_id, id_attrs = map_protein_ids(ids)
            if target_id is None:
                continue
            if type(target_id) is list:
                # go for a protein family:
                self.add_pfamily(list(map(lambda t: (t, pid), target_id)), pid)
                continue
            attrs = {
                self.source: {
                    'prefs': set([p['protein']]),
                    'pids': {
                        pid: {}
                    },
                    'refs': set(pubmeds),
                    'originals': set([])
                }
            }
            for original_id, id_a in iteritems(id_attrs):
                attrs[self.source]['originals'].add(original_id)
                for k, v in iteritems(id_a):
                    if k != 'id':
                        attrs[self.source]['pids'][pid][k] = set([])
                        attrs[self.source]['pids'][pid][k].add(v)
            protein = Protein(
                target_id, sources=set([self.source]), attrs=attrs)
            self.proteins_added += 1
            if target_id in self.proteins:
                self.proteins[target_id] += protein
            else:
                self.proteins[target_id] = protein
            self.rproteins[self.source][pid] = target_id

    def preprocess_seqmodvoc(self):
        self.seqmod_dict = {}
        if self.source in common.mod_keywords:
            kws = common.mod_keywords[self.source]
            aas = common.aanames
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
                return int(float(self.parser.seqsites[seqsite]))

        def get_residue(protein, isof, resnum, resname):
            if protein in self.seq:
                if isof is not None and isof in self.seq[protein].isof:
                    sresname = self.seq[protein].get(resnum, isoform=isof)
                    if sresname == resname or resname is None:
                        return sresname, isof
                for isof in sorted(self.seq[protein].isoforms()):
                    sresname = self.seq[protein].get(resnum, isoform=isof)
                    if sresname == resname or resname is None:
                        return sresname, isof
            return resname, isof

        self.fragfeatures_added = 0
        self.modfeatures_added = 0

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
                                            self.frags[frag].add_evidences(
                                                self.source
                                            )
                                        else:
                                            instance = \
                                                self.seq[protein].get(
                                                    start, end, isof)
                                            mot = intera.Motif(
                                                protein,
                                                start,
                                                end,
                                                isoform=isof,
                                                instance=instance,
                                                evidences=self.source,
                                            )
                                            self.frags[frag] = mot
                                        self.proteins[protein].update_attr([
                                            self.source, 'pids', pid, 'frags',
                                            {frag}
                                        ])
                                        self.fragfeatures_added += 1

                if modfea in self.parser.modfeas:
                    resnum = get_seqsite(self.parser.modfeas[modfea][0])
                    seqmodvoc = self.parser.modfeas[modfea][1]
                    if seqmodvoc in self.seqmod_dict:
                        typ, resname = self.seqmod_dict[seqmodvoc]
                        for protein, isof in proteins:
                            if protein in self.seq and resnum is not None:
                                resname, isof = \
                                    get_residue(protein, isof, resnum, resname)
                                mod = (protein, isof, resnum, resname, typ)
                                if mod in self.mods:
                                    self.mods[mod].add_evidences(self.source)
                                else:
                                    res = intera.Residue(
                                        resnum, resname, protein, isoform=isof)
                                    start, end, instance = self.seq[
                                        protein].get_region(
                                            resnum, isoform=isof)
                                    mot = intera.Motif(
                                        protein,
                                        start,
                                        end,
                                        isoform=isof,
                                        instance=instance)
                                    ptm = intera.Ptm(
                                        protein,
                                        motif=mot,
                                        residue=res,
                                        evidences=self.source,
                                        isoform=isof,
                                        typ=typ,
                                    )
                                    self.mods[mod] = ptm
                                try:
                                    self.proteins[protein].update_attr([
                                        self.source, 'pids', pid, 'mods',
                                        set([mod])
                                    ])
                                except:
                                    print(protein, pid)
                                self.modfeatures_added += 1

    def merge_pfamilies(self):

        if self.source not in self.rpfamilies:
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
                                iteritems(self.pfamilies[spf].attrs[
                                          self.source][spfid])
                            )
                        )
                    proteins.extend(spfmembs)

                self.add_pfamily(proteins, pfid)

            this_round = next_round
            next_round = []

    def add_pfamily(self, proteins, pfid):
        if self.source not in self.rpfamilies:
            self.rpfamilies[self.source] = {}

        members = sorted(common.unique_list(map(lambda p: p[0], proteins)))

        # this necessary if we add protein family because of
        # ambiguous id mapping; we want to make sure protein
        # exists for each member of the family.
        for m in members:
            if m not in self.proteins:
                p = Protein(m, sources=self.source)
                p.attrs[self.source]['pids'] = {}
                p.attrs[self.source]['pids'][pfid] = {}
                self.proteins[m] = p
        if len(members):
            pf = ProteinFamily(members, source=self.source, parent=self)
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
            self.pfamilies_added += 1

    def merge_complexes(self, this_round=None):
        """
        Merges complexes from the active ``BioPaxReader`` object.
        Protein families and subcomplexes are expanded, and all combinations
        are created as separate complexes. The complexes from the same ID
        are added to sets in the ``rcomplexes`` dict.
        """

        self.complexes_added = 0
        if self.source not in self.rcomplexes:
            self.rcomplexes[self.source] = {}

        no_protein = set([])
        this_round = set(list(self.parser.complexes.keys())) \
            if this_round is None else this_round
        next_round = []
        prev_round = -1

        while len(this_round) - prev_round != 0:
            prev_round = len(this_round)
            for cid in this_round:
                start_time = time.time()
                stois = self.parser.complexes[cid]

                if len(self.parser.stoichiometries):
                    pids = list(
                        map(lambda stoi: self.parser.stoichiometries[stoi],
                            stois))
                else:
                    pids = list(map(lambda comp: (comp, 1), stois))

                subc_unproc = \
                    any(
                        map(
                            lambda pid:
                                (pid[0] in self.parser.complexes or
                                    pid[0] in self.parser.cvariations)
                            and pid[0] not in self.rcomplexes[self.source]
                            and pid[0] not in no_protein,
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
                pfamilies = \
                    list(
                        map(
                            lambda pfid:
                                list(
                                    map(
                                        lambda memb:
                                            (memb[0], pfid[1], memb[1]['pid']),
                                        iteritems(
                                            self.pfamilies[
                                                self.rpfamilies
                                                [self.source][pfid[0]]]
                                                .attrs[self.source][pfid[0]
                                                                    ]
                                        )
                                    )
                                ),
                            filter(
                                lambda pfid:
                                    pfid[0] in self.rpfamilies[self.source],
                                pids
                            )
                        )
                    )

                pfnum = 0
                if len(pfamilies):
                    pfnum = reduce(lambda pf1l, pf2l: pf1l * pf2l,
                                   map(lambda pf: len(pf), pfamilies))

                if pfnum > self.max_complex_combinations:
                    self.huge_complexes[self.source][cid] = pfnum
                    continue

                subcplexs = \
                    list(
                        map(
                            lambda scid:
                                map(
                                    lambda memb:
                                        (memb, scid[1], scid[0]),
                                    self.rcomplexes[self.source][scid[0]]
                                ),
                            filter(
                                lambda scid:
                                    scid[0] in self.rcomplexes[self.source],
                                pids
                            )
                        )
                    )

                if len(subcplexs):

                    subcplexs = itertools.product(*subcplexs)
                    subcmembs = []

                    for this_subcplex in subcplexs:
                        for sckey, scstoi, scid in this_subcplex:
                            if scid not in no_protein:
                                sc = self.complexes[sckey]
                                scmembs = sc.get_stoichiometries(
                                    self.source, scid, with_pids=True)
                                scmembs = list(
                                    map(lambda p: (p[0], p[1] * scstoi, p[2]),
                                        scmembs))
                                subcmembs.append(scmembs)

                else:
                    subcmembs = [[]]

                if len(subcmembs) * pfnum > self.max_complex_combinations:
                    self.huge_complexes[self.source][cid] = \
                        len(subcmembs) * pfnum
                    continue

                if len(proteins) or len(pfamilies) or \
                        type(subcplexs) is not list:

                    if not len(pfamilies):
                        pfamilies = [[]]
                    else:
                        pfamilies = itertools.product(*pfamilies)

                    for pfamily in pfamilies:

                        for subc in subcmembs:

                            this_proteins = \
                                proteins + list(pfamily) + list(subc)

                            members = sorted(
                                common.unique_list(
                                    map(lambda p: p[0], this_proteins)))
                            if not len(members):
                                continue
                            cplex = Complex(
                                members, source=self.source, parent=self)

                            members = tuple(members)
                            cplex.attrs[self.source][cid] = {}
                            for protein, stoi, pid in this_proteins:
                                cplex.attrs[self.source][cid][protein] = {}
                                cplex.attrs[self.source][cid][protein][
                                    'pid'] = pid
                                cplex.attrs[self.source][cid][protein][
                                    'stoi'] = stoi
                            if members not in self.complexes:
                                self.complexes[members] = cplex
                            else:
                                self.complexes[members] += cplex
                            self.complexes_added += 1
                            if cid not in self.rcomplexes[self.source]:
                                self.rcomplexes[self.source][cid] = set([])
                            self.rcomplexes[self.source][cid].add(members)

                else:
                    no_protein.add(cid)

                elapsed = time.time() - start_time

                if elapsed > 5:
                    self.slow_complexes[self.source][cid] = elapsed

            this_round = next_round
            next_round = []

        return this_round

    def merge_cvariations(self):
        """
        This processes those complexes which are in fact a set of complex variations.
        As simple complexes also are always extended to complex variations because
        they might have not only simple proteins but protein families as members,
        here we only add new records to the attributes of already existing complexes.
        After ``merge_complexes`` will be called again, to process those simple
        complexes which have any of the complex variations processed here among their
        subcomplexes.
        """
        self.cvariations_added = 0
        for cvid, cv in iteritems(self.parser.cvariations):

            cplexes = \
                dict(
                    map(
                        lambda cid:
                            (cid, self.rcomplexes[self.source][cid]),
                        filter(
                            lambda cid:
                                cid in self.rcomplexes[self.source],
                            cv
                        )
                    )
                )

            for cid, ckeys in iteritems(cplexes):
                for ckey in ckeys:
                    c = self.complexes[ckey]
                    if cvid not in c.attrs[self.source] and cid in c.attrs[
                            self.source]:
                        c.attrs[self.source][cvid] = {'children': set([])}
                    if cid in c.attrs[self.source]:
                        c.attrs[self.source][cvid]['children'].add(cid)
                        c.attrs[self.source][cvid].update(c.attrs[self.source][
                            cid])

            if len(cplexes):
                self.rcomplexes[self.source][cvid] = \
                    reduce(
                        lambda c1, c2:
                            c1 | c2,
                        cplexes.values()
                )
                self.cvariations_added += 1

    def gen_cvariations(self):
        """
        Because one key from the BioPax file might represent more complexes,
        *complexvariations* are created to give a way to represent sets of
        combinations. These are created for all complexes, even with only
        one unambiguous constitution. The keys are the constitutions of
        all the combinations listed in alphabetic order, separated by ``|``.
        For example, ``A,B,C|A,B,D|A,B,E``.
        """

        self.rcvariations[self.source] = {}

        for cid, keys in iteritems(self.rcomplexes[self.source]):

            membs = map(lambda key: self.complexes[key], keys)

            cvar = ComplexVariations(membs, source=self.source, parent=self)
            cvar.attrs[self.source]['cids'] = set([cid])

            key = cvar.__str__()
            if key in self.cvariations:
                self.cvariations[key] += cvar
            else:
                self.cvariations[key] = cvar

            self.rcvariations[self.source][cid] = key

    def merge_reactions(self):
        self.reactions_added = 0
        self._merge_reactions(('reactions', 'reaction'))

    def merge_cassemblies(self):
        self.cassemblies_added = 0
        self._merge_reactions(('cassemblies', 'cassembly'))

    def _merge_reactions(self, rclass):
        """
        Merges reaction type entities from the active parser.
        Here protein families and complex variations are not expanded.
        """

        if self.source not in self.rreactions:
            self.rreactions[self.source] = {}

        def get_side(ids):
            members = []
            memb_ids = {}
            for _id in ids:
                for cls in ('proteins', 'pfamilies', 'cvariations'):
                    r = getattr(self, 'r%s' % cls)[self.source]
                    if _id in r:
                        e = getattr(self, cls)[r[_id]]
                        members.append(e)
                        memb_ids[e.key()] = {'id': _id, 'type': cls}
            return members, memb_ids

        for rid, reac in iteritems(getattr(self.parser, rclass[0])):

            left, l_ids = get_side(reac['left'])
            right, r_ids = get_side(reac['right'])
            left_attrs = {self.source: {rid: l_ids}}
            right_attrs = {self.source: {rid: r_ids}}

            nleft = \
                reduce(
                    lambda m1, m2:
                        m1 * m2,
                    map(
                        lambda m:
                            len(m.members) if hasattr(m, 'members') else 1,
                        left
                    ),
                    1
                )
            nright = \
                reduce(
                    lambda m1, m2:
                        m1 * m2,
                    map(
                        lambda m:
                            len(m.members) if hasattr(m, 'members') else 1,
                        right
                    ),
                    1
                )

            if len(left) or len(right):
                if nleft <= self.max_reaction_combinations and \
                        nright <= self.max_reaction_combinations:

                    reaction = Reaction(
                        left,
                        right,
                        left_attrs,
                        right_attrs,
                        source=self.source,
                        parent=self)
                    reaction.attrs[self.source][rid] = {}

                    this_refs = \
                        set(
                            list(
                                map(
                                    lambda r:
                                        self.rrefs[self.source][r],
                                    filter(
                                        lambda r:
                                            r in self.parser.pubrefs,
                                        reac['refs']
                                    )
                                )
                            )
                        )

                    reaction.attrs[self.source][rid]['refs'] = this_refs
                    reaction.attrs[self.source][rid]['type'] = rclass[1]

                    key = reaction.__str__()

                    if key in self.reactions:
                        # print(key, type(self.reactions[key]), self.reactions[key].__str__(), type(reaction), reaction.__str__())
                        self.reactions[key] += reaction
                    else:
                        self.reactions[key] = reaction

                    setattr(self, '%s_added' % rclass[0],
                            getattr(self, '%s_added' % rclass[0]) + 1)

                    self.rreactions[self.source][rid] = key

                else:
                    self.huge_reactions[self.source][rid] = max(nleft, nright)

    def merge_controls(self):
        self.controls_added = 0
        self._merge_controls(('controls', 'control'))

    def merge_catalyses(self):
        self.catalyses_added = 0
        self._merge_controls(('catalyses', 'catalysis'))

    def _merge_controls(self, cclass):

        if self.source not in self.rcontrols:
            self.rcontrols[self.source] = {}

        def get_party(_id):
            for cls in ['proteins', 'pfamilies', 'cvariations', 'reactions']:
                if _id in getattr(self, 'r%s' % cls)[self.source]:
                    key = getattr(self, 'r%s' % cls)[self.source][_id]
                    entity = getattr(self, cls)[key]
                    return (cls, key, entity)
            return None, None, None

        for cid, ctrl in iteritems(getattr(self.parser, cclass[0])):

            erclass, erkey, erent = get_party(ctrl['controller'])
            edclass, edkey, edent = get_party(ctrl['controlled'])

            # print('n = %u, erclass: %s, edclass: %s, er: %s, ed: %s' % (n, erclass, edclass, ctrl['controller'], ctrl['controlled']))

            if erent is not None and edent is not None:

                this_refs = \
                    set(
                        list(
                            map(
                                lambda r:
                                    self.rrefs[self.source][r],
                                filter(
                                    lambda r:
                                        r in self.parser.pubrefs,
                                    ctrl['refs']
                                )
                            )
                        )
                    ) \
                    if 'refs' in ctrl else set([])

                control = Control(
                    erent, edent, source=self.source, parent=self)
                control.attrs[self.source][cid] = {}
                control.attrs[self.source][cid]['refs'] = this_refs
                control.attrs[self.source][cid]['class'] = cclass[1]
                control.attrs[self.source][cid]['type'] = ctrl['type']

                key = control.__str__()

                if key in self.controls:
                    self.controls[key] += control
                else:
                    self.controls[key] = control
                setattr(self, '%s_added' % cclass[0],
                        getattr(self, '%s_added' % cclass[0]) + 1)

                self.rcontrols[self.source][cid] = key

    def basic_stats(self, exclude_empty=False):
        self.stats = {
            'proteins': {},
            'complexes': {},
            'mods': {},
            'reactions': {},
            'controls': {},
            'refs': {}
        }
        comb = []
        for n in xrange(1, len(self.sources) + 1):
            comb.extend(list(itertools.combinations(self.sources, n)))
        comb = \
            list(
                map(
                    lambda s:
                        (tuple(sorted(s)), set(s)),
                    comb
                )
            )
        for etyp in self.stats.keys():
            self.stats[etyp] = \
                dict(
                map(
                    lambda s:
                        (s[0], 0),
                    comb
                )
            )
            for e in getattr(self, etyp).values():
                for c in comb:
                    _sources = (
                        e.evidences.get_resource_names()
                            if hasattr(e, 'evidences') else
                        common.to_set(e.sources)
                    )
                    if c[1] <= _sources:
                        if \
                            not exclude_empty \
                                or (
                                    etyp not in
                                    ['complexes', 'reactions', 'controls']
                                ) or (
                                    etyp == 'complexes'
                                    and
                                    len(e.members) > 1
                                ) or (
                                    etyp == 'reactions'
                                    and
                                    len(e.left.members)
                                    and
                                    len(e.right.members)
                                ) or (
                                    etyp == 'controls'
                                    and (
                                        (
                                            (
                                                e.controller.__class__.__name__ == 'Complex'
                                                or
                                                e.controller.__class__.__name__ == 'ProteinFamily'
                                            ) and
                                            len(e.controller.members)
                                        ) or (
                                            e.controller.__class__.__name__ ==
                                'ComplexVariations'
                                            and
                                            any(map(lambda m: bool(len(m.members)),
                                                    e.controller.members))
                                        )
                                    ) and (
                                        len(e.controlled.left.members)
                                        and
                                        len(e.controlled.right.members)
                                    )
                                ):

                            self.stats[etyp][c[0]] += 1

    def simpson_stats(self):
        if not hasattr(self, 'stats'):
            self.basic_stats()
        self.simpson_sim = {
            'proteins': {},
            'complexes': {},
            'mods': {},
            'reactions': {},
            'controls': {},
            'refs': {}
        }
        for etyp in self.simpson_sim.keys():
            for s1 in self.sources:
                for s2 in self.sources:
                    if s1 != s2:
                        self.simpson_sim[etyp][(s1, s2)] = \
                            common.simpson_index_counts(
                                self.stats[etyp][tuple([s1])],
                                self.stats[etyp][tuple([s2])],
                                self.stats[etyp][tuple(sorted([s1, s2]))]
                        )

    def resource_graph_edges(self, etyp):
        if not hasattr(self, 'simpson_sim'):
            self.simpson_stats()
        stats = self.stats[etyp]
        sim = self.simpson_sim[etyp]
        edges = []
        nodes = {}
        for s1 in self.sources:
            nodes[s1] = stats[(s1, )]
            for s2 in self.sources:
                if s1 != s2 and sim[(s1, s2)] > 0.0:
                    edges.append([s1, s2, sim[(s1, s2)]])
        return edges, nodes

    def iterate_reactions(self):
        pass

    def load_sequences(self):
        if self.seq is None:
            self.seq = seq.swissprot_seq(
                self.ncbi_tax_id, isoforms=True)

    # interaction iterators from here

    def expand(self):

        def add_interactions(gen):
            for i in gen:
                key = (i[0], i[1])
                if key not in aggregate:
                    aggregate[key] = i
                    aggregate[key][2] = set([i[2]])
                else:
                    aggregate[key][4].update(i[4])
                    aggregate[key][5].update(i[5])
                    aggregate[key][2].add(i[2])
                    aggregate[key][3] = aggregate[key][3] or i[3]

        aggregate = {}
        add_interactions(self.in_same_component())
        add_interactions(self.co_control())
        add_interactions(self.interacts_with())
        add_interactions(self.state_change())
        self.interactions = list(aggregate.values())

    def expand_by_source(self):

        def add_interactions(gen):
            for i in gen:
                key = (i[0], i[1], i[4])
                if key not in aggregate:
                    aggregate[key] = i
                    aggregate[key][2] = set([i[2]])
                else:
                    aggregate[key][5].update(i[5])
                    aggregate[key][2].add(i[2])
                    aggregate[key][3] = aggregate[key][3] or i[3]

        aggregate = {}
        add_interactions(self.in_same_component(by_source=True))
        add_interactions(self.co_control(by_source=True))
        add_interactions(self.interacts_with(by_source=True))
        add_interactions(self.state_change(by_source=True))
        self.interactions_by_source = list(aggregate.values())
        pickle.dump(
            self.interactions_by_source,
            open(
                os.path.join(
                    self.cachedir,
                    'reaction_interactions_by_source.pickle',
                ),
                'wb'
            )
        )

    def in_same_component(self, by_source=False):
        """
        For all complexes connects all members of the complex with each other.
        """
        self.prg = progress.Progress(
            len(self.complexes), 'Expanding `in same component` interactions',
            1)
        aggregate_src = {}
        for c in self.complexes.values():
            self.prg.step()
            for i, p1 in enumerate(c):
                for p2 in list(c)[i + 1:]:
                    key = (p1, p2)
                    if key not in aggregate_src:
                        aggregate_src[key] = set([])
                    aggregate_src[key].update(c.sources)
        self.prg.terminate()
        for (p1, p2), s in iteritems(aggregate_src):
            if by_source:
                for ss in s:
                    yield [p1, p2, 'IN_SAME_COMPONENT', False, ss, set([])]
            else:
                yield [p1, p2, 'IN_SAME_COMPONENT', False, s, set([])]

    def protein_get_refs(self, source, protein_elem):
        if protein_elem in self.rproteins[source]:
            protein = self.proteins[self.rproteins[source][protein_elem]]
            elem = protein.attrs[source]['pids'][protein_elem]
            if 'refs' in elem:
                return elem['refs']
        return set([])

    def complex_get_refs(self, source, cplex_elem):
        refs = set([])
        for cplex_key in self.rcomplexes[source][cplex_elem]:
            cplex = self.complexes[cplex_key]
            if cplex_elem in cplex.attrs[source]:
                elem = cplex.attrs[source][cplex_elem]
                for protein, pdata in iteritems(elem):
                    if protein != 'children':
                        refs.update(
                            self.protein_get_refs(source, elem[protein][
                                'pid']))
        return refs

    def co_control(self, by_source=False):
        self.prg = progress.Progress(
            len(self.controls), 'Expanding `co-control` interactions', 1)
        aggregate_src = {}
        aggregate_ref = {}
        for co in self.controls.values():
            self.prg.step()
            if co.controller.type != 'pfamily':
                proteins = sorted(list(co.controller.proteins()))
                for i, p1 in enumerate(proteins):
                    for p2 in proteins[i + 1:]:
                        key = (p1, p2)
                        if key not in aggregate_src:
                            aggregate_src[key] = set([])
                        aggregate_src[key].update(co.sources)
                        for s, codata in iteritems(co.attrs):
                            for coid, d in iteritems(codata):
                                if by_source:
                                    if s not in aggregate_ref:
                                        aggregate_ref[s] = {}
                                    if key not in aggregate_ref[s]:
                                        aggregate_ref[s][key] = set([])
                                    aggregate_ref[s][key].update(d['refs'])
                                else:
                                    if key not in aggregate_ref:
                                        aggregate_ref[key] = set([])
                                    aggregate_ref[key].update(d['refs'])
        self.prg.terminate()
        for (p1, p2), s in iteritems(aggregate_src):
            if by_source:
                for ss in s:
                    yield [
                        p1, p2, 'CO_CONTROL', False, ss, aggregate_ref[ss][(
                            p1, p2)]
                    ]
            else:
                yield [p1, p2, 'CO_CONTROL', False, s, aggregate_ref[(p1, p2)]]

    def interacts_with(self, by_source=False):
        aggregate_src = {}
        aggregate_ref = {}
        self.prg = progress.Progress(
            len(self.reactions), 'Expanding `interacts with` interactions', 1)
        for rs in self.reactions.values():
            self.prg.step()
            isrc = \
                set(list(
                    map(
                        lambda s:
                            s[0],
                        filter(
                            lambda s:
                                len(
                                    list(
                                        filter(
                                            lambda rr:
                                                rr['type'] == 'interaction',
                                            s[1].values()
                                        )
                                    )
                                ),
                            iteritems(rs.attrs)
                        )
                    )
                ))
            if len(isrc):
                for r in rs.expand():
                    for i, p1 in enumerate(r.left.proteins()):
                        for p2 in list(r.right.proteins())[i + 1:]:
                            key = tuple(sorted([p1, p2]))
                            if key not in aggregate_src:
                                aggregate_src[key] = set([])
                            aggregate_src[key].update(isrc)
                            for s in isrc:
                                if by_source:
                                    if s not in aggregate_ref:
                                        aggregate_ref[s] = {}
                                    if key not in aggregate_ref[s]:
                                        aggregate_ref[s][key] = set([])
                                    for rid, rdata in iteritems(r.attrs[s]):
                                        if rdata['type'] == 'interaction':
                                            aggregate_ref[s][key].update(rdata[
                                                'refs'])
                                else:
                                    if key not in aggregate_ref:
                                        aggregate_ref[key] = set([])
                                    for rid, rdata in iteritems(r.attrs[s]):
                                        if rdata['type'] == 'interaction':
                                            aggregate_ref[key].update(rdata[
                                                'refs'])
        self.prg.terminate()
        for (p1, p2), s in iteritems(aggregate_src):
            if by_source:
                for ss in s:
                    yield [
                        p1, p2, 'INTERACTS_WITH', False, ss, aggregate_ref[ss][
                            (p1, p2)]
                    ]
            else:
                yield [
                    p1, p2, 'INTERACTS_WITH', False, s, aggregate_ref[(p1, p2)]
                ]

    def state_change(self, by_source=False):
        self.prg = progress.Progress(
            len(self.controls), 'Expanding `state change` interactions', 1)
        aggregate_src = {}
        aggregate_ref = {}
        for cos in self.controls.values():
            self.prg.step()
            for co in cos.expand():
                er_proteins = co.controller.proteins()
                for s in co.sources:
                    ldiff, rdiff = self.reaction_mod_diff(s, co.controlled)
                    for p2 in ldiff.keys():
                        if ldiff[p2] != rdiff[p2]:
                            for p1 in er_proteins:
                                key = (p1, p2)
                                if key not in aggregate_src:
                                    aggregate_src[key] = set([])
                                aggregate_src[key].add(s)
                                refs = set([])
                                for rid, rdata in iteritems(
                                        co.controlled.attrs[s]):
                                    refs.update(rdata['refs'])
                                if by_source:
                                    if s not in aggregate_ref:
                                        aggregate_ref[s] = {}
                                    if key not in aggregate_ref[s]:
                                        aggregate_ref[s][key] = set([])
                                    aggregate_ref[s][key].update(refs)
                                else:
                                    if key not in aggregate_ref:
                                        aggregate_ref[key] = set([])
                                    aggregate_ref[key].update(refs)
        self.prg.terminate()
        for (p1, p2), s in iteritems(aggregate_src):
            if by_source:
                for ss in s:
                    yield [
                        p1, p2, 'STATE_CHANGE', True, ss, aggregate_ref[ss][(
                            p1, p2)]
                    ]
            else:
                yield [
                    p1, p2, 'STATE_CHANGE', True, s, aggregate_ref[(p1, p2)]
                ]

    def reaction_mod_diff(self, source, reaction, by_rid=False):
        left = self.reaction_side_get_mods(source, reaction.left, by_rid)
        right = self.reaction_side_get_mods(source, reaction.right, by_rid)
        return common.dict_diff(left, right)

    def reaction_side_get_mods(self, source, rside, by_rid=False):
        mods = {}
        for rid, rd in iteritems(rside.attrs[source]):
            for ent_key, data in iteritems(rd):
                if data['type'] == 'proteins':
                    next_mods = {
                        ent_key: self.protein_get_mods(source, data['id'])
                    }
                elif data['type'] == 'complexes':
                    next_mods = self.complex_get_mods(source, data['id'])
                if by_rid:
                    if rid not in mods:
                        mods[rid] = {}
                    common.merge_dicts(mods[rid], {rid: next_mods})
                else:
                    common.merge_dicts(mods, next_mods)
        return mods

    def protein_get_mods(self, source, protein_elem):
        if protein_elem in self.rproteins[source]:
            protein = self.proteins[self.rproteins[source][protein_elem]]
            elem = protein.attrs[source]['pids'][protein_elem]
            if 'mods' in elem:
                return elem['mods']
        return set([])

    def complex_get_mods(self, source, complex_elem):
        mods = {}
        for cplex_key in self.rcomplexes[source][complex_elem]:
            cplex = self.complexes[cplex_key]
            if complex_elem in cplex.attrs[source]:
                elem = cplex.attrs[source][complex_elem]
                for protein, pdata in iteritems(elem):
                    if protein != 'children':
                        mods[protein] = self.protein_get_mods(
                            source, elem[protein]['pid'])
        return mods

# ## ## ## ## ## ##


class ReactionSide(AttributeHandler):
    def __init__(self, members, source=[], parent=None):
        super(ReactionSide, self).__init__()

        self.members = sorted(members)
        self.sources = set([])
        self.attrs = {}
        self.add_source(source)
        self.parent = parent
        self.is_expanded = False

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        map(lambda m: m.reload(), self.members)

    def __hash__(self):
        return hash(self.__str__())

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return 'ReactionSide: (%s)' % \
            ('+'.join(map(lambda m: m.__str__(), self.members)))

    def equality(self, one, two):
        return \
            all(
                map(
                    lambda m1:
                    any(
                        map(
                            lambda m2:
                                not m1.set.isdisjoint(m2.set),
                            filter(
                                lambda m2:
                                    m2.type == m1.type,
                                two.members
                            )
                        )
                    ),
                    filter(
                        lambda m1:
                        m1.type == 'cvariations' or m1.type == 'pfamily',
                        one.members
                    )
                )
            ) and \
            all(
                map(
                    lambda m1:
                    m1 in two.set,
                    filter(
                        lambda m1:
                        m1.type == 'protein',
                        one.members
                    )
                )
            )

    def __eq__(self, other):
        return self.equality(self, other) and self.equality(other, self)

    def __iter__(self):
        for m in self.members:
            yield m

    def expand(self):
        """
        Expands the ``ReactionSide`` by iterating over all combinations
        of all ``ComplexVariation`` and ``ProteinFamily`` members, so
        yields ``ReactionSide`` objects with only ``Protein`` and ``Complex``
        members. Yields tuple, because ``ReactionSide`` is initialized in
        ``Reaction``, the tuple is suitable to serve as ``members`` and
        ``attrs``.
        """
        # collecting protein attributes
        if self.is_expanded:
            for i in [1]:
                yield self.members, self.attrs
        else:
            try:
                pattrs = \
                    dict(
                        map(
                            lambda m:
                                (
                                    m.id,
                                    dict(
                                        map(
                                            lambda d1:
                                                (
                                                    d1[0],
                                                    dict(
                                                        map(
                                                            lambda d2:
                                                                (d2[0], d2[
                                                                 1][m.id]),
                                                            iteritems(d1[1])
                                                        )
                                                    )
                                                ),
                                            iteritems(self.attrs)
                                        )
                                    )
                                ),
                            filter(
                                lambda m:
                                    m.type == 'protein',
                                self.members
                            )
                        )
                    )
            except:
                print(self.attrs)
            for c in \
                    itertools.product(
                        *list(
                            map(
                                lambda m:
                                    list(
                                        zip(
                                            m.itermembers(),
                                            [m.key()] * (
                                                len(m.members)
                                                if hasattr(m, 'members')
                                                else 1)
                                        )
                                    ),
                                self.members
                            )
                        )
                    ):

                attrs = dict(map(lambda s: (s,
                                            dict(map(lambda rid: (rid, {}),
                                                     self.attrs[s].keys()))), self.sources))
                members = []
                for ((m, a), k) in c:
                    members.append(m)
                    if m.type == 'protein':
                        # if it was a protein, we just copy
                        if m.id in pattrs:
                            for s, d1 in iteritems(pattrs[m.id]):
                                for rid, d2 in iteritems(d1):
                                    attrs[s][rid][m.id] = d2
                        # if it is from a protein family
                        else:
                            # for each resource
                            for s, r in iteritems(attrs):
                                # for each original reaction id
                                for rid, d in iteritems(self.attrs[s]):
                                    # the key of the new entity (here: str,
                                    # uniprot id)
                                    if k in self.attrs[s][rid]:
                                        attrs[s][rid][m.key()] = (
                                            # the type is obvious, the id is from the `a` dict supplied
                                            # by the ProteinFamily object, and we look up the id belonging
                                            # to the key of the original entity
                                            {
                                                'type': 'proteins',
                                                'id': a[s][self.attrs[s][rid][
                                                    k]['id']]
                                            })
                    # if it is a complex from a complex variations
                    elif m.type == 'complex':
                        for s, r in iteritems(attrs):
                            for rid, d in iteritems(r):
                                if k in self.attrs[s][rid]:
                                    cid = self.attrs[s][rid][k]['id']
                                    attrs[s][rid][m.key()] = \
                                        {'type': 'complexes', 'id': cid}
                yield members, attrs

    def proteins(self):
        return \
            reduce(
                lambda m1, m2:
                    m1 | m2.proteins(),
                self.members,
                set([])
            )

    def __ror__(self, other):
        return self.proteins() | other

    def __rand__(self, other):
        return self.proteins() & other

    def __rsub__(self, other):
        return self.proteins() - other


class Reaction(AttributeHandler):
    def __init__(self,
                 left,
                 right,
                 left_attrs,
                 right_attrs,
                 source=[],
                 parent=None):
        super(Reaction, self).__init__()

        self.parent = parent
        self.left = ReactionSide(left, source, parent=self.parent)
        self.right = ReactionSide(right, source, parent=self.parent)
        self.left.merge_attrs(left_attrs)
        self.right.merge_attrs(right_attrs)
        self.attrs = {}
        self.sources = set([])
        self.add_source(source)
        self.is_expanded = False

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        self.left.reload()
        self.right.reload()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return 'Reaction: LEFT(%s) --> RIGHT(%s)' % \
            (self.left.__str__(), self.right.__str__())

    def __hash__(self):
        return hash(self.__str__())

    def __eq__(self, other):
        return self.left == other.left and self.right == other.right

    def __iadd__(self, other):
        self = AttributeHandler.__iadd__(self, other)
        self.left += other.left
        self.right += other.right
        return self

    def expand(self):
        if self.is_expanded:
            for i in [1]:
                yield self
        else:
            expanded = []

            lAllProteins = self.left.proteins()
            rAllProteins = self.right.proteins()
            diffAllProteins = rAllProteins ^ lAllProteins

            lefts = list(self.left.expand())
            rights = list(self.right.expand())

            diffs = set([])
            for t in xrange(2):
                if t == 1:
                    try:
                        minDiff = min(diffs)
                    except ValueError:
                        print('Empty sequence error: %s' % self.__str__())
                for left in lefts:
                    for right in rights:
                        r = Reaction(
                            left[0],
                            right[0],
                            left[1],
                            right[1],
                            source=self.sources,
                            parent=self.parent)
                        lProteins = r.left.proteins()
                        rProteins = r.right.proteins()
                        diffProteins = lProteins ^ rProteins

                        diff = len(diffProteins)
                        diffs.add(diff)
                        if t == 1 and diff == minDiff:
                            r.left.is_expanded = True
                            r.right.is_expanded = True
                            r.is_expanded = True
                            r.attrs = self.attrs

                            yield r

                        else:
                            del r
            #
            #minDiff = min(map(lambda e: e[1], expanded))

            # for r, d in expanded:
            # if d == minDiff:
            # yield r


class Control(AttributeHandler):
    def __init__(self, er, ed, source=[], parent=None):
        super(Control, self).__init__()

        self.controller = er
        self.controlled = ed
        self.attrs = {}
        self.sources = set([])
        self.add_source(source)
        self.parent = parent
        self.is_expanded = False

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def __str__(self):
        return 'Control: C.ER(%s) --> C.ED(%s)' % \
            (
                self.controller.__str__(),
                self.controlled.__str__()
            )

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.__str__())

    def __eq__(self, other):
        return self.controller == other.controller \
            and self.controlled == other.controlled

    def expand(self):
        if self.is_expanded:
            for i in [self]:
                yield self
        else:
            for ed in self.controlled.expand():
                for er, erattrs in self.controller.itermembers():
                    c = Control(
                        er, ed, source=self.sources, parent=self.parent)
                    c.attrs = self.attrs
                    yield c
