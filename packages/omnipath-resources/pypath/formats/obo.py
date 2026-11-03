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

import re
import importlib as imp
import collections
import itertools

import pypath.share.session as session
from pypath.share.downloads import download_and_open



class OboValue(
        collections.namedtuple(
            'OboValueBase',
            [
                'value',
                'modifiers',
                'comment',
            ],
        )
    ):
    
    def __new__(cls, value, modifiers = None, comment = None):
        
        value = value.strip('"')
        
        return super(OboValue, cls).__new__(cls, value, modifiers, comment)


class Obo(session.Logger):
    """
    Reader for an OBO (Open Biomedical Ontologies) file.
    See more details at
    https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html#S.1
    
    :param str obofile:
        Path or URL to an OBO file.
    :param list single_tags:
        A list of tag names which have a single occurence in each stanza and
        ideally should present in all stanzas.
    """
    
    _single_tags = [
        'stanza',
        'id',
        'name',
        'definition',
        'namespace',
    ]
    restanza = re.compile(r'\[(\w+)\]\s*')
    retag = re.compile(
        r'([\w_]+):\s?' # tag
        r'([^\s^"]+|".*")\s?' # value
        r'(?:([^!]+[^!^\s]))?' # modifiers
        r'(?:'
            r'\s*!?\s*' # comment separator
            r'([^!]+[^\s])' # comment
        r'?)\s*'
    )
    _disallowed_keys = {
        'def': 'definition',
    }
    
    def __init__(
            self,
            obofile,
            single_tags = None,
            name = None,
        ):
        
        session.Logger.__init__(self, name = 'obo')
        
        self.obofile = obofile
        self.single_tags = single_tags or self._single_tags
        self._single_tags_set = set(self.single_tags)
        self.name = name or 'OBO'
        self._set_record()
        self.open()
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def open(self):

        # Download and open file using download manager if it's a URL
        if self.obofile.startswith('http://') or self.obofile.startswith('https://'):
            self.reader = download_and_open(
                self.obofile,
                filename=self.obofile.split('/')[-1],
                subfolder='obo',
                large=True,
                encoding='utf-8',
            )
        else:
            # For local files, create a simple object with .result attribute
            # to match the Opener interface
            class LocalFileOpener:
                def __init__(self, filepath):
                    self.result = open(filepath, 'r')
            self.reader = LocalFileOpener(self.obofile)

        self._in_stanza = False

        self._log('OBO reader created, opened `%s`.' % self.obofile)
    
    
    def _set_record(self):
        
        fields = self.single_tags + ['attrs',]
        
        self.record = collections.namedtuple(
            '%sRecord' % self.name.capitalize(),
            fields,
        )
        self.record.__new__.__defaults__ = (None,) * len(fields)
    
    
    def __iter__(self):
        
        if hasattr(self.reader.result, 'seek'):
            
            self.reader.result.seek(0)
        
        self._reset_record()
        
        iterlines0, iterlines1 = itertools.tee(self.reader.result)
        
        for line, next_line in itertools.zip_longest(
            iterlines0,
            itertools.islice(iterlines1, 1, None),
        ):
            
            if self._in_stanza:
                
                if not next_line or next_line[0] == '[':
                    
                    yield self._create_record()
                    self._reset_record()
                    self._in_stanza = False
                    
                else:
                    
                    m = self.retag.match(line)

                    if m:

                        tag, *value = m.groups()
                        tag = self._disallowed_keys.get(tag) or tag

                        # Special handling for 'name' field: concatenate value and modifiers
                        # because name fields in OBO files are multi-word unquoted strings
                        # without actual modifiers
                        if tag == 'name' and value[1]:  # value[1] is modifiers
                            value = (f"{value[0]} {value[1]}", None, value[2])

                        value = OboValue(*value)

                        if tag in self._single_tags_set:

                            self._current_record[tag] = value

                        else:

                            self._current_record['attrs'][tag].add(value)
                
            elif line[0] == '[':
                
                m = self.restanza.match(line)
                self._current_record['stanza'] = m.groups()[0]
                self._in_stanza = True
    
    
    def _create_record(self):
        
        self._current_record['attrs'] = dict(self._current_record['attrs'])
        
        return self.record(**self._current_record)
    
    
    def _reset_record(self):
        
        self._current_record = {
            'attrs': collections.defaultdict(set)
        }
    
    
    def __repr__(self):
        
        return '<OBO file `%s`>' % self.obofile


    def parent_terms(self):
        """
        Retrieves is_a relations between ontology terms.
        Defines a dictionary: keys are the ontology term IDs and 
        values are list of their parent terms.
        """

        self.parents = dict(
                (term.id.value, 
                [is_a.value for is_a in term.attrs["is_a"]]
                        if "is_a" in term.attrs else 
                []
            )
                for term in self
                if term.stanza == "Term"
            )
