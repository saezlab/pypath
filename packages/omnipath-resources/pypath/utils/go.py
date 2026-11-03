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

from future.utils import iteritems
from past.builtins import xrange, range

import os
import sys
import importlib as imp
import re
from collections import Counter, OrderedDict
import numpy as np
import itertools

try:
    import cPickle as pickle
except:
    import pickle

import pypath.share.cache as cache
import pypath.inputs.go as go_input
import pypath.share.progress as progress
import pypath.share.common as common
from pypath.share.common import *
import pypath.share.session as session_mod
import pypath.share.settings as settings

# this is for GO terms parsing:
_reexprterm = re.compile(r'and|or|not|\(|\)|GO:[0-9]{7}')
_reexprname = re.compile(
    r'(?!\s)' # no space at the beginning
    r'(?:AND|OR|NOT|\(|\)|' # either AND, OR, NOT or parentheses
       r'(?:(?!OR|AND|NOT|\s{2:})(?:[-\w: ]))+)' # or something else
                                                 # (words with spaces)
    r'(?<!\s)' # no space at the end
)


ROOT_NODES = {
    'cellular_component': 'GO:0005575',
    'biological_process': 'GO:0008150',
    'molecular_function': 'GO:0003674',
}

ROOT_ACS = set(ROOT_NODES.values())


class GeneOntology(session_mod.Logger):

    all_relations = {
        'is_a',
        'part_of',
        'occurs_in',
        'regulates',
        'positively_regulates',
        'negatively_regulates',
    }


    def __init__(
            self,
            terms = None,
            ancestors = None,
            descendants = None,
            aspect = None,
            term = None,
            name = None,
        ):
        """
        Loads data about Gene Ontology terms and their relations.
        """

        session_mod.Logger.__init__(self, name = 'go')

        self._terms_provided = terms
        self._ancestors_provided = ancestors
        self._descendants_provided = descendants
        self._aspect_provided = aspect
        self._term_provided = term
        self._name_provided = name

        self._load()


    def reload(self):
        """Reloads the object from the module level."""

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def _load(self):

        self._log('Populating Gene Ontology: ontology.')

        self._load_terms()
        self._load_tree()
        self._set_aspect()
        self._set_name()
        self._set_term()

        # delattr(self, '_terms')

        self._log('Gene Ontology: ontology populated.')


    def _load_terms(self):

        self._terms = self._terms_provided or go_input.go_terms_quickgo()


    def _load_tree(self):

        self._log('Gene Ontology: building the ontology tree.')

        self.ancestors = (
            self._ancestors_provided or
            self._merge_aspects(
                go_input.go_ancestors_quickgo()
            )
        )
        self.descendants = (
            self._descendants_provided or
            self._merge_aspects(
                go_input.go_descendants_quickgo()
            )
        )


    def _set_aspect(self):

        self.aspect = (
            self._aspect_provided or
            dict(
                (term, asp)
                for asp, terms in iteritems(self._terms)
                for term in terms.keys()
            )
        )


    def _set_name(self):

        self._log('Collecting short names of GO terms.')

        self.name = (
            self._name_provided or
            dict(
                i
                for ii in self._terms.values()
                for i in iteritems(ii)
            )
        )


    def _set_term(self):

        self.term = (
            self._term_provided or
            dict(
                reversed(i)
                for i in iteritems(self.name)
            )
        )


    def is_term(self, term):
        """
        Tells if ``term`` is a GO accession number.
        """

        return term in self.name


    def is_name(self, name):
        """
        Tells if ``name`` is a GO term name.
        """

        return name in self.term


    def get_name(self, term):
        """
        For a GO accession number returns the name of the term.
        If ``term`` is already a GO term name returns it unchanged.
        """

        return (
            term
                if self.is_name(term) else
            None
                if term not in self.name else
            self.name[term]
        )


    def get_term(self, name):
        """
        For a GO term name returns its GO accession number.
        If ``name`` is a GO accession returns it unchanged.
        """

        result = (
            name
                if self.is_term(name) else
            None
                if name not in self.term else
            self.term[name]
        )

        if result is None:

            self._log('Could not find GO term name: `%s`.' % name)

        return result


    def terms_to_names(self, terms):
        """
        For a list of GO names returns a list of tuples with the terms
        and their names.
        """

        return [(term, self.get_name(term)) for term in terms]


    def terms_to_names_aspects(self, terms):
        """
        For a list of GO terms returns a list of tuples with the terms,
        their names and the ontology aspect.
        """

        return [
            (term, self.get_name(term), self.get_aspect(term))
            for term in terms
        ]


    def names_to_terms(self, names):
        """
        For a list of GO terms returns a list of tuples with the terms
        and their names.
        """

        return [(self.get_term(name), name) for name in names]


    def names_to_terms_aspects(self, names):
        """
        For a list of GO namess returns a list of tuples with the terms,
        their names and ontology aspects.
        """

        return [
            (self.get_term(name), name, self.aspect_from_name(name))
            for name in names
        ]


    def aspect_from_name(self, name):
        """
        Tells about a Gene Ontology term name which aspect does it belong to.
        """

        term = self.get_term(name)

        if term:

            return self.get_aspect(term)


    @staticmethod
    def _merge_aspects(dct):

        dct['P'].update(dct['C'])
        dct['P'].update(dct['F'])
        return dct['P']


    def subgraph_nodes(
            self,
            direction,
            terms,
            relations = None,
            include_seed = True,
        ):
        """
        Returns a set of all nodes either in the subgraph of ancestors
        or descendants of a single term or a set of terms.

        :param str direction:
            Possible values: `ancestors` or `descendants`.
        :param bool include_seed:
            Include ``terms`` in the subgraph or only the related nodes.
        """

        relations = relations or self.all_relations

        if isinstance(terms, str):

            terms = {terms}

        graph = getattr(self, direction)
        subgraph = set(terms) if include_seed else set()

        for term in terms:

            if term not in graph:

                if term not in ROOT_ACS:

                    self._log(
                        'GO term without known %ss: `%s`.' % (direction, term)
                    )

                continue

            for related, relation in graph[term]:

                if relation not in relations:

                    continue

                if related not in subgraph:

                    subgraph.update(
                        self.subgraph_nodes(direction, related, relations)
                    )
                    subgraph.add(related)

        return subgraph


    def get_all_ancestors(self, terms, relations = None, include_seed = True):
        """
        Returns a set of all ancestors of a single term or a set of terms.
        """

        terms = self.set_of_terms(terms)

        return self.subgraph_nodes(
            direction = 'ancestors',
            terms = terms,
            relations = relations,
            include_seed = include_seed,
        )


    def get_all_descendants(
            self,
            terms,
            relations = None,
            include_seed = True,
        ):
        """
        Returns a set of all descendants of a single term or a set of terms.
        """

        terms = self.set_of_terms(terms)

        return self.subgraph_nodes(
            direction = 'descendants',
            terms = terms,
            relations = relations,
            include_seed = include_seed,
        )


    def get_aspect(self, term):
        """
        For a GO term tells which aspect does it belong to.
        Returns `None` if the term is not in the ontology.
        """

        if term in self.aspect:

            return self.aspect[term]


    def all_from_aspect(self, aspect):
        """
        Returns the set of all GO terms of one aspect.
        """

        return set(
            term
            for term, asp in iteritems(self.aspect)
            if asp == aspect
        )


    def is_root(self, term):
        """
        Tells if a term is the root of the graph i.e. it has no ancestors.
        """

        return term in self.ancestors and bool(self.ancestors[term])


    def is_leaf(self, term):
        """
        Tells if a term is a leaf of the graph i.e. it has no descendants.
        """

        return (
            (
                term in self.ancestors and
                term not in self.descendants
            ) or
            not bool(self.descendants[term])
        )


    def lowest(self, terms, *args):
        """
        From a set of terms returns the lowest level ones, removing all which
        are parents of some others in the set.
        """

        return self.flatten(terms, *args)


    def highest(self, terms, *args):
        """
        From a set of terms returns the highest level ones, removing all
        which are descendants of some others in the set.
        """

        return self.flatten(terms, *args, lowest = False)


    def flatten(self, terms, *args, lowest = True):
        """
        Returns a set of terms by removing either all redundant ancestors or
        descendants from the provided set terms. By removing the ancestors
        you get the lowest level set of terms, by removing the descendants
        the result will be the highest level non-redundant terms.

        :param str direction:
            Either `lowest` or `highest`.
        """

        terms = self.set_of_terms(terms, *args)

        method = getattr(
            self,
            'get_all_%s' % (
                'ancestors' if lowest else 'descendants'
            )
        )

        return (
            terms -
            set.union(
                *(
                    method(term, include_seed = False)
                    for term in terms
                )
            )
        )


    def set_of_terms(self, terms_names, *args):
        """
        Converts anything to a set of terms. ``terms_names`` can be either a
        single term or name or an iterable of terms and names.
        """

        return self.set_of(terms_names, *args)


    def set_of_names(self, terms_names, *args):
        """
        Converts anything to a set of names. ``terms_names`` can be either a
        single term or name or an iterable of terms and names.
        """

        return self.set_of(terms_names, *args, to_terms = False)


    def set_of(self, terms_names, *args, to_terms = True):
        """
        Converts anything to a set of terms or names.
        ``terms_names`` can be either a single term or name or an iterable
        of terms and names.

        :param bool to_terms:
            The target identifier type is `term`; if ``False`` the target
            will be `name`.
        """

        if isinstance(terms_names, str):

            terms_names = {terms_names}

        elif not isinstance(terms_names, set):

            terms_names = set(terms_names)

        if args:

            terms_names.update(set(args))

        method = getattr(
            self,
            'get_term' if to_terms else 'get_name'
        )

        return {
            method(term) for term in terms_names
        }


class GOAnnotation(session_mod.Logger):

    aspects = ('C', 'F', 'P')


    def __init__(
            self,
            organism = 9606,
            ontology = None,
            pickle_file = None,
            use_pickle_cache = True,
        ):
        """
        For one organism loads Gene Ontology annotations, in addition it
        accepts or creates a ``GeneOntology`` object.
        """

        session_mod.Logger.__init__(self, name = 'go')

        self.organism = organism
        self._pickle_file = pickle_file
        self._use_pickle_cache = use_pickle_cache

        if self._pickle_cache_load_hook():

            return

        self.ontology = ontology or GeneOntology()

        self._log(
            'Populating Gene Ontology: '
            'annotations for organism `%u`.' % organism
        )

        annot = go_input.go_annotations_goa(organism = organism)
        self.c = annot['C']
        self.f = annot['F']
        self.p = annot['P']

        self._ancestors_annotate()
        self._merge_annotations()

        self._pickle_cache_save_hook()


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def _ancestors_annotate(self):

        self._log('Creating ancestors lookup dictionary.')

        for asp in self.aspects:

            setattr(
                self,
                '%s_full' % asp.lower(),
                dict(
                    (
                        uniprot,
                        self.ontology.get_all_ancestors(annot)
                    )
                    for uniprot, annot in
                    iteritems(getattr(self, asp.lower()))
                )
            )


    def _merge_annotations(self):

        self._log('Creating complete lookup dictionary.')

        uniprots = self.all_uniprots()

        self.all = dict(
            (
                uniprot,
                set.union(*(
                    self.get_annot(uniprot, asp)
                    for asp in self.aspects
                ))
            )
            for uniprot in uniprots
        )

        self.all_full = dict(
            (
                uniprot,
                set.union(*(
                    self.get_annot_ancestors(uniprot, asp)
                    for asp in self.aspects
                ))
            )
            for uniprot in uniprots
        )


    def get_name(self, term):
        """
        For a GO accession number returns the name of the term.
        """

        return self.ontology.get_name(term)


    def get_term(self, name):
        """
        For a GO term name returns its GO accession number.
        """

        return self.ontology.get_term(name)


    def get_annot(self, uniprot, aspect):
        """
        For a UniProt ID returns its direct annotations from one aspect
        of Gene Ontology.
        Returns set.
        """

        annot = getattr(self, aspect.lower())
        return annot[uniprot] if uniprot in annot else set()


    def get_annots(self, uniprot):
        """
        For a UniProt ID returns its direct annotations from all aspects
        of Gene Ontology.
        Returns set.
        """

        return self.all[uniprot] if uniprot in self.all else set()


    def get_annot_ancestors(self, uniprot, aspect):
        """
        For a UniProt ID returns its annotations including lowest level
        terms and their ancestors from one aspect of Gene Ontology.
        Returns set.
        """

        annot = getattr(self, '%s_full' % aspect.lower())
        return annot[uniprot] if uniprot in annot else set()


    def get_annots_ancestors(self, uniprot):
        """
        For a UniProt ID returns its annotations including lowest level
        terms and their ancestors from all aspects of Gene Ontology.
        Returns set.
        """

        return self.all_full[uniprot] if uniprot in self.all_full else set()


    def has_term(self, uniprot, term):
        """
        Tells if an UniProt ID is annotated with a GO term.
        """

        return uniprot in self.all_full and term in self.all_full[uniprot]


    def has_any_term(self, uniprot, terms):
        """
        Tells if an UniProt ID is annotated with any of a set of GO terms.
        """

        return uniprot in self.all_full and term & self.all_full[uniprot]


    def all_uniprots(self):
        """
        Returns all UniProt IDs having annotations.
        """

        return set.union(*(
            set(getattr(self, asp.lower()).keys())
            for asp in self.aspects
        ))


    def i_select_by_term(self, term, uniprots = None):
        """
        Accepts a list of UniProt IDs and one or more gene ontology terms
        and returns a set of indices of those UniProts which are annotated
        with any of the terms.
        If no UniProts given all annotated UniProts considered.

        :param str,set term:
            A single GO term or set of terms.
        """

        uniprots = uniprots or sorted(self.all_uniprots())

        method = self.has_any_term if isinstance(term, set) else self.has_term

        return set(
            i
            for i, uniprot in enumerate(uniprots)
            if method(uniprot, term)
        )


    def select_by_name(self, name, uniprots = None, return_uniprots = False):
        """
        Accepts a list of UniProt IDs and one or more gene ontology names
        and returns the UniProts which are annotated with any of the names.
        If no UniProts given all annotated UniProts returned.

        :param str,set name:
            A single GO term name or set of names.
        :param bool return_uniprots:
            By default returns list of indices; if ``True`` returns a set of
            the selected UniProt IDs.
        """

        if isinstance(name, str):

            term = self.ontology.get_term(name)

        else:

            term = set(i[0] for i in self.ontology.names_to_terms(name))

        return self.select(
            term,
            uniprots = uniprots,
            return_uniprots = return_uniprots,
        )


    def select_by_term(self, term, uniprots = None):
        """
        Accepts a list of UniProt IDs and one or more gene ontology terms
        and returns the UniProts which are annotated with any of the terms.
        If no UniProts given all annotated UniProts returned.

        :param str,set term:
            A single GO term or set of terms.
        """

        uniprots = uniprots or sorted(self.all_uniprots())

        return set(
            np.array(uniprots)[
                list(self.i_select_by_term(term, uniprots))
            ]
        )


    def expr_names_to_terms(self, expr):
        """
        Processes an expression built by names to expressions of terms.

        :arg str expr:
            An expression using Gene Ontology names, parentheses and logical
            operators.
        """

        not_name = {'(', ')', 'AND', 'OR', 'NOT'}

        tokens_names = _reexprname.findall(expr)
        tokens_terms = []

        if tokens_names:

            for t in tokens_names:

                t = t.strip()

                if not t:

                    continue

                tokens_terms.append((
                    t
                        if t[:3] == 'GO:' else
                    t.lower()
                        if t in not_name else
                    self.get_term(t)
                ))

        return tokens_terms


    def select_by_expr(
            self,
            expr,
            uniprots = None,
            return_uniprots = False,
        ):
        """
        Selects UniProts based on an expression of Gene Ontology terms.
        Operator precedence not considered, please use parentheses.
        Return indices of the selected elements in the ``uniprots`` list
        or the set of selected UniProt IDs.

        :param str expr:
            An expression of Gene Ontology terms and names. E.g.
            ``'(GO:0005576 and not GO:0070062) or GO:0005887'``. Parentheses
            and operators ``and``, ``or`` and ``not`` can be used.
            Another example:
            ``hormone binding AND (cell surface OR GO:0009897)``.
        :param bool return_uniprots:
            By default returns list of indices; if ``True`` returns a set of
            the selected UniProt IDs.
        """

        expr_terms = self.expr_names_to_terms(expr)

        result = self.select_by_expr_terms(
            expr = expr_terms,
            uniprots = uniprots,
            return_uniprots = return_uniprots,
        )

        if any(e is None for e in result):

            self._log(
                'Could not process Gene Ontology expression: `%s`. '
                'Please check if the expression consists only of '
                'GO terms and/or ACs, the operators `AND`, `OR` and '
                '`NOT`, and braces. Whitespaces and newlines are OK. '
                'If you think the expression is correct please open '
                'an issue for `pypath`.' % expr,
                -9,
            )

        return result


    def select_by_expr_terms(
            self,
            expr,
            uniprots = None,
            return_uniprots = False,
        ):
        """
        Selects UniProts based on an expression of Gene Ontology terms.
        Operator precedence not considered, please use parentheses.
        Return indices of the selected elements in the ``uniprots`` list
        or the set of selected UniProt IDs.

        :param str expr:
            An expression of Gene Ontology terms. E.g.
            ``'(GO:0005576 and not GO:0070062) or GO:0005887'``. Parentheses
            and operators ``and``, ``or`` and ``not`` can be used.
        :param bool return_uniprots:
            By default returns list of indices; if ``True`` returns a set of
            the selected UniProt IDs.
        """

        ops = {
            'and': 'intersection',
            'or':  'union',
        }

        # if no UniProts provided does not make sense to return indices
        return_uniprots = return_uniprots or uniprots is None

        uniprots = uniprots or sorted(self.all_uniprots())

        if isinstance(expr, str):

            # tokenizing expression if it is a string
            # (method is recursive)
            expr = _reexprterm.findall(expr)

            if any(e is None for e in expr):

                self._log(
                    'Could not process Gene Ontology expression: `%s`. '
                    'Please check if the expression consists only of '
                    'GO terms and/or ACs, the operators `AND`, `OR` and '
                    '`NOT`, and braces. Whitespaces and newlines are OK. '
                    'If you think the expression is correct please open '
                    'an issue for `pypath`.' % expr,
                    -9,
                )
                return set()

        # initial values
        result   = set()
        stack    = []
        sub      = False
        negate   = False
        op       = None
        this_set = None

        for it in expr:

            # processing expression by tokens

            # we are in a sub-selection part
            if sub:

                if it == ')':

                    # token is a closing parenthesis
                    # execute sub-selection
                    this_set = self.select_by_expr_terms(
                        expr = stack,
                        uniprots = uniprots,
                    )
                    # empty stack
                    stack = []
                    sub = False

                else:

                    # token is something else
                    # add to sub-selection stack
                    stack.append(it)

            elif it is None:

                self._log(
                    'One part of a Gene Ontology the expression failed to '
                    'translate to GO AC. Substituting with empty set, this '
                    'will alter your results. Check for more specific '
                    'information earlier in the log.'
                )
                this_set = set()

            # we do actual processing of the expression
            elif it.lower() == 'not':

                # token is negation
                # turn on negation for the next set
                negate = True
                continue

            # open a sub-selection part
            elif it == '(':

                # token is a parenthesis
                # start a new sub-selection
                sub = True
                continue

            elif it[:3] == 'GO:':

                # token is a GO term
                # get the vertex selection by the single term method
                this_set = self.i_select_by_term(it, uniprots = uniprots)

                if negate:

                    # take the inverse of the current set
                    this_set = set(xrange(len(uniprots))) - this_set
                    # set negation again to False
                    negate = False

            elif it.lower() in ops:

                # token is an operator
                # set it for use at the next operation
                op = ops[it.lower()]

            # we found a set
            if this_set is not None:

                # and an operator
                if op is not None:

                    result = getattr(result, op)(this_set)

                # this normally happens only at the first set
                else:

                    result = this_set

                this_set = None
                op       = None

        return self._uniprot_return(result, uniprots, return_uniprots)


    def select(self, terms, uniprots = None, return_uniprots = False):
        """
        Retrieves the UniProt IDs annotated with any Gene Ontology terms or
        their descendants, or evaluates string expression
        (see ``select_by_expr``).
        Returns indices of the selected elements in the ``uniprots`` list
        or the set of selected UniProt IDs.

        :param str,set terms:
            A single GO term, a set of GO terms or an expression with
            GO terms.
        :param bool return_uniprots:
            By default returns list of indices; if ``True`` returns a set of
            the selected UniProt IDs.
        """

        return_uniprots = return_uniprots or uniprots is None

        uniprots = uniprots or sorted(self.all_uniprots())

        # this is not an individual term but an expression
        if (
            isinstance(terms, str) and
            not terms.startswith('GO')
        ):

            result = self.select_by_expr(terms, uniprots = uniprots)

        # either one term or a set of terms
        else:

            result = self.i_select_by_term(terms, uniprots = uniprots)

        return self._uniprot_return(result, uniprots, return_uniprots)


    def select_by_all(self, terms, uniprots = None, return_uniprots = False):
        """
        Selects the nodes annotated by all GO terms in ``terms``.

        Returns indices of the selected elements in the ``uniprots`` list
        or the set of selected UniProt IDs.

        :param list terms:
            List, set or tuple of GO terms.
        :param bool return_uniprots:
            By default returns list of indices; if ``True`` returns a set of
            the selected UniProt IDs.
        """

        return_uniprots = return_uniprots or uniprots is None

        uniprots = uniprots or sorted(self.all_uniprots())

        idx = set.intersection(*[self.select_by_term(term) for term in terms])

        return self._uniprot_return(idx, uniprots, return_uniprots)


    def _uniprot_return(self, idx, uniprots, return_uniprots):

        if return_uniprots:

            return set(np.array(uniprots)[list(idx)])

        return idx


    def _pickle_cache_load_hook(self):

        if not self._use_pickle_cache:

            return

        self._set_pickle_path()
        self.load_from_pickle()

        return os.path.exists(self._pickle_file)


    def _pickle_cache_save_hook(self):

        self._set_pickle_path()
        self.save_to_pickle()


    def _set_pickle_path(self):

        self._pickle_file = (
            self._pickle_file or
            os.path.join(
                cache.get_cachedir(),
                settings.get('go_pickle_cache_fname') % self.organism,
            )
        )


    def save_to_pickle(self, pickle_file = None):

        pickle_file = pickle_file or self._pickle_file

        if not isinstance(pickle_file, str):

            self._log(
                'Pickle file path must be a string: `%s`.' % str(pickle_file)
            )
            return

        self._log('Saving to pickle `%s`.' % pickle_file)

        with open(pickle_file, 'wb') as fp:

            pickle.dump(
                obj = (
                    self.c_full,
                    self.p_full,
                    self.f_full,
                    self.all_full,
                    self.c,
                    self.p,
                    self.f,
                    self.all,
                    self.ontology._terms,
                    self.ontology.ancestors,
                    self.ontology.descendants,
                    self.ontology.term,
                    self.ontology.name,
                ),
                file = fp,
            )

        self._log('Saved to pickle `%s`.' % pickle_file)


    def load_from_pickle(self, pickle_file = None):

        pickle_file = pickle_file or self._pickle_file

        if not os.path.exists(pickle_file):

            self._log('Pickle file does not exist: `%s`.' % str(pickle_file))
            return

        self._log('Loading from pickle `%s`.' % pickle_file)

        with open(pickle_file, 'rb') as fp:

            (
                self.c_full,
                self.p_full,
                self.f_full,
                self.all_full,
                self.c,
                self.p,
                self.f,
                self.all,
                ontology_terms,
                ontology_ancestors,
                ontology_descendants,
                ontology_term,
                ontology_name,
            ) = pickle.load(fp)

        self.ontology = GeneOntology(
            terms = ontology_terms,
            ancestors = ontology_ancestors,
            descendants = ontology_descendants,
            term = ontology_term,
            name = ontology_name,
        )

        self._log('Loaded from pickle `%s`.' % pickle_file)


class GOCustomAnnotation(session_mod.Logger):


    def __init__(
            self,
            categories,
            go_annot = None,
            ncbi_tax_id = 9606,
        ):
        """
        Provides annotations by a custom set of GO terms or expressions
        built from multiple terms.

        :arg dict categories:
            A dict with custom category labels as keys and single GO terms
            or names or complex expressions as values.
            Alternatively a set of GO terms, in this case the term names will
            be used as labels.
        :arg pypath.go.GOAnnotation go_annot:
            A :class:``pypath.go.GOAnnotation`` object.
        """

        session_mod.Logger.__init__(self, name = 'go')

        self.go_annot = go_annot or get_db() # TODO: consider ncbi_tax_id at
                                             # selection DB

        self._categories = categories
        self.process_categories()


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def process_categories(self):
        """
        Translates GO term names listed in categories to GO terms ACs.
        """

        # if the categories are grouped by aspects
        if (
            isinstance(self._categories, dict) and
            not set(self._categories.keys()) - set(self.go_annot.aspects)
        ):

            if isinstance(list(self._categories.values())[0], set):

                self._categories = set.union(*self._categories.values())

            elif isinstance(list(self._categories.values())[0], dict):

                self._categories = dict(
                    itertools.chain(
                        *(d.items() for d in self._categories.values())
                    )
                )

        # if a set provided we use names as keys
        # and accessions as values
        if isinstance(self._categories, set):

            self._categories = dict(
                (
                    self.go_annot.get_name(cat),
                    self.go_annot.get_term(cat)
                )
                for cat in self._categories
            )

        self.categories = self._categories


    def get_annotation(self, category, uniprots = None):
        """
        For a category name returns a set of UniProt IDs annotated with
        the corresponding Gene Ontology terms or expression.

        :arg str category:
            The category name, should be a key in the ``categories`` dict.
        :arg set uniprots:
            A set or list of UniProt IDs. If ``None``, annotations based on
            all UniProts in GO annotation will be returned.
        """

        return self.go_annot.select(
            self.categories[category],
            uniprots = uniprots,
            return_uniprots = True,
        )


    def get_annotations(self, uniprots = None):
        """
        Returns a dict with set of UniProt IDs for each category.

        :arg set uniprots:
            A set or list of UniProt IDs. If ``None``, annotations based on
            all UniProts in GO annotation will be returned.
        """

        return dict(
            (
                category,
                self.get_annotation(category, uniprots = uniprots)
            )
            for category in self.categories.keys()
        )


def annotate(graph, organism = 9606, aspects = ('C', 'F', 'P')):
    """
    Adds Gene Ontology annotations to the nodes of a graph.

    :param igraph.Graph graph:
        Any ``igraph.Graph`` object with uniprot IDs
        in its ``name`` vertex attribute.
    """

    aspects = aspects if type(aspects) in {list, tuple} else (aspects, )

    graph.vs['go'] = [
        {'C': set(), 'F': set(), 'P': set()}
        for _ in xrange(graph.vcount())
    ]

    terms, annot = go_input.go_annotations_goa(organism = organism)

    prg = progress.Progress(graph.vcount(), 'Loading GO annotations', 9)

    for v in graph.vs:

        prg.step()

        for asp in aspects:

            if v['name'] in annot[asp]:

                v['go'][asp] = annot[asp][v['name']]

    prg.terminate()


# old name as synonym
load_go = annotate


def init_db(organism = 9606, pickle_file = None, use_pickle_cache = True):
    """
    Initializes or reloads the GO annotation database.
    The database will be assigned to the ``db`` attribute of this module.
    """

    if 'db' not in globals():

        globals()['db'] = {}

    globals()['db'][organism] = GOAnnotation(
        organism,
        pickle_file = pickle_file,
        use_pickle_cache = use_pickle_cache,
    )


def get_db(organism = 9606, pickle_file = None, use_pickle_cache = True):
    """
    Retrieves the current database instance and initializes it if does
    not exist yet.
    """

    # TODO: consider organism
    # TODO: delete the DB if not used in order to free memory
    # TODO: introduce pickle cache to make it load quicker

    if 'db' not in globals() or organism not in globals()['db']:

        init_db(
            organism,
            pickle_file = pickle_file,
            use_pickle_cache = use_pickle_cache,
        )

    return globals()['db'][organism]
