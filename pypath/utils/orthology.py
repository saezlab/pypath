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

from __future__ import annotations

from future.utils import iteritems
from past.builtins import xrange, range
from typing import Literal

import os
import sys
import itertools
import functools
import collections
import importlib as imp
import re
import time
import datetime
import json
import pickle
import copy
import abc
import inspect
import types as _types

import timeloop
import pandas as pd

import pypath.utils.mapping as mapping
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.data as _data
import pypath.internals.intera as intera
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.uniprot_db as uniprot_db
import pypath.inputs.homologene as homologene_input
import pypath.inputs.oma as oma_input
import pypath.inputs.biomart as biomart
import pypath.utils.seq as _se
import pypath.share.session as session
import pypath.share.settings as settings
import pypath.utils.taxonomy as taxonomy
import pypath.share.cache as cache_mod

_orthology_cleanup_timeloop = timeloop.Timeloop()
_orthology_cleanup_timeloop.logger.setLevel(9999)

_logger = session.Logger(name = 'orthology')
_log = _logger._log


class OrthologMeta(abc.ABCMeta):


    def __new__(
            cls,
            name,
            supercls,
            attrs,
            resource: str,
            fields: tuple[str] = ('id',),
        ):

        Base = collections.namedtuple(f'{name}Base', fields)

        def new_method(cls, *args, **kwargs):

            return Base.__new__(cls, *args, **kwargs)

        sig = inspect.signature(new_method)
        new_params = [
            inspect.Parameter(arg, inspect.Parameter.POSITIONAL_OR_KEYWORD)
            for arg in fields
        ]
        new_method.__signature__ = sig.replace(
            parameters = list(itertools.chain(
                (sig.parameters['cls'],),
                new_params,
            ))
        )

        clazz = _types.new_class(name, (Base,))
        setattr(clazz, '__new__', new_method)
        setattr(clazz, '_resource', resource)
        clazz.__module__ = __name__

        return clazz


class OrthologBase:


    def __str__(self):

        return self.id


    def __repr__(self):

        return f'<Ortholog {self.id} ({self._resource})>'


    def __eq__(self, other):

        return self.__str__() == other.__str__()


    def __hash__(self):

        return self.id.__hash__()


class OmaOrtholog(
        OrthologBase,
        metaclass = OrthologMeta,
        resource = 'OMA',
        fields = ('id', 'rel_type', 'score'),
    ): pass



class EnsemblOrtholog(
        OrthologBase,
        metaclass = OrthologMeta,
        resource = 'Ensembl',
        fields = ('id', 'types', 'hc'),
    ): pass


class HomologeneOrtholog(
        OrthologBase,
        metaclass = OrthologMeta,
        resource = 'HomoloGene',
    ): pass


OrthologyTableKey = collections.namedtuple(
    'OrthologyTableKey',
    ('source', 'target', 'only_swissprot', 'resource', 'id_type'),
)

class OrthologyManager(session.Logger):

    TRANSLATION_PARAM = (
        'oma',
        'homologene',
        'oma_rel_type',
        'oma_score',
        'ensembl',
        'ensembl_hc',
        'ensembl_types',
    )
    RESOURCE_PARAM = {
        'oma': ('rel_type', 'score'),
        'ensembl': ('hc', 'types'),
        'homologene': (),
    }

    def __init__(
            self,
            cleanup_period: int = 10,
            lifetime: int = 300,
            **kwargs
        ):

        session.Logger.__init__(self, name = 'orthology')


        @_orthology_cleanup_timeloop.job(
            interval = datetime.timedelta(
                seconds = cleanup_period
            )
        )
        def _cleanup():

            self._remove_expired()


        _orthology_cleanup_timeloop.start(block = False)

        self.lifetime = lifetime
        self.tables = {}
        self.expiry = {}
        self._param = {k: kwargs.get(k, None) for k in self.TRANSLATION_PARAM}
        self._log('OrthologyManager has been created.')


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def which_table(
            self,
            target: str | int,
            source: str | int = 9606,
            only_swissprot: bool = True,
            resource: Literal['oma', 'homologene', 'ensembl'] = 'oma',
            id_type: str = 'uniprot',
        ):

        loc = locals()
        key = OrthologyTableKey(**{
            f: loc[f]
            for f in OrthologyTableKey._fields
        })

        self.expiry[key] = time.time()

        if key not in self.tables:

            self.load(key)

        if key in self.tables:

            return self.tables[key]


    def load(self, key):

        self.tables[key] = globals()[f'{key.resource.capitalize()}Orthology'](
            target = key.target,
            source = key.source,
            only_swissprot = key.only_swissprot,
            id_type = key.id_type,
        )


    @common.ignore_unhashable
    @functools.lru_cache(maxsize = int(1e5))
    def translate(
            self,
            identifiers: str | Iterable[str],
            target: str | int,
            source: str | int = 9606,
            id_type: str = 'uniprot',
            only_swissprot: bool = True,
            oma: bool = None,
            homologene: bool = None,
            ensembl: bool = None,
            oma_rel_type: (
                set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
                None
            ) = None,
            oma_score: float | None = None,
            ensembl_hc: bool = True,
            ensembl_types: (
                list[Literal['one2one', 'one2many', 'many2many']] |
                None
            ) = None,
            full_records: bool = False,
        ):
        """
        Translate one or more identifiers by orthologous gene pairs.

        Args:
            identifiers:
                One or more identifers of the source organism, of ID type
                `id_type`.
            target:
                Name or NCBI Taxonomy ID of the target organism.
            source:
                Name or NCBI Taxonomy ID of the source organism.
            id_type:
                The identifier type to use.
            only_swissprot:
                Use only SwissProt IDs.
            oma
                Use orthology information from the Orthologous Matrix (OMA).
                Currently this is the recommended source for orthology data.
            homologene:
                Use orthology information from NCBI HomoloGene.
            ensembl:
                Use orthology information from Ensembl.
            oma_rel_type:
                Restrict relations to certain types.
            oma_score:
                Lower threshold for similarity metric.
            ensembl_hc:
                Use only the high confidence orthology relations from Ensembl.
            ensembl_types:
                Ensembl orthology relation types to use. Possible values are
                `one2one`, `one2many` and `many2many`. By default only
                `one2one` is used.
            full_records:
                Include not only the identifiers, but also some properties of
                the orthology relationships.

        Returns:
            Set of identifiers of orthologous genes or proteins in the
            target taxon.
        """

        target = taxonomy.ensure_ncbi_tax_id(target)
        source = taxonomy.ensure_ncbi_tax_id(source)
        param = self._translation_param(locals())
        proc = (lambda x: x) if full_records else (lambda x: x.id)

        result = set()

        for resource, keys in self.RESOURCE_PARAM.items():

            if not param[resource]:

                continue

            table = self.which_table(
                target = target,
                source = source,
                only_swissprot = only_swissprot,
                id_type = id_type,
                resource = resource,
            )

            result.update(
                table.translate(
                    identifiers,
                    full_records = full_records,
                    **{k: v for k, v in param.items() if k in keys},
                )
            )

        return result


    def get_dict(
            self,
            target: str | int,
            source: str | int = 9606,
            id_type: str = 'uniprot',
            only_swissprot: bool = True,
            oma: bool = None,
            homologene: bool = None,
            ensembl: bool = None,
            oma_rel_type: (
                set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
                None
            ) = None,
            oma_score: float | None = None,
            ensembl_hc: bool = True,
            ensembl_types: (
                list[Literal['one2one', 'one2many', 'many2many']] |
                None
            ) = None,
            full_records: bool = False,
        ) -> dict[str, set[OrthologBase]]:
        """
        Create a dictionary for one source organism and ID type.

        Args:
            target:
                Name or NCBI Taxonomy ID of the target organism.
            source:
                Name or NCBI Taxonomy ID of the source organism.
            id_type:
                The identifier type to use.
            only_swissprot:
                Use only SwissProt IDs.
            oma
                Use orthology information from the Orthologous Matrix (OMA).
                Currently this is the recommended source for orthology data.
            homologene:
                Use orthology information from NCBI HomoloGene.
            ensembl:
                Use orthology information from Ensembl.
            oma_rel_type:
                Restrict relations to certain types.
            oma_score:
                Lower threshold for similarity metric.
            ensembl_hc:
                Use only the high confidence orthology relations from Ensembl.
            ensembl_types:
                Ensembl orthology relation types to use. Possible values are
                `one2one`, `one2many` and `many2many`. By default only
                `one2one` is used.
            full_records:
                Include not only the identifiers, but also some properties of
                the orthology relationships.

        Returns:
            A dict with identifiers of the source organism as keys, and
            sets of their orthologs as values.
        """

        target = taxonomy.ensure_ncbi_tax_id(target)
        source = taxonomy.ensure_ncbi_tax_id(source)
        param = self._translation_param(locals())

        result = collections.defaultdict(set)

        for resource, keys in self.RESOURCE_PARAM.items():

            if not param[resource]:

                continue

            table = self.which_table(
                target = target,
                source = source,
                only_swissprot = only_swissprot,
                id_type = id_type,
                resource = resource,
            )
            dct = table.asdict(
                full_records = full_records,
                **{
                    p: v
                    for p, v in param.items()
                    if p in keys
                }
            )

            for s, o in dct.items():

                result[s].update(o)

        return dict(result)


    def get_df(
            self,
            target: str | int,
            source: str | int = 9606,
            id_type: str = 'uniprot',
            only_swissprot: bool = True,
            oma: bool = None,
            homologene: bool = None,
            ensembl: bool = None,
            oma_rel_type: (
                set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
                None
            ) = None,
            oma_score: float | None = None,
            ensembl_hc: bool = True,
            ensembl_types: (
                list[Literal['one2one', 'one2many', 'many2many']] |
                None
            ) = None,
            full_records: bool = False,
            **kwargs
        ) -> pd.DataFrame:
        """
        Create a data frame for one source organism and ID type.

        Args:
            target:
                Name or NCBI Taxonomy ID of the target organism.
            source:
                Name or NCBI Taxonomy ID of the source organism.
            id_type:
                The identifier type to use.
            only_swissprot:
                Use only SwissProt IDs.
            oma
                Use orthology information from the Orthologous Matrix (OMA).
                Currently this is the recommended source for orthology data.
            homologene:
                Use orthology information from NCBI HomoloGene.
            ensembl:
                Use orthology information from Ensembl.
            oma_rel_type:
                Restrict relations to certain types.
            oma_score:
                Lower threshold for similarity metric.
            ensembl_hc:
                Use only the high confidence orthology relations from Ensembl.
            ensembl_types:
                Ensembl orthology relation types to use. Possible values are
                `one2one`, `one2many` and `many2many`. By default only
                `one2one` is used.
            full_records:
                Include not only the identifiers, but also some properties of
                the orthology relationships.
            kwargs:
                Ignored.

        Returns:
            A data frame with pairs of orthologous identifiers,
            in two columns: "source" and "target".
        """

        target = taxonomy.ensure_ncbi_tax_id(target)
        source = taxonomy.ensure_ncbi_tax_id(source)
        param = self._translation_param(locals())

        result = []

        for resource, keys in self.RESOURCE_PARAM.items():

            if not param[resource]:

                continue

            table = self.which_table(
                target = target,
                source = source,
                only_swissprot = only_swissprot,
                id_type = id_type,
                resource = resource,
            )

            result.append(
                table.df(
                    full_records = full_records,
                    **{
                        p: v
                        for p, v in param.items()
                        if p in keys
                    }
                )
            )

        return pd.concat(result)


    def translate_df(
            self,
            df: pd.DataFrame,
            target: str | int,
            source: str | int = 9606,
            cols: str | list[str] | dict[str, str] | None = None,
            id_type: str = 'uniprot',
            only_swissprot: bool = True,
            oma: bool = None,
            homologene: bool = None,
            ensembl: bool = None,
            oma_rel_type: (
                set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
                None
            ) = None,
            oma_score: float | None = None,
            ensembl_hc: bool = True,
            ensembl_types: (
                list[Literal['one2one', 'one2many', 'many2many']] |
                None
            ) = None,
            **kwargs: str | tuple[str, str]
        ) -> pd.DataFrame:
        """
        Translate columns in a data frame.

        Args:
            df:
                A data frame.
            cols:
                One or more columns to be translated. It can be a single
                column name, an iterable of column names or a dict where
                keys are column names and values are ID types. Except this
                last case, identifiers are assumed to be `id_type`.
            target:
                Name or NCBI Taxonomy ID of the target organism.
            source:
                Name or NCBI Taxonomy ID of the source organism.
            id_type:
                The default identifier type to use, will be used for all
                columns where ID type is not specified.
            only_swissprot:
                Use only SwissProt IDs.
            oma
                Use orthology information from the Orthologous Matrix (OMA).
                Currently this is the recommended source for orthology data.
            homologene:
                Use orthology information from NCBI HomoloGene.
            ensembl:
                Use orthology information from Ensembl.
            oma_rel_type:
                Restrict relations to certain types.
            oma_score:
                Lower threshold for similarity metric.
            ensembl_hc:
                Use only the high confidence orthology relations from Ensembl.
            ensembl_types:
                Ensembl orthology relation types to use. Possible values are
                `one2one`, `one2many` and `many2many`. By default only
                `one2one` is used.
            kwargs:
                Same as providing a dict to ``cols``, but beware, keys
                (column names) can not match existing argument names of
                this function.

        Returns:
            A data frame with the same column layout as the input, and the
            identifiers translated as demanded. Rows that could not be
            translated are omitted.
        """

        if not isinstance(cols, dict):

            cols = dict((col, id_type) for col in common.to_list(cols))

        kwargs.update(cols)
        id_types = set(kwargs.values())

        for _id_type in set(cols.values()):

            args = locals().copy()
            args.pop('self')
            args['id_type'] = _id_type
            ortho_df = self.get_df(**args)

            table = self.which_table(
                target = target,
                source = source,
                only_swissprot = only_swissprot,
                id_type = _id_type,
                resource = 'oma',
            )

            df = table.translate_df(
                df = df,
                cols = [c for c, i in cols.items() if i == _id_type],
                ortho_df = ortho_df,
            )

        return df


    def _translation_param(self, loc: dict) -> dict:

        param = {}

        for resource, keys in self.RESOURCE_PARAM.items():

            enabled = common.first_value(
                loc[resource],
                self._param[resource],
                settings.get(f'orthology_{resource}'),
            )
            param[resource] = enabled

            if enabled:

                for key in keys:

                    param[key] = common.first_value(
                        loc[f'{resource}_{key}'],
                        self._param[f'{resource}_{key}'],
                        settings.get(f'orthology_{resource}_{key}'),
                    )

        return param


    def _remove_expired(self):

        for key, last_used in list(self.expiry.items()):

            if time.time() - last_used > self.lifetime and key in self.tables:

                self._log(
                    'Removing orthology table from taxon %u to %u '
                    '(only SwissProt: %s; resource: %s; ID type: %s)' % key
                )

                del self.tables[key]
                del self.expiry[key]


    def __del__(self):

        if hasattr(_orthology_cleanup_timeloop, 'stop'):

            _orthology_cleanup_timeloop.stop()


class SequenceContainer(session.Logger):


    def __init__(self, preload_seq = [], isoforms = True):
        """
        This is an object to store sequences of multiple
        organisms and select the appropriate one.
        """

        if not hasattr(self, '_logger'):

            session.Logger.__init__(self, name = 'orthology')

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


    def __init__(
            self,
            preload_prot: list[int] | None = None,
            only_swissprot: bool = True,
        ):

        if not hasattr(self, '_taxonomy'):

            self._taxonomy = {}
            self._up_taxonomy = {}
            self._proteomes = {}

        self.only_swissprot = only_swissprot
        self.load_taxonomy()

        for taxon in (preload_prot or ()):

            self.load_proteome(taxon)


    def load_proteome(self, taxon: int, only_swissprot: bool | None = None):

        only_swissprot = (
            self.only_swissprot
                if only_swissprot is None else
            only_swissprot
        )

        key = (taxon, only_swissprot)

        if key not in self._proteomes:

            self._proteomes[key] = (
                set(uniprot_db.all_uniprots(*key))
            )

            for protein in self._proteomes[key]:

                self._taxonomy[protein] = key

            if not only_swissprot:

                self.load_proteome(taxon, True)


    def get_taxon(self, protein, only_swissprot = True):

        ncbi_tax_id = self.get_taxon_trembl(protein)

        if (
            only_swissprot and
            ncbi_tax_id and
            not uniprot_db.is_swissprot(protein, organism = ncbi_tax_id)
        ):

            ncbi_tax_id = None

        return ncbi_tax_id


    def get_taxon_trembl(self, protein):

        return self._up_taxonomy.get(protein, None)


    def has_protein(self, protein):

        return protein in self._taxonomy


    def is_swissprot(self, protein):

        return bool(self.get_taxon(protein, only_swissprot = True))


    def load_taxonomy(self):

        self._up_taxonomy = uniprot_input.uniprot_taxonomy(ncbi_tax_ids = True)


class ProteinOrthology(Proteomes):

    _param = ('id',)


    def __init__(
            self,
            target: str | int,
            source: str | int | None = 9606,
            id_type: str = 'uniprot',
            only_swissprot: bool = True,
            **kwargs
        ):
        """
        This class translates between homologous UniProt IDs of two organisms
        based on NCBI HomoloGene and Ensembl data. In case of HomoloGene,
        the UniProt-UniProt translation table is created by translating the
        source organism UniProts to RefSeq and Entrez IDs, finding the
        homologues (orthologues) for these IDs, and then translating them
        to the target organism UniProt IDs. In case of Ensembl, we obtain
        data with Ensembl protein identifiers and translate those to UniProt.

        Args
            target:
                Name or NCBI Taxonomy ID of the target organism.
            source:
                Name or NCBI Taxonomy ID of the source organism.
            id_type:
                The identifier type to use.
            only_swissprot:
                Use only SwissProt IDs.
            kwargs:
                Resource specific parameters.
        """

        self.data = {}
        self.target = taxonomy.ensure_ncbi_tax_id(target)
        self.source = taxonomy.ensure_ncbi_tax_id(source)
        self.id_type = id_type
        self._resource_l = self.resource.lower()
        Proteomes.__init__(self, only_swissprot = only_swissprot)
        self.load_proteome(self.source)
        self._set_param(kwargs, *self._param)
        self.load()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def load(self, source = None):

        pass


    def translate(
            self,
            identifier: str | Iterable[str],
            full_records: bool = False,
            **kwargs
        ) -> set[str]:
        """
        For one UniProt ID of the source organism returns all orthologues
        from the target organism.

        Args:
            identifier:
                An identifier corresponding to the ID type and source organism
                of the instance.
            full_records:
                Include not only the identifiers, but also some properties of
                the orthology relationships.
            kwargs:
                Resource specific translation parameters.

        Returns:
            A set of identifiers of orthologues in the target taxon.
        """

        identifier = (
            (identifier,)
                if hasattr(identifier, 'components') else
            common.to_list(identifier)
        )

        result = set.union(*(self.data.get(i, set()) for i in identifier))

        if not full_records:

            result = {o.id for o in result}

        return result


    def asdict(
            self,
            full_records: bool = False,
            **kwargs
        ) -> dict[str, set[OrthologBase]]:
        """
        Create a dictionary from the translation table.

        Args:
            full_records:
                Include not only the identifiers, but also some properties of
                the orthology relationships.
            kwargs:
                Resource specific filtering criteria.

        Returns:
            A dict with identifiers of the source organism as keys, and
            sets of their orthologs as values.
        """

        proc = (lambda x: x) if full_records else (lambda x: x.id)

        return {
            s: {proc(o) for o in orthologs if self.match(o, **kwargs)}
            for s, orthologs in self.data.items()
        }


    def df(self, full_records: bool = False, **kwargs) -> pd.DataFrame:
        """
        Orthologous pairs as data frame.

        Args:
            full_records:
                Include not only the identifiers, but also some properties of
                the orthology relationships.
            kwargs:
                Resource specific filtering criteria.

        Returns:
            A data frame with pairs of orthologous identifiers,
            in two columns: "source" and "target".
        """

        _log(
            'Creating translation data frame between '
            f'organisms `{self.source}` and `{self.target}`, '
            f'ID type `{self.id_type}`.'
        )

        df = (
            pd.DataFrame(
                self.asdict(full_records = full_records, **kwargs).items(),
                columns = ['source', 'target'],
            ).
            explode('target', ignore_index = True).
            dropna().
            reset_index(drop = True)
        )

        if full_records:

            # some beautiful pandas code again
            df = (
                pd.concat(
                    [
                        df.source,
                        pd.DataFrame(df.target.tolist()),
                    ],
                    axis = 1,
                ).
                rename(columns = {'id': 'target'})
            )

        return df


    def translate_df(
            self,
            df: pd.DataFrame,
            cols: str | list[str] | None = None,
            ortho_df: pd.DataFrame | None = None,
            **kwargs
        ):
        """
        Translate columns in a data frame.

        Args:
            df:
                A data frame.
            cols:
                One or more columns to be translated. It can be a single
                column name, an iterable of column names or a dict where
                keys are column names and values are ID types. Except this
                last case, identifiers are assumed to be UniProt.
            ortho_df:
                Override the translation data frame. If provided, the
                parameters in `kwargs` won't have an effect. Must have
                columns "source" and "target".
            kwargs:
                Resource specific translation parameters.

        Returns:
            A data frame with the same column layout as the input, and the
            identifiers translated as demanded. Rows that could not be
            translated are omitted.
        """

        _log(
            f'Translating data frame column(s) from '
            f'organism `{self.source}` to `{self.target}`.'
        )

        ortho_df = (
            (self.df(**kwargs) if ortho_df is None else ortho_df).
            rename(
                {
                    'source': 'pypath_internal_source',
                    'target': 'pypath_internal_target',
                },
                axis = 1,
            )
        )

        col_order = df.columns
        cols = common.to_list(cols)

        for col in cols:

            _log(
                f'Translating `{self.id_type}` IDs of organism `{self.source}` '
                f'in column `{col}` to organism `{self.target}`.'
            )

            df = (
                df.merge(
                    ortho_df.rename({'pypath_internal_source': col}, axis = 1),
                    on = col,
                    how = 'inner',
                ).
                drop(col, axis = 1).
                rename({'pypath_internal_target': col}, axis = 1)
            )

        return df[col_order]


    def _translation_param(self, loc: dict) -> dict:

        return dict(
            (p, loc[p])
            for p in OrthologyManager.TRANSLATION_PARAM
        )


    def _set_param(self, loc: dict, *params: str):

        for param in params:

            key = f'orthology_{self._resource_l}_{param}'
            setattr(
                self,
                param,
                common.first_value(loc.get(param, None), settings.get(key)),
            )


    def match(self, ortholog: OrthologBase, **kwargs) -> bool:

        return True


    def _from_pickle(self) -> bool:

        if (
            settings.get('orthology_cache') and
            os.path.exists(self.pickle_path)
        ):

            with open(self.pickle_path, 'rb') as fp:

                self.data = pickle.load(fp)
                _log(
                    'Orthology table from taxon %u to %u (only SwissProt: %s; '
                    'resource: %s; ID type: %s) has been loaded from `%s`.' % (
                        self.key + (self.pickle_path,)
                    )
                )
                return True

        return False


    def _to_pickle(self):

        with open(self.pickle_path, 'wb') as fp:

            pickle.dump(self.data, fp)
            _log(
                'Orthology table from taxon %u to %u (only SwissProt: %s; '
                'resource: %s; ID type: %s) has been saved to `%s`.' % (
                    self.key + (self.pickle_path,)
                )
            )


    @property
    def key(self):

        return OrthologyTableKey(
            source = self.source,
            target = self.target,
            only_swissprot = self.only_swissprot,
            resource = self._resource_l,
            id_type = self.id_type,
        )


    @property
    def pickle_path(self):

        return os.path.join(
            cache_mod.get_cachedir(),
            f'{common.md5(json.dumps(self.key))}.pickle',
        )


    def __len__(self):

        return sum(map(len, self.data.values()))


    def __repr__(self):

        return (
            f'<{self.resource} Orthology table from {self.source} to '
            f'{self.target}: {self.id_type} IDs, {len(self)} relationships>'
        )


class HomologeneOrthology(ProteinOrthology):

    resource = 'HomoloGene'


    def load(self):
        """
        Load orthology data from NCBI HomoloGene.

        Builds orthology translation table as dict based on NCBI HomoloGene
        data. If the `id_type` is supported by HomoloGene (Gene Symbol, RefSeq,
        Entrez, GI), the data will be simply loaded. For other ID types it
        translates HomoloGene Gene Symbol, RefSeq and Entrez tables to UniProt
        and then translates the orthologous UniProt pairs to the desired
        ID type.
        """

        if self._from_pickle():

            return

        if self.id_type in ('genesymbol', 'refseq', 'refseqp', 'entrez', 'gi'):

            data = homologene_input.homologene_dict(
                self.source,
                self.target,
                self.id_type,
            )
            self.data = {
                s: {HomologeneOrtholog(t) for t in target_ids}
                for s, target_ids in data.items()
            }
            return

        hg = {
            id_type: homologene_input.homologene_dict(
                self.source,
                self.target,
                id_type,
            )
            for id_type in ('genesymbol', 'refseq', 'entrez')
        }

        _log(
            'Loading orthology data from NCBI HomoloGene '
            f'between organisms `{self.source}` and `{self.target}`.'
        )

        self.data = collections.defaultdict(set)

        for u in self._proteomes[(self.source, self.only_swissprot)]:

            target_uniprots = set()

            for id_type, hgdata in hg.items():

                hg_source_ids = mapping.map_name(
                    u,
                    'uniprot',
                    id_type,
                    ncbi_tax_id = self.source,
                )

                if not hg_source_ids: continue

                hg_target_ids = set.union(*(
                    hgdata.get(s, set())
                    for s in hg_source_ids)
                )

                if not hg_target_ids: continue

                target_uniprots.update(
                    mapping.map_names(
                        hg_target_ids,
                        id_type,
                        'uniprot',
                        ncbi_tax_id = self.target,
                    )
                )

            if self.id_type == 'uniprot':

                source_ids = (u,)
                target_ids = target_uniprots

            else:

                source_ids = mapping.map_name(
                    u,
                    'uniprot',
                    self.id_type,
                    ncbi_tax_id = self.source,
                )
                target_ids = mapping.map_names(
                    target_uniprots,
                    'uniprot',
                    self.id_type,
                    ncbi_tax_id = self.target,
                )

            for s in source_ids:

                self.data[s].update(
                    {
                        HomologeneOrtholog(t)
                        for t in target_ids
                    }
                )

        self.data = dict(self.data)
        self._to_pickle()


class EnsemblOrthology(ProteinOrthology):

    _param = ('hc', 'types')
    resource = 'Ensembl'

    def __init__(
            self,
            target: int | str,
            source: int | str = 9606,
            id_type: str = 'uniprot',
            only_swissprot: bool = None,
            hc: bool = None,
            types: list[Literal[
                'one2one', 'one2many', 'many2many'
            ]] = None,
        ):
        """
        Orthology translation with Ensembl data.

        Args
            target:
                Name or NCBI Taxonomy ID of the target organism.
            source:
                Name or NCBI Taxonomy ID of the source organism.
            id_type:
                The identifier type to use.
            only_swissprot:
                Use only SwissProt IDs.
            hc:
                Use only high confidence orthology relations
                from Ensembl. By default it is True. You can also set it
                by the `ensembl_hc` attribute.
            types:
                The Ensembl orthology relationship types
                to use. Possible values are `one2one`, `one2many` and
                `many2many`. By default only `one2one` is used. You can
                also set this parameter by the `ensembl_types` attribute.
        """

        ProteinOrthology.__init__(**locals())


    def load(self):

        target_organism = taxonomy.ensure_ensembl_name(self.target)
        source_organism = taxonomy.ensure_ensembl_name(self.source)

        _log(
            'Loading orthology data from Ensembl '
            f'between organisms `{self.source}` and `{self.target}`.'
        )

        if self._from_pickle():

            return

        if not target_organism or not source_organism:

            _log(
                'No Ensembl orthology data available between '
                f'organisms `{self.source}` and `{self.target}`.'
            )
            return

        target_prefix = f'{target_organism}_homolog_'

        attr_target_ensp = f'{target_prefix}ensembl_peptide'
        attr_conf = f'{target_prefix}orthology_confidence'
        attr_type = f'{target_prefix}orthology_type'

        ensembl_data = biomart.biomart_homology(
            source_organism = self.source,
            target_organism = self.target,
        )

        _id_types = {
            'target': {
                'genesymbol': f'{target_prefix}associated_gene_name',
                'ensp': f'{target_prefix}ensembl_peptide',
                'ensg': f'{target_prefix}ensembl_gene',
            },
            'source': {
                'genesymbol': 'external_gene_name',
                'ensp': 'ensembl_peptide_id',
                'ensg': 'ensembl_gene_id',
            },
        }
        attr_tgt_id = _id_types['target'].get(
            self.id_type,
            f'{target_prefix}ensembl_peptide',
        )
        attr_src_id = _id_types['source'].get(
            self.id_type,
            'ensembl_peptide_id',
        )

        self.data = collections.defaultdict(set)

        if self.id_type in _id_types['target']:

            for r in ensembl_data:

                self.data[getattr(r, attr_src_id)].add(
                    EnsemblOrtholog(
                        id = getattr(r, attr_tgt_id),
                        hc = getattr(r, attr_conf) == '1',
                        types = getattr(r, attr_type).split('_')[-1],
                    )
                )

        for r in ensembl_data:

            ids = {}

            for side, attr_id in (
                    ('source', attr_src_id),
                    ('target', attr_tgt_id)
                ):

                uniprots = mapping.map_name(
                    getattr(r, attr_id),
                    'ensp',
                    'uniprot',
                    ncbi_tax_id = getattr(self, side),
                )
                ids[side] = mapping.map_names(
                    uniprots,
                    'uniprot',
                    self.id_type,
                    ncbi_tax_id = getattr(self, side),
                    uniprot_cleanup = False,
                )

                if not ids[side]:

                    continue

            for s in ids['source']:

                self.data[s].update(
                    {
                        EnsemblOrtholog(
                            id = t,
                            hc = getattr(r, attr_conf) == '1',
                            types = getattr(r, attr_type).split('_')[-1],
                        )
                        for t in ids['target']
                    }
                )

        self.data = dict(self.data)
        self._to_pickle()


    def match(self, ortholog: OrthologBase, **kwargs) -> bool:
        """
        Check an ortholog against filtering criteria.

        Args:
            ortholog:
                An ortholog record.
            kwargs:
                Override default filtering parameters.

        Returns:
            True if the ortholog meets the criteria.
        """

        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        hc = kwargs.get('hc', self.hc)
        types = kwargs.get('types', self.types)

        return (
            (not hc or ortholog.hc) and
            (not types or ortholog.types in types)
        )


class OmaOrthology(ProteinOrthology):

    _param = ('rel_type', 'score')
    resource = 'OMA'

    def __init__(
            self,
            target: int | str,
            source: int | str = 9606,
            id_type: str = 'uniprot',
            only_swissprot: bool = None,
            rel_type: (
                set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
                None
             ) = None,
            score: float | None = None,
        ):
        """
        Orthology translation with Ensembl data.

        Args
            target:
                Name or NCBI Taxonomy ID of the target organism.
            source:
                Name or NCBI Taxonomy ID of the source organism.
            id_type:
                The identifier type to use.
            only_swissprot:
                Use only SwissProt IDs.
            rel_type:
                Restrict relations to certain types.
            score:
                Lower threshold for similarity metric.
        """

        ProteinOrthology.__init__(**locals())


    def load(self):

        _log(
            'Loading orthology data from OMA '
            f'between organisms `{self.source}` and `{self.target}`.'
        )

        if self._from_pickle():

            return

        oma_data = oma_input.oma_orthologs(
            organism_a = self.source,
            organism_b = self.target,
            id_type = self.id_type,
        )
        self.data = collections.defaultdict(set)

        for rec in oma_data:

            self.data[rec.a.id].add(
                OmaOrtholog(
                    id = rec.b.id,
                    score = rec.score,
                    rel_type = rec.rel_type,
                )
            )

        self.data = dict(self.data)
        self._to_pickle()


    def match(self, ortholog: OrthologBase, **kwargs) -> bool:
        """
        Check an ortholog against filtering criteria.

        Args:
            ortholog:
                An ortholog record.
            kwargs:
                Override default filtering parameters.

        Returns:
            True if the ortholog meets the criteria.
        """

        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        score = kwargs.get('score', self.score)
        rel_type = kwargs.get('rel_type', self.rel_type)

        return (
            (score is None or ortholog.score >= score) and
            (not rel_type or ortholog.rel_type in rel_type)
        )


class PtmOrthology(Proteomes, SequenceContainer):


    def __init__(
            self,
            target: str | int,
            source: str | int | None = None,
            only_swissprot: bool = True,
            strict: bool = True,
            orthology_args: dict | None = None,
        ):

        if not hasattr(self, '_logger'):

            session.Logger.__init__(self, name = 'orthology')

        self.manager = get_manager()
        SequenceContainer.__init__(self)
        Proteomes.__init__(self, only_swissprot = only_swissprot)
        self.source = taxonomy.ensure_ncbi_tax_id(source)
        self.target = taxonomy.ensure_ncbi_tax_id(target)
        self.load_seq(taxon = self.target)
        self.reptm = re.compile(r'([A-Z\d]{6,10})_([A-Z])(\d*)')
        self.strict = strict
        self.orthology_args = orthology_args or {}
        self.id_type = 'uniprot'
        self.ptm_orthology()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def translate_site(
            self,
            protein: str | intera.Protein,
            res: str,
            offset: int,
            isoform: int = 1,
            typ: str = 'phosphorylation',
            source_organism: str | int | None = None,
        ) -> set[tuple]:
        """
        Translates one PTM site.

        Args:
            protein:
                A protein identifier or an intera.Protein object.
            res:
                Single letter code of the residue.
            offset:
                Sequence offset of the site.
            isoform:
                Sequence isoform.
            typ:
                Modification type.
            source_organism:
                Name or NCBI Taxonomy ID of the source organism.

        Returns:
            A list of tuples with the identifier, isoform, residue, offset,
            taxon and modification type of the orthologous PTM sites.
        """

        result = set()
        source = self._get_source(source_organism)
        protein_id = getattr(protein, 'identifier', protein)
        sourceptm = (protein_id, isoform, res, offset, source, typ)

        if self.get_taxon(protein_id) == self.target:
            result.add(sourceptm)
            return result

        if sourceptm in self.ptmortho:

            if self.target in self.ptmortho[sourceptm]:

                result = self.ptmortho[sourceptm]

        if not result and not self.strict:

            tsubs = self.manager.translate(
                identifiers = protein_id,
                target = self.target,
                source = source,
                only_swissprot = self.only_swissprot,
                id_type = self.id_type,
                **self.orthology_args
            )

            for tsub in tsubs:

                se = self.get_seq(tsub, taxon = self.target)

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
                                typ,
                            ))

                    if result:
                        break

        return result


    def translate_domain(self, domain: intera.Domain) -> list[intera.Domain]:

        return [
            intera.Domain(
                protein = target_id,
                ncbi_tax_id = self.target,
            )
            for target_id in self.manager.translate(
                identifiers = domain.protein.identifier,
                target = self.target,
                id_type = self.id_type,
                source = self.get_source(domain.ncbi_tax_id),
                only_swissprot = self.only_swissprot,
                **self.orthology_args
            )
        ]


    def translate_ptm(self, ptm: intera.Ptm) -> list[intera.Ptm]:

        tptms = self.translate_site(
            ptm.protein,
            ptm.residue.name,
            ptm.residue.number,
            ptm.residue.isoform,
            ptm.typ,
        )

        result = []

        for x in tptms:

            se = self.get_seq(x[0], taxon = self.target)

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


    def translate_domain_motif(
            self,
            dmotif: intera.DomainMotif,
        ) -> list[intera.DomainMotif]:

        ds = self.translate_domain(dmotif.domain)
        ps = self.translate_ptm(dmotif.ptm)

        return [
            intera.DomainMotif(
                x[0],
                x[1],
                evidences = dmotif.evidences,
            )
            for x in itertools.product(ds, ps)
        ]


    def translate_residue(
            self,
            residue: intera.Residue,
        ) -> list[intera.Residue]:

        return [
            intera.Residue(r[3], r[2], r[0], isoform = r[1])
            for r in self.translate_site(
                residue.protein,
                residue.name,
                residue.number,
                residue.isoform,
            )
        ]


    def translate(self, x, return_strings = False, **kwargs):
        """
        Translates anything: string notation, intera objects, tuples.

        - one PTM provided as tuple of (UniProt, amino acid, offest)
        - one PTM provided as string (e.g. `P00533_S231`)
        - instance from pypath.intera: DomainMotif, Domain or Ptm

        Additional arguments can be isoform and typ (modification type).
        """

        result = []

        if type(x) is tuple:

            result = self.translate_site(*x, **kwargs)

        elif type(x) in _const.CHAR_TYPES:

            ptm = self.reptm.match(x)

            if ptm is not None:

                result = self.translate_site(
                    ptm[1],
                    ptm[2],
                    int(ptm[3]),
                    **kwargs
                )

        if return_strings:

            result = ['%s_%s%u' % (r[0], r[2], r[3]) for r in result]

        elif isinstance(x, intera.Ptm):

            result = self.translate_ptm(x)

        elif isinstance(x, intera.Domain):

            result = self.translate_domain(x)

        elif isinstance(x, intera.DomainMotif):

            result = self.translate_domain_motif(x)

        return result


    def ptm_orthology(self):
        """
        Load PTM orthology data from PhosphoSite.

        Creates an orthology translation dict of phosphosites
        based on phosphorylation sites table from PhosphoSitePlus.
        In the result all PTMs represented by a tuple of the following
        6 elements: UniProt ID, isoform (int), residue one letter code,
        residue number (int), NCBI Taxonomy ID (int), modification type.
        """

        self.ptmortho = {}

        nondigit = re.compile(r'[^\d]+')

        unknown_taxa = set()

        for typ in _data.common_load('psite_mod_types'):

            groups = {}

            url = urls.urls['psite_%s' % typ[0]]['url']
            c = curl.Curl(url, silent = False, large = True)

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

                        if site1 not in self.ptmortho:

                            self.ptmortho[site1] = {}

                        if site2[4] not in self.ptmortho[site1]:

                            self.ptmortho[site1][site2[4]] = set([])

                        self.ptmortho[site1][site2[4]].add(site2)

        if len(unknown_taxa):

            self._log(
                'Unknown taxa encountered: %s' % (
                    ', '.join(sorted(unknown_taxa))
                )
            )

    def _get_source(self, source: str | int | None) -> int:
        """
        Returns the NCBI Taxonomy ID of the source taxon.
        """

        ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(source) or self.source

        if not ncbi_tax_id:

            msg = (
                f'No source taxon provided (argument: `{source}`, '
                f'instance: `{self.source}`)'
            )
            self._log(msg)
            raise ValueError(msg)

        return ncbi_tax_id


    def __len__(self):

        return len(getattr(self, 'ptmortho', ()))


    def __repr__(self):

        return f'<PTM Orthology: {len(self)} sites>'


def init():
    """
    Initialize the orthology manager.

    Creates an instance of the orthology manager. Stores it in the module
    namespace.
    """

    globals()['manager'] = OrthologyManager()


def get_manager():
    """
    Access the orthology manager.

    Returns the orthology manager, an object which loads and unloads the
    orthology lookup tables as necessary, and provides the interface for
    querying the orthology data. Normally an instance of the manager
    belongs to the module, and if it does not exist yet, will be created
    automatically.
    """

    if 'manager' not in globals():

        init()

    return globals()['manager']


def translate(
        identifiers: str | Iterable[str],
        target: str | int,
        source: str | int = 9606,
        id_type: str = 'uniprot',
        only_swissprot: bool = True,
        oma: bool = None,
        homologene: bool = None,
        ensembl: bool = None,
        oma_rel_type: (
            set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
            None
        ) = None,
        oma_score: float | None = None,
        ensembl_hc: bool = True,
        ensembl_types: (
            list[Literal['one2one', 'one2many', 'many2many']] |
            None
        ) = None,
        full_records: bool = False,
    ):
    """
    Translate one or more identifiers by orthologous gene pairs.

    Args:
        identifiers:
            One or more identifers of the source organism, of ID type
            `id_type`.
        target:
            Name or NCBI Taxonomy ID of the target organism.
        source:
            Name or NCBI Taxonomy ID of the source organism.
        id_type:
            The identifier type to use.
        only_swissprot:
            Use only SwissProt IDs.
        oma
            Use orthology information from the Orthologous Matrix (OMA).
            Currently this is the recommended source for orthology data.
        homologene:
            Use orthology information from NCBI HomoloGene.
        ensembl:
            Use orthology information from Ensembl.
        oma_rel_type:
            Restrict relations to certain types.
        oma_score:
            Lower threshold for similarity metric.
        ensembl_hc:
            Use only the high confidence orthology relations from Ensembl.
        ensembl_types:
            Ensembl orthology relation types to use. Possible values are
            `one2one`, `one2many` and `many2many`. By default only
            `one2one` is used.
        full_records:
            Include not only the identifiers, but also some properties of
            the orthology relationships.

    Returns:
        Set of identifiers of orthologous genes or proteins in the
        target taxon.
    """

    manager = get_manager()

    args = locals().copy()
    args.pop('manager')

    return manager.translate(**args)


def get_dict(
        target: str | int,
        source: str | int = 9606,
        id_type: str = 'uniprot',
        only_swissprot: bool = True,
        oma: bool = None,
        homologene: bool = None,
        ensembl: bool = None,
        oma_rel_type: (
            set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
            None
        ) = None,
        oma_score: float | None = None,
        ensembl_hc: bool = True,
        ensembl_types: (
            list[Literal['one2one', 'one2many', 'many2many']] |
            None
        ) = None,
        full_records: bool = False,
    ) -> dict[str, set[OrthologBase]]:
    """
    Create a dictionary for one source organism and ID type.

    Args:
        target:
            Name or NCBI Taxonomy ID of the target organism.
        source:
            Name or NCBI Taxonomy ID of the source organism.
        id_type:
            The identifier type to use.
        only_swissprot:
            Use only SwissProt IDs.
        oma
            Use orthology information from the Orthologous Matrix (OMA).
            Currently this is the recommended source for orthology data.
        homologene:
            Use orthology information from NCBI HomoloGene.
        ensembl:
            Use orthology information from Ensembl.
        oma_rel_type:
            Restrict relations to certain types.
        oma_score:
            Lower threshold for similarity metric.
        ensembl_hc:
            Use only the high confidence orthology relations from Ensembl.
        ensembl_types:
            Ensembl orthology relation types to use. Possible values are
            `one2one`, `one2many` and `many2many`. By default only
            `one2one` is used.
        full_records:
            Include not only the identifiers, but also some properties of
            the orthology relationships.

    Returns:
        A dict with identifiers of the source organism as keys, and
        sets of their orthologs as values.
    """

    manager = get_manager()

    args = locals().copy()
    args.pop('manager')

    return manager.get_dict(**args)


def get_df(
        target: str | int,
        source: str | int = 9606,
        id_type: str = 'uniprot',
        only_swissprot: bool = True,
        oma: bool = None,
        homologene: bool = None,
        ensembl: bool = None,
        oma_rel_type: (
            set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
            None
        ) = None,
        oma_score: float | None = None,
        ensembl_hc: bool = True,
        ensembl_types: (
            list[Literal['one2one', 'one2many', 'many2many']] |
            None
        ) = None,
        full_records: bool = False,
        **kwargs
    ) -> pd.DataFrame:
    """
    Create a data frame for one source organism and ID type.

    Args:
        target:
            Name or NCBI Taxonomy ID of the target organism.
        source:
            Name or NCBI Taxonomy ID of the source organism.
        id_type:
            The identifier type to use.
        only_swissprot:
            Use only SwissProt IDs.
        oma
            Use orthology information from the Orthologous Matrix (OMA).
            Currently this is the recommended source for orthology data.
        homologene:
            Use orthology information from NCBI HomoloGene.
        ensembl:
            Use orthology information from Ensembl.
        oma_rel_type:
            Restrict relations to certain types.
        oma_score:
            Lower threshold for similarity metric.
        ensembl_hc:
            Use only the high confidence orthology relations from Ensembl.
        ensembl_types:
            Ensembl orthology relation types to use. Possible values are
            `one2one`, `one2many` and `many2many`. By default only
            `one2one` is used.
        full_records:
            Include not only the identifiers, but also some properties of
            the orthology relationships.
        kwargs:
            Ignored.

    Returns:
        A data frame with pairs of orthologous identifiers,
        in two columns: "source" and "target".
    """

    manager = get_manager()

    args = locals().copy()
    args.pop('manager')
    args.pop('kwargs')

    return manager.get_df(**args)


def translate_df(
        df: pd.DataFrame,
        target: str | int,
        source: str | int = 9606,
        cols: str | list[str] | dict[str, str] | None = None,
        id_type: str = 'uniprot',
        only_swissprot: bool = True,
        oma: bool = None,
        homologene: bool = None,
        ensembl: bool = None,
        oma_rel_type: (
            set[Literal['1:1', '1:n', 'm:1', 'm:n']] |
            None
        ) = None,
        oma_score: float | None = None,
        ensembl_hc: bool = True,
        ensembl_types: (
            list[Literal['one2one', 'one2many', 'many2many']] |
            None
        ) = None,
        **kwargs: str | tuple[str, str]
    ) -> pd.DataFrame:
    """
    Translate columns in a data frame.

    Args:
        df:
            A data frame.
        target:
            Name or NCBI Taxonomy ID of the target organism.
        source:
            Name or NCBI Taxonomy ID of the source organism.
        cols:
            One or more columns to be translated. It can be a single
            column name, an iterable of column names or a dict where
            keys are column names and values are ID types. Except this
            last case, identifiers are assumed to be `id_type`.
        id_type:
            The default identifier type to use, will be used for all
            columns where ID type is not specified.
        only_swissprot:
            Use only SwissProt IDs.
        oma
            Use orthology information from the Orthologous Matrix (OMA).
            Currently this is the recommended source for orthology data.
        homologene:
            Use orthology information from NCBI HomoloGene.
        ensembl:
            Use orthology information from Ensembl.
        oma_rel_type:
            Restrict relations to certain types.
        oma_score:
            Lower threshold for similarity metric.
        ensembl_hc:
            Use only the high confidence orthology relations from Ensembl.
        ensembl_types:
            Ensembl orthology relation types to use. Possible values are
            `one2one`, `one2many` and `many2many`. By default only
            `one2one` is used.
        kwargs:
            Same as providing a dict to ``cols``, but beware, keys
            (column names) can not match existing argument names of
            this function.

    Returns:
        A data frame with the same column layout as the input, and the
        identifiers translated as demanded. Rows that could not be
        translated are omitted.
    """

    manager = get_manager()

    args = locals().copy()
    args.pop('manager')
    args.pop('kwargs')

    return manager.translate_df(**args, **kwargs)
