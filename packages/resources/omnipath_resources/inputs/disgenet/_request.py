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

from typing import Any, Iterable, Literal, Sequence

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.inputs.disgenet._auth as _auth
import pypath.inputs.disgenet._valid as _valid


class DisgenetRequest(_auth.DisgenetAuth):

    _query_param_defaults = {
        'format': 'json',
    }

    def __init__(
            self,
            query_type: _valid.Querytype | str,
            by: _valid.By | str,
            id_type: _valid.Idtype | str | None,
            identifiers: str | Iterable[str] | None,
            evidences: bool = False,
            query_param: dict[str, Any] | None = None,
            email: str | None = None,
            password: str | None = None,
            fail_on_invalid: bool = True,
            **kwargs
        ):

        super().__init__(email = email, password = password)
        query_param = {**(query_param or {}), **kwargs}
        query_param.update(self._query_param_defaults)
        query_param.pop('__class__')
        del kwargs, email, password
        self._param = locals()
        self._param.pop('self')
        self._param.pop('__class__')


    def __getitem__(self, key: str) -> Any:

        return self._param.get(key)


    @property
    def url(self) -> str:

        if not hasattr(self, '_url'):

            self.prepare()

        return self._url


    def prepare(self) -> None:

        self.validate()
        self._root()
        self._query()
        self._evidences()
        self._by()
        self._disease()
        self._id_type()
        self._identifiers()


    def validate(self) -> bool:
        """
        Validity checks on a DisGeNet query.

        Note: This is not a complete check, does not guarantee the query will
        be successful.

        Return:
            True if the query is valid, False if invalid but raising
            exception is disabled.

        Raises:
            ValueError:
                If any of the validations failed.
        """

        issues = []

        for argname, values in self._param.items():

            valids_name = argname.replace('_', '').lower().capitalize()
            valids = getattr(_valid, valids_name, None)

            if not valids:

                continue

            invalids = [
                val
                for val in common.to_list(values)
                if not self._is_valid(val, valids)
            ]

            if invalids:

                msg = f'Invalid value(s) for {argname}: {", ".join(invalids)}.'
                self._log(msg)
                issues.append(msg)

        for k in self._query_param.keys():

            if not self._is_valid(k, _valid.Get):

                msg = f'Invalid query string (HTTP GET) parameter: `{k}`.'
                self._log(msg)
                issues.append(msg)

        if issues and self['fail_on_invalid']:

            issues = ' '.join(issues)
            raise ValueError(f'Failed to process DisGeNet query: {issues}')

        return bool(issues)


    @staticmethod
    def _is_valid(value: str, valids: _valids.enum.StrEnum | None) -> bool:

        return valids is None or value in valids._value2member_map_


    def _root(self) -> None:

        self._url = self._api_url


    def _append(self, part: str | None) -> None:

        if part:

            self._url = self._url.strip('/') + '/' + part.strip('/')


    def _query(self) -> None:

        self._append(self['query_type'])


    def _evidences(self) -> None:

        if self['evidences']:

            self._append('evidences')


    def _by(self) -> None:

        plural = 's' if self._dda else ''
        self._append(self['by'] + plural)

    def _disease(self) -> None:

        if self._dda:

            self._append('disease')


    def _id_type(self) -> None:

        self._append(self['id_type'])


    def _identifiers(self) -> None:

        if self['identifiers']:

            ids = common.to_list(self['identifiers'])
            limit = self._query_param.get('limit', None)

            if limit and len(ids) > limit:

                ids = ids[:limit]
                self._log(f'Limiting to {limit} identifiers.')

            self._append(','.join(ids))


    @property
    def _dda(self) -> bool:

        return self['query_type'] == 'dda'


    @property
    def _query_param(self) -> dict[str, str]:

        return self['query_param'] or {}


    @property
    def get(self) -> dict[str, str]:

        return {
            k: ','.join(map(str, common.to_list(v)))
            for k, v in self._query_param.items()
        }


    @_auth.DisgenetAuth._if_authenticated
    def retrieve(self, url: str | None) -> curl.Curl:

        headers = ['Accept: */*', f'Authorization: Bearer {self._api_key}']
        url = url or self.url

        return curl.Curl(
            url,
            get = self.get,
            large = True,
            req_headers = headers,
        )


    def __iter__(self) -> Generator[dict]:

        url = self.url

        while True:

            if not url:

                break

            c = self.retrieve(url)
            url = None

            if c.status != 200 or not hasattr(c, 'fileobj'):

                msg = (
                    'Transaction to DisGeNet API failed, '
                    'check the log for details.'
                )
                self._log(msg)
                raise RuntimeError(msg)

            result = json.load(c.fileobj)

            if isinstance(result, dict):

                url = result.get('url', None)
                result = result.get('results', result)

            for record in result:

                yield record
