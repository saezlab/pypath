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

import json

import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.settings as settings
import pypath.resources.urls as urls

_logger = session.Logger(name = 'disgenet_input')
_log = _logger._log


class DisgenetAuth(session.Logger):

    _name = 'DisGeNET'
    _api_url = urls.urls['disgenet']['api_url']
    _authenticated: bool = False
    _api_key: str = None

    def __init__(self, email: str | None, password: str | None):
        """
        Handle user authentication to the DisGeNet API.
        """

        if not hasattr(self, '_log_name'):

            session.Logger.__init__(self, name = 'disgenet')

        self.email = email
        self.password = password


    def authenticate(self) -> bool:
        """
        Starts an authorization process to the DisGeNET API.
        Returns a boolean which is success of authentication.
        """

        if self._authenticated and self._api_key is not None:

            return True

        self._log(f'Authorizing to the {self._name} API.')
        email: str = self.email or settings.get('disgenet_email')
        password: str = self.password or settings.get('disgenet_password')

        if not email or not password:

            self._log(
                'Email or password missing: '
                'unable to authenticate to the DisGeNet API.'
            )
            return False

        url: str = f'{self._api_url}/auth/'
        post_params: dict[str, str] = {'email': email, 'password': password}
        headers: list[str] = [
            'Accept: */*',
            'Content-Type: application/x-www-form-urlencoded',
        ]

        c = curl.Curl(url = url, post = post_params, req_headers = headers)

        if c.result:

            try:
                self._api_key = json.loads(c.result).get('token')
                self._log('DisGeNet API Authentication successful.')

            except Exception as e:

                self._log(
                    'Failed to process response from '
                    'DisGeNet API authentication:'
                )
                self._log_traceback()

        self._authenticated = self._api_key is not None

        return self._authenticated


    def _if_authenticated(f):
        """
        Decorator to ensure DisGeNet API authentication.
        """

        def wrapper(self, *args, **kwargs):

            if self.authenticate():

                return f(self, *args, **kwargs)

            else:

                err = (
                    'Unable to connect DisGeNet API in lack of authorization. '
                    'Please check your credentials.'
                )
                self._log(err)
                raise RuntimeError(err)

        return wrapper


    def _delete_cache(f):
        """
        Decorator which calls under `cache_delete_on` context.
        """

        def wrapper(*args, **kwargs):

            with curl.cache_delete_on():

                return f(*args, **kwargs)

        return wrapper
