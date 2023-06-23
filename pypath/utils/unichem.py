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

import json
import os
import sys
import textwrap

import bs4

import pypath.share.progress as progress
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.inputs.unichem as unichem_input


class Unichem(object):
    """
    Client for the UniChem drug compound identifier translation service
    (https://www.ebi.ac.uk/unichem/).
    """

    def __init__(self):

        sys.stdout.write(
            '\n\tType `Unichem_instance.usage()` to get help.\n\n'
        )
        sys.stdout.flush()

        # from unichem id to db name
        self.uc_dict = unichem_input.unichem_sources()
        # from db name to unichem id
        self.name_dict = common.swap_dict(self.uc_dict)
        self.url_stem = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id'
        self.inchi_stem = 'https://www.ebi.ac.uk/unichem/rest/inchikey/%s'
        self.chembl_url = (
            'http://www.ebi.ac.uk/chemblws/compounds/smiles/{0}.json'
        )
        self.cpd_search = 'http://www.ebi.ac.uk/unichem/rest/{0}/{1}/{2}{3}'
        self.result = {}



    def usage(self):
        """
        Prints usage information and examples to the standard output.
        """

        msg = '''
        List of identifier types can be read above.
        To query UniChem, give either names or numbers of the
        ID types you wish to translate from and to.
        E.g.
        >>> u = unichem.Unichem()
        >>> u.translate('pubchem', 'chembl', list_of_pubchems)

        For example, the PubChem CID of Aspirin is 2244. Translate it to
        ChEMBL:
        >>> u.translate('pubchem', 'chembl', '2244')
        >>> u.result
        {'2244': ['CHEMBL25']}

        We can translate multiple identifiers the same way:
        >>>
        >>> u.translate('pubchem', 'chembl', ['2244', '4091'])
        >>> u.result
        {'2244': ['CHEMBL25'], '4091': ['CHEMBL1431']}

        Additional ways of translation are from SMILEs to ChEMBL IDs, and
        from InChiKeys to any ID. These are translated not by the UniChem,
        but by the ChEMBL webservice.
        >>> u.translate('smiles', 'chembl', list_of_smiles)
        >>> u.translate('inchikey', 'chembl', list_of_inchikeys)

        Other option to search is connectivity search from UniChem.
        A-G parameters can be defined optionally. See description at
        https://www.ebi.ac.uk/unichem/info/widesearchInfo
        >>> u.connectivity_search(list_of_zincs, 'zinc', parameters=[1,0,0,0,0,1,0])

        InChiKeys can be used in connectivity search too:
        >>> u.connectivity_search(list_of_inchikeys, 'inchikey', parameters=[1,0,0,0,0,1,0])

        You can also call directly functions accessing ChEMBL webservice, with the
        same result as you would call `translate()` or `connectivity_search()`:
        >>> u.smiles2chembl(list_of_smiles)
        >>> u.inchikey2anything('chembl', list_of_inchikeys)

        Find the dict in `u.result`. Untranslated items have value `None`.
        Every call overwrites previous result!

        For an up to date list of identifier types see
        https://www.ebi.ac.uk/unichem/ucquery/listSources or
        call `Unichem.info(<source>)`:
        >>> Unichem.info('chembl')
        '''

        sys.stdout.write(os.linesep)

        id_types = sorted(
            self.uc_dict.items(),
            key = lambda x: int(x[0])
        )

        if len(id_types) % 2:

            id_types.append(('',) * 2)

        nrows = len(id_types) // 2

        for i in xrange(nrows):

            sys.stdout.write(
                ''.join((
                    ' ' * 8,
                    id_types[i][0].rjust(2),
                    ' ' * 3,
                    id_types[i][1].ljust(20),
                    id_types[i + nrows][0].rjust(2),
                    ' ' * 3,
                    id_types[i + nrows][1].ljust(20),
                    os.linesep,
                ))
            )

        sys.stdout.write(msg + os.linesep)
        sys.stdout.flush()


    @staticmethod
    def info(source):
        """
        Print information about one source.

        Args
            source (int,str): The numeric or string ID of one source.
        """

        unichem_input.info(source)


    def translate(self, source, target, lst):
        """
        Translate one drug compound identifier to another identifier type
        using the UniChem web service. For an up to date list of identifier
        types see https://www.ebi.ac.uk/unichem/ucquery/listSources.

        Args
            source (str,int): The source ID type, either as a string label
                or as a number, as used in UniChem.
            target (str,int): The target ID type, either as a string label
                or as a number, as used in UniChem.
            lst (str,set): One or more identifiers to translate.

        Returns
            Returns None, the results are stored in the `result` attribute
            of this object.
        """

        lst = common.to_set(lst)
        self.result = {}

        if source == 'inchikey':

            self.inchikey2anything(target, lst)
            return None

        if source == 'smiles':

            self.smiles2chembl(lst)
            return None

        source = (
            str(source)
                if type(source) is int else
            self.name_dict[source]
        )
        target = (
            str(target)
                if type(target) is int else
            self.name_dict[target]
        )

        prg = progress.Progress(
            total=len(lst),
            name='Translating compound identifiers',
            interval=1,
        )

        for comp in lst:

            url = '/'.join([self.url_stem, comp, source, target])
            c = curl.Curl(url, large = False)
            result = c.result
            self.result[comp] = []

            if result is not None:

                data = json.loads(result)

                for d in data:

                    self.result[comp].append(d['src_compound_id'])

            prg.step()

        prg.terminate()


    def inchikey2anything(self, target, lst):
        """
        Translate InChi keys to another identifier type using the ChEMBL
        web service.

        Args
            target (str,int): The target ID type, either as a string label
                or as a number, as used in UniChem.
            lst (str,set): One or more InChi keys.

        Returns
            Returns None, the results are stored in the `result` attribute
            of this object.
        """

        lst = common.to_set(lst)

        self.result = {}
        target = (
            str(target)
                if type(target) is int else
            self.name_dict[target]
        )
        prg = progress.Progress(
            total=len(lst),
            name='Translating InChi-Keys',
            interval=1,
        )

        for inchik in lst:

            url = self.inchi_stem % inchik
            c = curl.Curl(url, large = False)
            result = c.result

            if result is not None:

                data = json.loads(result)
                self.result[inchik] = [
                    d['src_compound_id']
                    for d in data
                    if d['src_id'] == target
                ]
            prg.step()

        prg.terminate()


    def smiles2chembl(self, smiles):
        """
        Translate SMILES to ChEMBL ID using the ChEMBL web service.

        Args
            smiles (str,list): One or more SMILES.

        Returns
            Returns None, the results are stored in the `result` attribute
            of this object.
        """

        smiles = common.to_set(smiles)

        self.result = {}
        prg = progress.Progress(
            total=len(smiles),
            name='Translating SMILEs',
            interval=1
        )

        for sml in smiles:

            url = self.chembl_url.format(sml)
            c = curl.Curl(url, large = False)
            result = c.result
            self.result[sml] = []

            if result is not None:

                try:

                    data = json.loads(result)

                    for d in data['compounds']:

                        this_smile = d['smiles']
                        this_chembl = d['chemblId']
                        # if this_smile == sml:
                        self.result[sml].append(this_chembl)

                except ValueError:

                    soup = bs4.BeautifulSoup(result)
                    compounds = soup.find_all('compound')

                    if compounds is not None:

                        for compound in compounds:

                            this_smile = compound.find('smiles').text
                            this_chembl = compound.find('chemblid').text
                            # if this_smile == sml:
                            self.result[sml].append(this_chembl)
            prg.step()

        prg.terminate()


    def connectivity_search(
            self,
            id_list: str | set,
            id_type: str | int,
            parameters: list[int] = [1, 0, 0, 0, 0, 1, 0]
        ):
        """
        Search for structurally and chemically similar compounds based on
        cheminformatics similarity metrics. Read more at
        https://www.ebi.ac.uk/unichem/info/widesearchInfo.

        Args
            id_list:
                One or more identifiers to query.
            id_type:
                Type of the identifiers, either as a string
                label or a number as used by UniChem. SMILES is not
                available in this type of query.
            parameters:
                A list of parameters A-H as described in
                https://www.ebi.ac.uk/unichem/info/widesearchInfo.

        Returns
            Returns None, the results are stored in the `result` attribute
            of this object.
        """

        id_list = common.to_set(id_list)

        parameters.append(1)  # H parameter must be 1 to process the result
        parameters = [str(i) for i in parameters]
        self.result = {}

        if id_type == 'inchikey':

            id_type = ''
            method = 'key_search'

        elif id_type == 'smiles':

            return None

        else:

            id_type = (
                str(id_type)
                    if type(id_type) is int else
                self.name_dict[id_type]
            )
            id_type = '%s/' % id_type
            method = 'cpd_search'

        prg = progress.Progress(
            total=len(id_list),
            name='Connectivity search',
            interval=1
        )

        for i in id_list:

            prg.step()
            url = self.cpd_search.format(
                method,
                i,
                id_type,
                '/'.join(parameters)
            )
            c = curl.Curl(url, large = False)
            result = c.result
            self.result[i] = []

            if result is not None:

                data = json.loads(result)

                for k, v in iteritems(data):

                    for j in xrange(1, len(v)):

                        self.result[i].append(v[j][0])

            self.result[i] = list(set(self.result[i]))

        prg.terminate()
