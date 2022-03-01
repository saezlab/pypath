#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
This file collects all available models from the biomodels repository
using the bioservices python module, then downloads the individual
models to parse for the relevant information to enter into pypath.
"""

import bioservices.biomodels as biom
import xmltodict

import pypath.share.session as session
import pypath.share.beaker as bk

_logger = session.Logger(name = 'biomodels_input')
_log = _logger._log
cache = bk.get_cache_manager()

bm = biom.BioModels()

@cache.region('long_term', 'get_single_biomodels_model_information')
def get_single_model_information(model_id):
    """
    Get information on single BioModel using bioservices or load from
    cache.
    
    Args:
        model_id (str): ID of model

    Returns:
        dict: dictionary containing model specifics, eg name/id,
        description, associated files
    """

    model = bm.get_model(model_id)

    _log(f"Returning single BioModel: {model_id}.")
    return model

@cache.region('long_term', 'get_all_biomodels_models')
def get_all_models():
    """
    Fetch list of available models using bioservices or load from cache.

    Returns: 
        dict: A dictionary of models with model identifiers as keys and
        model attributes as values. Model attributes include format
        (SMBL being most common), model id, name,
        submission/modification date and author.
    """

    # CAUTION: Below function has a bug causing an infinite loop (Mar22)
    models = bm.get_all_models()

    _log("Returning all BioModels models.")
    return models

@cache.region('long_term', 'download_single_biomodels_model_main_file')
def get_single_model_main_file(model_id, filename):
    """
    Download or load from cache the main file of a single model to
    extract relevant data for pypath integration. Parse original XML to 
    OrderedDict.

    Args:
        model_id (str): Model identifier
    
        filename (str): Model filename

    Returns:
        OrderedDict: complete model data

    """

    if not "xml" in filename.lower():
        _log("Warning: it appears you are attempting to download a non-xml " 
        "file. This function is designed to retrieve the main BioModels "
        "XML file for a given model.")
        return False

    params = {}
    if filename:
        params["filename"] = filename

    dl = bm.http_get(
            f"model/download/{model_id}", params=params
        )

    _log(f"Returning BioModels download: model {model_id}, file {filename}.")
    return xmltodict.parse(dl.text)
    

