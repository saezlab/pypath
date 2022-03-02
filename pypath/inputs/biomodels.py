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

def get_single_model_information(model_id, invalidate = False):
    """
    Get information on single BioModel using bioservices or load from
    cache.
    
    Args:
        model_id (str): ID of model
        
        invalidate (bool): Whether to invalidate the cache of the given 
        model

    Returns:
        dict: dictionary containing model specifics, eg name/id,
        description, associated files
    """

    @cache.region('long_term', 'get_single_biomodels_model_information')
    def _get_single_model_information(model_id):

        model = bm.get_model(model_id)

        return model

    if invalidate:
        cache.region_invalidate(
                _get_single_model_information, 
                None, 
                "get_single_biomodels_model_information", 
                model_id
            )
    
    _log(f"Returning single BioModel: {model_id}.")
    return _get_single_model_information(model_id)

def get_all_models(invalidate=False):
    """
    Fetch list of available models using bioservices or load from cache.

    Args:
        invalidate (bool): Whether to invalidate the cache file

    Returns: 
        dict: A dictionary of models with model identifiers as keys and
        model attributes as values. Model attributes include format
        (SMBL being most common), model id, name,
        submission/modification date and author.
    """
    @cache.region('long_term', 'get_all_biomodels_models')
    def _get_all_models():

        # CAUTION: Below function has a bug causing an infinite loop (Mar22)
        models = bm.get_all_models()

        return models

    if invalidate:
        cache.region_invalidate(
                _get_all_models, 
                None, 
                "get_all_biomodels_models"
            )

    _log("Returning all BioModels models.")
    return _get_all_models()
    

def get_single_model_main_file(model_id, filename, version, invalidate = False):
    """
    Download or load from cache the main file of a single model to
    extract relevant data for pypath integration. Parse original XML to
    OrderedDict.

    Args:
        model_id (str): Model identifier
    
        filename (str): Model filename

        version (int): Current version of model derived from single
        model information; necessary to generate the correct download
        URL

        invalidate (bool): Whether to invalidate the cache file for a
        given model and filename

    Returns:
        OrderedDict: complete model data
    """

    @cache.region('long_term', 'download_single_biomodels_model_main_file')
    def _get_single_model_main_file(model_id, filename, version):

        if not "xml" in filename.lower():
            _log("Warning: it appears you are attempting to download a non-xml " 
            "file. This function is designed to retrieve the main BioModels "
            "XML file for a given model.")
            return False

        if version > 1:
            model_id = f"{model_id}.{version}"

        params = {}
        if filename:
            params["filename"] = filename

        dl = bm.http_get(
                f"model/download/{model_id}", params=params
            )

        
        return xmltodict.parse(dl.text)
    
    if invalidate:
        cache.region_invalidate(
                _get_single_model_main_file, 
                None, 
                "download_single_biomodels_model_main_file",
                model_id,
                filename,
                version
            )

    _log(f"Returning BioModels download: model {model_id}, file {filename}.")
    return _get_single_model_main_file(model_id, filename, version)
    

def parse_biomodels(invalidate=False):
    """
    Get the entire list of BioModels, download each model's main file, 
    and parse the main file for integration data.
    """

    # @cache.region("long_term", "parse_biomodels")
    def _parse_biomodels():

        all_models = get_all_models(invalidate)

        all_model_info = [
            get_single_model_information(model_id, invalidate) 
            for model_id in
            [model["id"] for model in all_models]
        ]

        def get_model_file_version_tuples(model):
            """
            Get model id, online file name, and current version.

            Args:
                model (str): dict of model info from
                `get_single_model_information()`  

            Returns:
                tuple of str: model id, filename, and current version
            """
            id = model.get("publicationId") or model.get("submissionId")
            history = model.get("history")
            if history and model["files"].get("main"):
                version = history.get("revisions")[-1].get("version")
                filename = model["files"]["main"][0]["name"]
                return (id, filename, version)
            else:
                _log(f"Model {id} does not appear to have an associated main "
                    "file and/or version.")
                return (id, None, None)

        # parse and extract TODO
        for model in all_model_info:
            model_id, filename, version = get_model_file_version_tuples(model)
            f = get_single_model_main_file(
                model_id, 
                filename, 
                version, 
                invalidate
            )

            d = f.get("sbml").get("model")

            # get RDF description
            rdf_desc = d["annotation"]["rdf:RDF"]["rdf:Description"]
            # disease 
            rdf_desc.get("")
            
            ["bqbiol:hasProperty"][0]["rdf:Bag"]["rdf:li"]["@rdf:resource"]

            # tissue

            # compound
            

    # if invalidate:
    #     cache.region_invalidate(
    #             _parse_biomodels, 
    #             None, 
    #             "parse_biomodels",
    #             invalidate
    #         )

    return _parse_biomodels()