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

from curses import KEY_A1
import bioservices.biomodels as biom
from sbmlutils.io import read_sbml
import xmltodict, json
import pandas as pd

import pypath.share.session as session
import pypath.share.beaker as bk

_logger = session.Logger(name = 'biomodels_input')
_log = _logger._log
cache = bk.get_cache_manager()

bm = biom.BioModels()


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
    

def get_single_model_main_file(model_id, filename, version, invalidate = False):
    """
    Download or load from cache the main file of a single model to
    extract relevant data for pypath integration. Parse original XML to
    dict.

    Args:
        model_id (str): Model identifier
    
        filename (str): Model filename

        version (int): Current version of model derived from single
        model information; necessary to generate the correct download
        URL

        invalidate (bool): Whether to invalidate the cache file for a
        given model and filename

    Returns:
        dict: complete model data
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

        # as xml (.text)
        dl = bm.http_get(
                f"model/download/{model_id}", params=params
            )

        # as dict
        d = xmltodict.parse(dl.text, dict_constructor=dict)
        # as json
        # j = json.dumps(d)

        return d
        
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
            [model["id"] for model in all_models[0:1]]
            # TODO revert to full range
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
        for model_info in all_model_info:
            model_id, filename, version = get_model_file_version_tuples(model_info)
            f = get_single_model_main_file(
                model_id, 
                filename, 
                version, 
                invalidate = True
            )

            ## read using SBMLutils directly from xml string
            # annot = read_sbml(f)
            # this does not return annotations ATM ...

            ## read from parsed dict and find relevant annotations
            """
            We are interested in bqbiol annotations
            (https://co.mbine.org/standards/qualifiers), which can refer
            to several relevant ontologies: - uniprot, - cl, - go, -
            kegg (kegg.compound, kegg.drug), - chebi, - taxonomy, -
            pubmed, - ..?

            A model has general annotation (in the model root), but
            model components (eg species) also have their own
            annotations which describe the individual component.

            We need to parse the entire dict for relevant information
            that should be connected to pypath constituents.

            The structure of annotation is always: (to be confirmed) 

            - "annotation" 

            - "rdf:RDF" 

            - "rdf:Description" 

            - specific descriptors (eg "bqbiol:is",
              "bqbiol:hasProperty") 
            
            - one "rdf:Bag" or a list containing "rdf:Bag" entities 
            
            - each bagcontains (always one?) "rdf:li" which contains the
              identifier url
                
            - the identifier url gives the single identifier and the
              ontology type in the ultimate and penultimate component
            """
            # get model
            model = f.get("sbml").get("model")
            # get RDF description
            parse = parse_single_biomodel(model)

            return parse
            

    # if invalidate:
    #     cache.region_invalidate(
    #             _parse_biomodels, 
    #             None, 
    #             "parse_biomodels"
    #         )

    return _parse_biomodels()

def deconstruct_url(url: str):
    """
    Get last and second to last part of url (id and ontology).
    """
    l = url.split("/")
    return {l[-2]: l[-1]}

def parse_description(component: dict):
    try:
        desc = component.get("annotation").get("rdf:RDF").get("rdf:Description")
        res = {}
        for k, v in desc.items():
            if "bqbiol" in k:
                # single bag
                if isinstance(v, dict):
                    id_on = deconstruct_url(
                        v.get("rdf:Bag").get("rdf:li").get("@rdf:resource")
                    )
                    res.update({k: id_on})
                elif isinstance(v, list):
                    ids = [
                        i.get("rdf:Bag").get("rdf:li").get("@rdf:resource")
                        for i in v
                    ]
                    id_l = [
                        deconstruct_url(id)
                        for id in ids
                    ]
                    # id_d = {k: v for d in id_l for k, v in d.items()}
                    # same key can occur more than once; needs to remain list
                    res.update({k: id_l})
        
        return res
    except AttributeError:
        return None

def parse_single_biomodel(model):
    res = {}
    # root annotation
    ra = parse_description(model)

    res.update({"main": ra})

    # single components
    for k, v in model.items():
        if "listOf" in k:
            # TODO generalise? which lists are relevant?
            if "Species" in k:
                species = v.get("species")
                if isinstance(species, list):
                    rc_d = {}
                    for i, s in enumerate(species):
                        d = parse_description(s)
                        if d is not None:
                            rc_d.update({i: d})
                    res.update({"species": rc_d})

    return res