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
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections
import csv
import json
from getpass import getpass

from typing import List, Union, NamedTuple, Dict, Tuple, Optional

import pypath.share.curl as curl
import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping

_logger = session.Logger(name='disgenet_input')
_log = _logger._log


class DisgenetApi:

    _name = 'DisGeNET'
    _api_url = urls.urls['disgenet']['api_url']
    _authenticated: bool = False
    _api_key: str = None

    def authenticate(self) -> bool:
        '''
        Starts an authorization process in DisGeNET API.
        Returns a boolean which is success of authentication.
        '''

        if self._authenticated and self._api_key != None:

            return True

        print(f'Authorizing in {self._name} API...')
        e_mail: str = input('E-mail: ')
        password: str = getpass('Password: ')

        url: str = f'{self._api_url}/auth/'
        post_params: dict[str, str] = {'email': e_mail, 'password': password}
        headers: dict[str, str] = {
            'accept: */*',
            'Content-Type: application/x-www-form-urlencoded',
        }

        c = curl.Curl(url=url, post=post_params, req_headers=headers)
        response: int = c.status

        if response == 200 or response == 0:

            result: dict[str, str] = json.loads(c.result)
            self._authenticated = True
            self._api_key = result['token']
            print('Authorization successful.')

        elif response == 404:

            self._authenticated = False
            self._api_key = None

        return self._authenticated and (self._api_key != None)

    def _if_authenticated(f):
        '''
        Simple wrapper to get rid of the burden of
        checking authentication status and authenticating
        if already haven't.
        '''

        def wrapper(self, *args, **kwargs):

            if self.authenticate():

                return f(self, *args, **kwargs)

            else:

                print('Failure in authorization, check your credentials.')

        return wrapper

    def _delete_cache(f):
        '''
        A necessary wrapper as databases may be updated
        after the initial download of a particular data.
        '''

        def wrapper(*args, **kwargs):

            with curl.cache_delete_on():

                return f(*args, **kwargs)

        return wrapper

    def get_ddas_that_share_genes(
        self,
        disease: Union[str, List[str]],
        vocabulary: str = None,
        source: str = None,
        p_value: float = None,
        limit: int = 10,
    ) -> NamedTuple(
        'DiseaseDiseaseAssociation',
        [
            ('disease1_name', str),
            ('disease2_name', str),
            ('disease1_ngenes', int),
            ('disease2_ngenes', int),
            ('disease1_disease_class', Tuple[str]),
            ('disease2_disease_class', Tuple[str]),
            ('disease1_disease_class_name', Tuple[str]),
            ('disease2_disease_class_name', Tuple[str]),
            ('jaccard_genes', float),
            ('pvalue_jaccard_genes', str),
            ('source', str),
            ('ngenes1', int),
            ('ngenes2', int),
            ('ngenes', int),
            ('nvariants1', int),
            ('nvariants2', int),
            ('diseaseid1', str),
            ('diseaseid2', str),
        ],
    ):
        '''
        Returns Disease-Disease Associations with given query.

        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10,MeSH, OMIM, DO, EFO,
                NCI, HPO, MONDO, or ORDO identifier) or list of disease ids up to 100.
            else:
                Disease id (UMLS CUI) or list of disease ids up to 100.
        @vocabulary: str
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo,
            nci, hpo, mondo, ordo
        @source: str
            Source of the DDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE,
            CGI, CLINGEN, CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND,
            GWASCAT, GWASDB, HPO, LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @p_value: float
            p value associated to the Jaccard Index based on the shared genes.
        @limit: int
            Number of associated diseases to retrieve.
            Default value : 10
        '''

        return self._get_ddas(
            disease=disease,
            share='genes',
            vocabulary=vocabulary,
            source=source,
            p_value=p_value,
            limit=limit,
        )

    def get_ddas_that_share_variants(
        self,
        disease: Union[str, List[str]],
        vocabulary: str = None,
        source: str = None,
        p_value: float = None,
        limit: int = 10,
    ) -> NamedTuple(
        'DiseaseDiseaseAssociation',
        [
            ('disease1_name', str),
            ('disease2_name', str),
            ('disease1_nvariants', int),
            ('disease2_nvariants', int),
            ('disease1_disease_class', Tuple[str]),
            ('disease2_disease_class', Tuple[str]),
            ('disease1_disease_class_name', Tuple[str]),
            ('disease2_disease_class_name', Tuple[str]),
            ('jaccard_variants', float),
            ('pvalue_jaccard_variants', str),
            ('source', str),
            ('ngenes1', int),
            ('ngenes2', int),
            ('nvariants', int),
            ('nvariants1', int),
            ('nvariants2', int),
            ('diseaseid1', str),
            ('diseaseid2', str),
        ],
    ):
        '''
        Returns Disease-Disease Associations with given query.

        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10,MeSH, OMIM, DO, EFO,
                NCI, HPO, MONDO, or ORDO identifier) or list of disease ids up to 100.
            else:
                Disease id (UMLS CUI) or list of disease ids up to 100.
        @vocabulary: str
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo,
            nci, hpo, mondo, ordo
        @source: str
            Source of the DDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE,
            CGI, CLINGEN, CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND,
            GWASCAT, GWASDB, HPO, LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @p_value: float
            p value associated to the Jaccard Index based on the shared genes.
        @limit: int
            Number of associated diseases to retrieve.
            Default value : 10
        '''

        return self._get_ddas(
            disease=disease,
            share='variants',
            vocabulary=vocabulary,
            source=source,
            p_value=p_value,
            limit=limit,
        )

    def get_vdas_by_variants(
        self,
        variant: Union[str, List[str]],
        gene: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'VariantDiseaseAssociation',
        [
            ('variantid', str),
            ('gene_symbol', str),
            ('variant_dsi', float),
            ('variant_dpi', float),
            ('variant_consequence_type', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ): 
        '''
        Returns Variant-Disease Associations by variant(s).

        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants.
        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids to filter the results.
            else:
                Disease id (UMLS CUI) or list of diseases to filter the results.
        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
        @source: str
            Source of the VDA.
            Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
        @min_score: float
            Min value of the variant-disease score range.
        @max_score: float
            Max value of the variant-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the variant.
        @max_dsi: float
            Max value of the DSI range for the variant.
        @min_dpi: float
            Min value of the DPI range for the variant.
        @max_dpi: float
            Max value of the DPI range for the variant.
        @limit: int
            Number of VDAs to retrieve.
        '''

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID')
        variant = self._list_to_str(variant, 'Variant ID', limit=100)

        return self._get_vdas(
            gene=gene,
            disease=disease,
            variant=variant,
            vocabulary=None,
            by='variant',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            limit=limit,
        )

    def get_vdas_by_genes(
        self,
        gene: Union[str, List[str]],
        disease: Union[str, List[str]] = None,
        variant: Union[str, List[str]] = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'VariantDiseaseAssociation',
        [
            ('variantid', str),
            ('gene_symbol', str),
            ('variant_dsi', float),
            ('variant_dpi', float),
            ('variant_consequence_type', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Variant-Disease Associations by gene(s).

        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes.
        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants to filter the results.
        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids to filter the results.
            else:
                Disease id (UMLS CUI) or list of diseases to filter the results.
        @source: str
            Source of the VDA.
            Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
        @min_score: float
            Min value of the variant-disease score range.
        @max_score: float
            Max value of the variant-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the variant.
        @max_dsi: float
            Max value of the DSI range for the variant.
        @min_dpi: float
            Min value of the DPI range for the variant.
        @max_dpi: float
            Max value of the DPI range for the variant.
        @limit: int
            Number of VDAs to retrieve.
        '''

        gene = self._list_to_str(gene, 'Gene ID', limit=100)
        disease = self._list_to_str(disease, 'Disease ID')
        variant = self._list_to_str(variant, 'Variant ID')

        return self._get_vdas(
            gene=gene,
            disease=disease,
            variant=variant,
            vocabulary=None,
            by='gene',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            limit=limit,
        )

    def get_vdas_by_diseases(
        self,
        disease: Union[str, List[str]],
        gene: Union[str, List[str]] = None,
        variant: Union[str, List[str]] = None,
        vocabulary: str = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'VariantDiseaseAssociation',
        [
            ('variantid', str),
            ('gene_symbol', str),
            ('variant_dsi', float),
            ('variant_dpi', float),
            ('variant_consequence_type', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Variant-Disease Associations disease(s).

        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids up to 100.
            else:
                Disease id (UMLS CUI) or list of diseases up to 100.
        @vocabulary: str
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants to filter the results.
        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
        @source: str
            Source of the VDA.
            Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
        @min_score: float
            Min value of the variant-disease score range.
        @max_score: float
            Max value of the variant-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the variant.
        @max_dsi: float
            Max value of the DSI range for the variant.
        @min_dpi: float
            Min value of the DPI range for the variant.
        @max_dpi: float
            Max value of the DPI range for the variant.
        @limit: int
            Number of VDAs to retrieve.
        '''

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID', limit=100)
        variant = self._list_to_str(variant, 'Variant ID')

        return self._get_vdas(
            gene=gene,
            disease=disease,
            variant=variant,
            vocabulary=vocabulary,
            by='disease',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            limit=limit,
        )

    def get_vdas_by_source(
        self,
        source: str,
        gene: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        variant: Union[str, List[str]] = None,
        vocabulary: str = None,
        by: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'VariantDiseaseAssociation',
        [
            ('variantid', str),
            ('gene_symbol', str),
            ('variant_dsi', float),
            ('variant_dpi', float),
            ('variant_consequence_type', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Variant-Disease Associations by source.

        @source: str
            Source of the VDA.
            Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
        @disease: Union[str, List[str]]
            Disease id (UMLS CUI) or list of diseases to filter the results..
        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants to filter the results..
        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results..
        @min_score: float
            Min value of the variant-disease score range.
        @max_score: float
            Max value of the variant-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the variant.
        @max_dsi: float
            Max value of the DSI range for the variant.
        @min_dpi: float
            Min value of the DPI range for the variant.
        @max_dpi: float
            Max value of the DPI range for the variant.
        @limit: int
            Number of VDAs to retrieve.
        '''


        disease = self._list_to_str(disease, 'Disease ID')
        variant = self._list_to_str(variant, 'Variant ID')
        gene = self._list_to_str(gene, 'Gene ID')

        return self._get_vdas(
            gene=gene,
            disease=disease,
            variant=variant,
            vocabulary=vocabulary,
            by='source',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            limit=limit,
        )

    def get_gdas_by_genes(
        self,
        gene: Union[str, List[str]],
        disease: Union[str, List[str]] = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        min_pli: float = None,
        max_pli: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'GeneDiseaseAssociation',
        [
            ('geneid', int),
            ('gene_symbol', str),
            ('uniprotid', str),
            ('gene_dsi', float),
            ('gene_dpi', float),
            ('gene_pli', float),
            ('protein_class', str),
            ('protein_class_name', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('el', str),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Gene-Disease Associations by gene(s).

        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes up to 100.
        @disease: Union[str, List[str]]
            Disease id (UMLS CUI) or list of disease ids up to 100.
        @source: str
            Source of the GDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
            CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
            LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @min_score: float
            Min value of the gene-disease score range.
        @max_score: float
            Max value of the gene-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the gene.
        @max_dsi: float
            Max value of the DSI range for the gene.
        @min_dpi: float
            Min value of the DPI range for the gene.
        @max_dpi: float
            Max value of the DPI range for the gene.
        @min_pli: float
            Min value of the pLI range.
        @max_pli: float
            Max value of the pLI range.
        @limit: int
            Number of GDAs to retrieve.
        '''

        gene = self._list_to_str(gene, 'Gene ID', limit=100)
        disease = self._list_to_str(disease, 'Disease ID', limit=100)

        return self._get_gdas(
            gene=gene,
            disease=disease,
            uniprot=None,
            vocabulary=None,
            by='gene',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            min_pli=min_pli,
            max_pli=max_pli,
            limit=limit,
        )

    def get_gdas_by_diseases(
        self,
        disease: Union[str, List[str]],
        gene: Union[str, List[str]] = None,
        vocabulary: str = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        min_pli: float = None,
        max_pli: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'GeneDiseaseAssociation',
        [
            ('geneid', int),
            ('gene_symbol', str),
            ('uniprotid', str),
            ('gene_dsi', float),
            ('gene_dpi', float),
            ('gene_pli', float),
            ('protein_class', str),
            ('protein_class_name', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('el', str),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Gene-Disease Associations by disease(s).

        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids up to 100.
            else:
                Disease id (UMLS CUI) or list of diseases up to 100.
        @vocabulary: str
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
        @by: str
            Return associations by:
            Avaliable values : genes, disease, uniprot, source
        @source: str
            Source of the GDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
            CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
            LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @min_score: float
            Min value of the gene-disease score range.
        @max_score: float
            Max value of the gene-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the gene.
        @max_dsi: float
            Max value of the DSI range for the gene.
        @min_dpi: float
            Min value of the DPI range for the gene.
        @max_dpi: float
            Max value of the DPI range for the gene.
        @min_pli: float
            Min value of the pLI range.
        @max_pli: float
            Max value of the pLI range.
        @limit: int
            Number of GDAs to retrieve.
        '''

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID', limit=100)

        return self._get_gdas(
            gene=gene,
            disease=disease,
            uniprot=None,
            vocabulary=vocabulary,
            by='disease',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            min_pli=min_pli,
            max_pli=max_pli,
            limit=limit,
        )

    def get_gdas_by_uniprots(
        self,
        uniprot: Union[str, List[str]],
        disease: Union[str, List[str]] = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        min_pli: float = None,
        max_pli: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'GeneDiseaseAssociation',
        [
            ('geneid', int),
            ('gene_symbol', str),
            ('uniprotid', str),
            ('gene_dsi', float),
            ('gene_dpi', float),
            ('gene_pli', float),
            ('protein_class', str),
            ('protein_class_name', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('el', str),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Gene-Disease Associations by UniProt Accession(s).

        @uniprot: Union[str, List[str]]
            UniProt identifier or list of UniProt identifiers up to 100.
        @disease: Union[str, List[str]]
            Disease id (UMLS CUI) or list of diseases to filter the results.
        @source: str
            Source of the GDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
            CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
            LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @min_score: float
            Min value of the gene-disease score range.
        @max_score: float
            Max value of the gene-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the gene.
        @max_dsi: float
            Max value of the DSI range for the gene.
        @min_dpi: float
            Min value of the DPI range for the gene.
        @max_dpi: float
            Max value of the DPI range for the gene.
        @min_pli: float
            Min value of the pLI range.
        @max_pli: float
            Max value of the pLI range.
        @limit: int
            Number of GDAs to retrieve.
        '''

        uniprot = self._list_to_str(uniprot, 'Uniprot ID', limit=100)
        disease = self._list_to_str(disease, 'Disease ID')

        return self._get_gdas(
            gene=None,
            disease=disease,
            uniprot=uniprot,
            vocabulary=None,
            by='uniprot',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            min_pli=min_pli,
            max_pli=max_pli,
            limit=limit,
        )

    def get_gdas_by_source(
        self,
        source: str,
        gene: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        min_pli: float = None,
        max_pli: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'GeneDiseaseAssociation',
        [
            ('geneid', int),
            ('gene_symbol', str),
            ('uniprotid', str),
            ('gene_dsi', float),
            ('gene_dpi', float),
            ('gene_pli', float),
            ('protein_class', str),
            ('protein_class_name', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('el', str),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Gene-Disease Associations by source.

        @source: str
            Source of the GDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
            CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
            LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
        @disease: Union[str, List[str]]
            Disease id (UMLS CUI) or list of diseases to filter the results.
        @min_score: float
            Min value of the gene-disease score range.
        @max_score: float
            Max value of the gene-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the gene.
        @max_dsi: float
            Max value of the DSI range for the gene.
        @min_dpi: float
            Min value of the DPI range for the gene.
        @max_dpi: float
            Max value of the DPI range for the gene.
        @min_pli: float
            Min value of the pLI range.
        @max_pli: float
            Max value of the pLI range.
        @limit: int
            Number of GDAs to retrieve.
        '''

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID')

        return self._get_gdas(
            gene=gene,
            disease=disease,
            uniprot=None,
            vocabulary=None,
            by='source',
            source=source,
            min_score=min_score,
            max_score=max_score,
            min_ei=min_ei,
            max_ei=max_ei,
            disease_type=disease_type,
            disease_class=disease_class,
            min_dsi=min_dsi,
            max_dsi=max_dsi,
            min_dpi=min_dpi,
            max_dpi=max_dpi,
            min_pli=min_pli,
            max_pli=max_pli,
            limit=limit,
        )

    def _get_ddas(
        self,
        disease: Union[str, List[str]],
        share: str = None,
        vocabulary: str = None,
        source: str = None,
        p_value: float = None,
        limit: int = 10,
    ) -> NamedTuple(
        'DiseaseDiseaseAssociation',
        [
            ('disease1_name', str),
            ('disease2_name', str),
            ('disease1_n{share}', int),
            ('disease2_n{share}', int),
            ('disease1_disease_class', Tuple[str]),
            ('disease2_disease_class', Tuple[str]),
            ('disease1_disease_class_name', Tuple[str]),
            ('disease2_disease_class_name', Tuple[str]),
            ('jaccard_{share}', float),
            ('pvalue_jaccard_{share}', str),
            ('source', str),
            ('ngenes1', int),
            ('ngenes2', int),
            ('n{share}', int),
            ('nvariants1', int),
            ('nvariants2', int),
            ('diseaseid1', str),
            ('diseaseid2', str),
        ],
    ):
        '''
        Returns Disease-Disease Associations with given query.

        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10,MeSH, OMIM, DO, EFO,
                NCI, HPO, MONDO, or ORDO identifier) or list of disease ids
            else:
                Disease id (UMLS CUI) or list of disease ids
        @share: str
            Return associations that share:
            Avaliable values : genes, variants
        @vocabulary: str
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo,
            nci, hpo, mondo, ordo
        @source: str
            Source of the DDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE,
            CGI, CLINGEN, CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND,
            GWASCAT, GWASDB, HPO, LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @p_value: float
            p value associated to the Jaccard Index based on the shared genes.
        @limit: int
            Number of associated diseases to retrieve.
            Default value : 10
        '''

        disease = self._list_to_str(disease, 'disease ID', limit=100)

        url = f'{self._api_url}/dda/{share}/disease/'

        if vocabulary != None:

            url += f'{vocabulary}/'

        url += disease
        get_params = dict()
        # headers = dict()

        if source != None:

            get_params['source'] = source

        if p_value != None:

            get_params['pvalue'] = str(p_value)

        get_params['limit'] = str(max(1, min(limit, 100)))

        result = self._retrieve_data(url, get_params)

        if result == None:

            return None

        DiseaseDiseaseAssociation = collections.namedtuple(
            'DiseaseDiseaseAssociation',
            [
                'disease1_name',
                'disease2_name',
                f'disease1_n{share}',
                f'disease2_n{share}',
                'disease1_disease_class',
                'disease2_disease_class',
                'disease1_disease_class_name',
                'disease2_disease_class_name',
                f'jaccard_{share}',
                f'pvalue_jaccard_{share}',
                'source',
                'ngenes1',
                'ngenes2',
                f'n{share}',
                'nvariants1',
                'nvariants2',
                'diseaseid1',
                'diseaseid2',
            ],
        )

        for index, entry in enumerate(result):

            result[index] = DiseaseDiseaseAssociation(
                self._get_string(entry['disease1_name']),
                self._get_string(entry['disease2_name']),
                self._get_int(entry[f'disease1_n{share}']),
                self._get_int(entry[f'disease2_n{share}']),
                self._get_tuple(entry['disease1_disease_class'], ';'),
                self._get_tuple(entry['disease2_disease_class'], ';'),
                self._get_tuple(entry['disease1_disease_class_name'], ';'),
                self._get_tuple(entry['disease2_disease_class_name'], ';'),
                self._get_float(entry[f'jaccard_{share}']),
                self._get_string(entry[f'pvalue_jaccard_{share}']),
                self._get_string(entry['source']),
                self._get_int(entry['ngenes1']),
                self._get_int(entry['ngenes2']),
                self._get_int(entry[f'n{share}']),
                self._get_int(entry['nvariants1']),
                self._get_int(entry['nvariants2']),
                self._get_string(entry['diseaseid1']),
                self._get_string(entry['diseaseid2']),
            )

        return result

    def _get_vdas(
        self,
        gene: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        variant: Union[str, List[str]] = None,
        vocabulary: str = None,
        by: str = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'VariantDiseaseAssociation',
        [
            ('variantid', str),
            ('gene_symbol', str),
            ('variant_dsi', float),
            ('variant_dpi', float),
            ('variant_consequence_type', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Variant-Disease Associations with given query.

        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes.
        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids.
            else:
                Disease id (UMLS CUI) or list of diseases.
        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants.
        @vocabulary: str
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
        @by: str
            Return associations by:
            Avaliable values : gene, disease, variant, source
        @source: str
            Source of the VDA.
            Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
        @min_score: float
            Min value of the variant-disease score range.
        @max_score: float
            Max value of the variant-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the variant.
        @max_dsi: float
            Max value of the DSI range for the variant.
        @min_dpi: float
            Min value of the DPI range for the variant.
        @max_dpi: float
            Max value of the DPI range for the variant.
        @limit: int
            Number of VDAs to retrieve.
        '''

        url = f'{self._api_url}/vda/'
        get_params = dict()

        if by == 'gene' and gene != None:

            url += f'gene/{gene}'

            if disease != None:

                get_params['disease'] = disease

            if variant != None:

                get_params['variant'] = variant

        elif by == 'disease' and disease != None:

            if vocabulary != None:

                url += f'disease/{vocabulary}/{disease}'

            else:

                url += f'disease/{disease}'

            if gene != None:

                get_params['gene'] = gene

            if variant != None:

                get_params['variant'] = variant

        elif by == 'variant' and variant != None:

            url += f'variant/{variant}'

            if disease != None:

                get_params['disease'] = disease

            if gene != None:

                get_params['gene'] = gene

        elif by == 'source':

            url += f'source/{source}'

            if disease != None:

                get_params['disease'] = disease

            if variant != None:

                get_params['variant'] = variant

            if gene != None:

                get_params['gene'] = gene
        else:

            print('Problem in function call. Check arguments.')

            return None

        if source != None:

            get_params['source'] = source

        if min_score != None:

            get_params['min_score'] = str(min_score)

        if max_score != None:

            get_params['max_score'] = str(max_score)

        if min_ei != None:

            get_params['min_ei'] = str(min_ei)

        if max_ei != None:

            get_params['max_ei'] = str(max_ei)

        if min_score != None:

            get_params['min_score'] = str(min_score)

        if disease_type != None:

            get_params['type'] = disease_type

        if disease_class != None:

            get_params['disease_class'] = disease_class

        if min_dsi != None:

            get_params['min_dsi'] = str(min_dsi)

        if max_dsi != None:

            get_params['max_dsi'] = str(max_dsi)

        if min_dpi != None:

            get_params['min_dpi'] = str(min_dpi)

        if max_dpi != None:

            get_params['max_dpi'] = str(max_dpi)

        if limit != None:

            get_params['limit'] = str(max(1, limit))

        result = self._retrieve_data(url, get_params)

        if result == None:

            return None

        VariantDiseaseAssociation = collections.namedtuple(
            'VariantDiseaseAssociation',
            [
                'variantid',
                'gene_symbol',
                'variant_dsi',
                'variant_dpi',
                'variant_consequence_type',
                'diseaseid',
                'disease_name',
                'disease_class',
                'disease_class_name',
                'disease_type',
                'disease_semantic_type',
                'score',
                'ei',
                'year_initial',
                'year_final',
                'source',
            ],
        )

        for index, entry in enumerate(result):

            result[index] = VariantDiseaseAssociation(
                self._get_string(entry['variantid']),
                self._get_string(entry['gene_symbol']),
                self._get_float(entry['variant_dsi']),
                self._get_float(entry['variant_dpi']),
                self._get_string(entry['variant_consequence_type']),
                self._get_string(entry['diseaseid']),
                self._get_string(entry['disease_name']),
                self._get_tuple(entry['disease_class'], ';'),
                self._get_tuple(entry['disease_class_name'], ';'),
                self._get_string(entry['disease_type']),
                self._get_string(entry['disease_semantic_type']),
                self._get_float(entry['score']),
                self._get_float(entry['ei']),
                self._get_int(entry['year_initial']),
                self._get_int(entry['year_final']),
                self._get_string(entry['source']),
            )

        return result

    def _get_gdas(
        self,
        gene: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        uniprot: Union[str, List[str]] = None,
        vocabulary: str = None,
        by: str = None,
        source: str = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        disease_type: str = None,
        disease_class: Union[str, List[str]] = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        min_pli: float = None,
        max_pli: float = None,
        limit: int = None,
    ) -> NamedTuple(
        'GeneDiseaseAssociation',
        [
            ('geneid', int),
            ('gene_symbol', str),
            ('uniprotid', str),
            ('gene_dsi', float),
            ('gene_dpi', float),
            ('gene_pli', float),
            ('protein_class', str),
            ('protein_class_name', str),
            ('diseaseid', str),
            ('disease_name', str),
            ('disease_class', Tuple[str]),
            ('disease_class_name', Tuple[str]),
            ('disease_type', str),
            ('disease_semantic_type', str),
            ('score', float),
            ('ei', float),
            ('el', str),
            ('year_initial', int),
            ('year_final', int),
            ('source', str),
        ],
    ):
        '''
        Returns Gene-Disease Associations with given query.

        @gene: Union[str, List[str]]
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes.
        @disease: Union[str, List[str]]
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids.
            else:
                Disease id (UMLS CUI) or list of diseases.
        @uniprot: Union[str, List[str]]
            Disease id (UMLS CUI) or list of disease ids.
        @vocabulary: str
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
        @by: str
            Return associations by:
            Avaliable values : gene, disease, uniprot, source
        @source: str
            Source of the GDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
            CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
            LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        @min_score: float
            Min value of the gene-disease score range.
        @max_score: float
            Max value of the gene-disease score range.
        @min_ei: float
            Min value of the evidence index range.
        @max_ei: float
            Max value of the evidence index range.
        @disease_type: str
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        @disease_class: Union[str, List[str]]
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        @min_dsi: float
            Min value of the DSI range for the gene.
        @max_dsi: float
            Max value of the DSI range for the gene.
        @min_dpi: float
            Min value of the DPI range for the gene.
        @max_dpi: float
            Max value of the DPI range for the gene.
        @min_pli: float
            Min value of the pLI range.
        @max_pli: float
            Max value of the pLI range.
        @limit: int
            Number of GDAs to retrieve.
        '''

        url = f'{self._api_url}/gda/'
        get_params = dict()

        if by == 'gene' and gene != None:

            url += f'gene/{gene}'

            if disease != None:

                get_params['disease'] = disease

        elif by == 'disease' and disease != None:

            if vocabulary != None:

                url += f'disease/{vocabulary}/{disease}'

            else:

                url += f'disease/{disease}'

            if gene != None:

                get_params['gene'] = gene

        elif by == 'uniprot' and uniprot != None:

            url += f'gene/uniprot/{uniprot}'

            if disease != None:

                get_params['disease'] = disease

        elif by == 'source' and source != None:

            url += f'source/{source}'

            if gene != None:

                get_params['gene'] = gene

            if disease != None:

                get_params['disease'] = disease

        else:

            print('Problem in function call. Check arguments.')

            return None

        if by != 'source' and source != None:

            get_params['source'] = source

        if min_score != None:

            get_params['min_score'] = str(min_score)

        if max_score != None:

            get_params['max_score'] = str(max_score)

        if min_ei != None:

            get_params['min_ei'] = str(min_ei)

        if max_ei != None:

            get_params['max_ei'] = str(max_ei)

        if min_score != None:

            get_params['min_score'] = str(min_score)

        if disease_type != None:

            get_params['type'] = disease_type

        if disease_class != None:

            get_params['disease_class'] = disease_class

        if min_dsi != None:

            get_params['min_dsi'] = str(min_dsi)

        if max_dsi != None:

            get_params['max_dsi'] = str(max_dsi)

        if min_dpi != None:

            get_params['min_dpi'] = str(min_dpi)

        if max_dpi != None:

            get_params['max_dpi'] = str(max_dpi)

        if min_pli != None:

            get_params['min_pli'] = str(min_pli)

        if max_pli != None:

            get_params['max_pli'] = str(max_pli)

        if limit != None:

            get_params['limit'] = str(max(1, limit))

        result = self._retrieve_data(url, get_params)

        if result == None:

            return None

        GeneDiseaseAssociation = collections.namedtuple(
            'GeneDiseaseAssociation',
            [
                'geneid',
                'gene_symbol',
                'uniprotid',
                'gene_dsi',
                'gene_dpi',
                'gene_pli',
                'protein_class',
                'protein_class_name',
                'diseaseid',
                'disease_name',
                'disease_class',
                'disease_class_name',
                'disease_type',
                'disease_semantic_type',
                'score',
                'ei',
                'el',
                'year_initial',
                'year_final',
                'source',
            ],
        )

        for index, entry in enumerate(result):

            result[index] = GeneDiseaseAssociation(
                self._get_int(entry['geneid']),
                self._get_string(entry['gene_symbol']),
                self._get_string(entry['uniprotid']),
                self._get_float(entry['gene_dsi']),
                self._get_float(entry['gene_dpi']),
                self._get_float(entry['gene_pli']),
                self._get_string(entry['protein_class']),
                self._get_string(entry['protein_class_name']),
                self._get_string(entry['diseaseid']),
                self._get_string(entry['disease_name']),
                self._get_tuple(entry['disease_class'], ';'),
                self._get_tuple(entry['disease_class_name'], ';'),
                self._get_string(entry['disease_type']),
                self._get_string(entry['disease_semantic_type']),
                self._get_float(entry['score']),
                self._get_float(entry['ei']),
                self._get_string(entry['el']),
                self._get_int(entry['year_initial']),
                self._get_int(entry['year_final']),
                self._get_string(entry['source']),
            )

        return result

    def _list_to_str(self, list_obj: List[str], name: str, limit: Optional[int]=None) -> List[str]:
        '''
        Joins the list object and returns a string

        @list_obj : List[str]
            List object to be processed
        @name: str
            Name of items, like Gene ID or something
        @limit: Optional[int]
            Maximum number of list items
        '''

        if isinstance(list_obj, list):

            if limit != None and len(list_obj) > limit:

                print(f'Maximum length of {name}\'s are {limit}.')
                print(f'First {limit} {name}\'s will be used.')

                return ','.join(list_obj[:limit])

            return ','.join(list_obj)

        return list_obj

    @_if_authenticated
    @_delete_cache
    def _retrieve_data(self, url: str, get_params: Union[List[str], Dict[str, str]]) -> List[Dict[str, str]]:
        '''
        Retrieves the data with given request body

        @url : str
            Query url
        @get_params: Union[List[str], Dict[str, str]]
            GET request parameters
        '''

        headers = ['accept: */*', f'Authorization: Bearer {self._api_key}']

        get_params_extend = list()

        try:

            if isinstance(get_params['disease_class'], list):

                get_params_extend = [
                    f'disease_class={value}' for value in get_params['disease_class']
                ]

                del get_params['disease_class']

        except KeyError:

            pass

        get_params['format'] = 'json'
        get_params = [f'{key}={value}' for key, value in get_params.items()]

        if get_params_extend:

            get_params.extend(get_params_extend)

        c = curl.Curl(url=url, get=get_params, req_headers=headers)

        if c.status == 0 or c.status == 200:

            result = c.result
            result = json.loads(result)

            return result

        print(f'An error occurred with the code {c.status}')

    def _get_int(self, str_obj) -> int:
        '''
        Returns an int if the object is not None

        @str_obj : str
            String to be processed
        '''

        if str_obj != None and not isinstance(str_obj, int):

            return int(str_obj)

        return str_obj

    def _get_float(self, str_obj) -> float:
        '''
        Returns a float if the object is not None

        @str_obj : str
            String to be processed
        '''

        if str_obj != None and not isinstance(str_obj, float):

            return float(str_obj)

        return str_obj

    def _get_string(self, obj) -> str:
        '''
        Returns a string if the object is not None

        @obj : object
            Object to be processed
        '''

        if obj != None and not isinstance(obj, str):

            return str(obj).strip()

        return obj

    def _get_tuple(self, str_obj: str, delim: str) -> Tuple[str]:
        '''
        Returns a splitted tuple with given delimiter

        @str_obj : str
            String object to be processed
        @delim : str
            Char to split the str_obj
        '''

        if str_obj != None and not isinstance(str_obj, tuple):

            return tuple([item.strip() for item in str_obj.split(delim)])

        return str_obj


@DisgenetApi._delete_cache
def variant_gene_mappings() -> Dict[str, NamedTuple(
    'VariantGeneMapping',
    [
        ('geneId', str),
        ('geneSymbol', str),
        ('sourceIds', Tuple[str]),
    ],
    )
    ]:
    '''
    Downloads and processes variant-gene mappings.
    Returns a dict where the \'snpId\' is the key.
    '''

    url = urls.urls['disgenet']['variant_gene_mappings']
    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding='utf-8',
        default_mode='r',
    )
    reader = csv.DictReader(c.result, delimiter='\t')
    mapping = dict()

    for rec in reader:

        snpId = rec.pop('snpId')

        try:

            match = False

            for index, entry in enumerate(mapping[snpId]):

                if (
                    rec['geneId'] == entry['geneId']
                    and rec['geneSymbol'] == entry['geneSymbol']
                ):

                    match = True

                    if isinstance(mapping[snpId][index]['sourceId'], list):

                        mapping[snpId][index]['sourceId'].append(rec['sourceId'])

                    else:

                        mapping[snpId][index]['sourceId'] = [
                            mapping[snpId][index]['sourceId'],
                            rec['sourceId'],
                        ]

                    break

            if not match:

                mapping[snpId].append(rec)

        except KeyError:

            mapping[snpId] = [rec]

    VariantGeneMapping = collections.namedtuple(
        'VariantGeneMapping',
        [
            'geneId',
            'geneSymbol',
            'sourceIds',
        ],
    )

    for key, values in mapping.items():

        for index, value in enumerate(values):

            mapping[key][index] = VariantGeneMapping(
                value['geneId'],
                value['geneSymbol'],
                tuple(value['sourceId']),
            )

    return mapping


@DisgenetApi._delete_cache
def disease_id_mappings() -> dict[str, NamedTuple(
        'DiseaseIdMapping',
        [
            ('name', str),
            ('vocabularies', Tuple[NamedTuple(
                            'Vocabulary',
                            [
                                ('vocabulary', str),
                                ('code', str),
                                ('vocabularyName', str),
                            ],
                        )
                    ]
                ),
            ],
        ),
    ]:
    '''
    Downloads and processes disease-id mappings.
    Returns a dict where the \'diseaseId\' is the key.
    '''

    url = urls.urls['disgenet']['disease_id_mappings']
    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding='utf-8',
        default_mode='r',
    )
    reader = csv.DictReader(c.result, delimiter='\t')
    mapping = dict()

    Vocabulary = collections.namedtuple(
        'Vocabulary',
        [
            'vocabulary',
            'code',
            'vocabularyName',
        ],
    )

    for rec in reader:

        diseaseId = rec.pop('diseaseId')
        name = rec.pop('name')
        rec = Vocabulary(
            rec['vocabulary'],
            rec['code'],
            rec['vocabularyName'],
        )
        try:

            mapping[diseaseId]['vocabularies'].append(rec)

        except KeyError:
            mapping[diseaseId] = dict()
            mapping[diseaseId]['name'] = name
            mapping[diseaseId]['vocabularies'] = [rec]

    DiseaseIdMapping = collections.namedtuple(
        'DiseaseIdMapping',
        [
            'name',
            'vocabularies',
        ],
    )

    for key, value in mapping.items():

        mapping[key] = DiseaseIdMapping(
            value['name'],
            value['vocabularies'],
        )

    return mapping


@DisgenetApi._delete_cache
def disgenet_annotations(dataset='curated'):
    '''
    Downloads and processes the list of all human disease related proteins
    from DisGeNet.
    Returns dict of dicts.

    @dataset : str
        Name of DisGeNet dataset to be obtained:
        `curated`, `literature`, `befree` or `all`.
    '''

    DisGeNetAnnotation = collections.namedtuple(
        'DisGeNetAnnotation',
        [
            'disease',
            'type',
            'score',
            'dsi',
            'dpi',
            'nof_pmids',
            'nof_snps',
            'source',
        ],
    )

    url = urls.urls['disgenet']['annotations'] % dataset
    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding='utf-8',
        default_mode='r',
    )
    reader = csv.DictReader(c.result, delimiter='\t')
    data = collections.defaultdict(set)

    for rec in reader:

        uniprots = mapping.map_name(
            rec['geneSymbol'],
            'genesymbol',
            'uniprot',
        )

        if not uniprots:

            continue

        for uniprot in uniprots:

            data[uniprot].add(
                DisGeNetAnnotation(
                    disease=rec['diseaseName'],
                    type=rec['diseaseType'],
                    score=float(rec['score']),
                    dsi=float(rec['DSI']) if rec['DSI'] else None,
                    dpi=float(rec['DPI']) if rec['DPI'] else None,
                    nof_pmids=int(rec['NofPmids']),
                    nof_snps=int(rec['NofSnps']),
                    source=tuple(x.strip() for x in rec['source'].split(';')),
                )
            )

    return dict(data)
