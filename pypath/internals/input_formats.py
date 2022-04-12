#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Enables ID translations and mapping
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

import copy

import pypath.share.settings as settings
import pypath.share.session as session
import pypath.share.constants as constants

_logger = session.Logger(name = 'input_formats')

__all__ = [
    'FileMapping',
    'PickleMapping',
    'NetworkInput',
    'ReadList',
    'Reference',
    'UniprotListMapping',
    'ProMapping',
    'ArrayMapping',
    'BiomartMapping',
]


AC_QUERY = {
    'genesymbol': ('genes', 'PREFERRED'),
    'genesymbol-syn': ('genes', 'ALTERNATIVE'),
    'hgnc': ('database', 'HGNC'),
    'embl': ('database', 'embl'),
    'entrez': ('database', 'geneid'),
    'refseqp': ('database', 'refseq'),
    'enst': ('database', 'ensembl'),
    'uniprot-entry': ('entry name', None),
    'protein-name': ('protein names', None),
    'ec': ('ec', None),
}

AC_MAPPING = {
    'uniprot': 'ACC',
    'uniprot_id': 'ID',
    'embl': 'EMBL',
    'embl_id': 'EMBL_ID',
    'pir': 'PIR',
    'entrez': 'P_ENTREZGENEID',
    'gi': 'P_GI',
    'refseqp': 'P_REFSEQ_AC',
    'refseqn': 'REFSEQ_NT_ID',
    'ensembl': 'ENSEMBL_ID',
    'ensp': 'ENSEMBL_PRO_ID',
    'enst': 'ENSEMBL_TRS_ID',
    'ensg': 'ENSEMBLGENOME_ID',
    'ensgp': 'ENSEMBLGENOME_PRO_ID',
    'ensgt': 'ENSEMBLGENOME_TRS_ID',
    'hgnc': 'HGNC_ID',
    'ensp_string': 'STRING_ID',
}

BIOMART_MAPPING = {
    'hgnc_symbol': 'hgnc_symbol',
    'rnacentral': 'rnacentral',
    'hgnc_trans_name': 'hgnc_trans_name',
    'wikigene_name': 'wikigene_name',
    'gene_name': 'external_gene_name',
    'genesymbol': 'external_gene_name',
    'transcript_name': 'external_transcript_name',
    'gene_description': 'description',
    'gene_synonym': 'external_synonym',
    'interpro_description': 'interpro_description',
    'interpro': 'interpro',
    'interpro_short_description': 'interpro_short_description',
    'enst_biomart': 'ensembl_transcript_id',
    'ensg_biomart': 'ensembl_gene_id',
    'ensp_biomart': 'ensembl_peptide_id',
    'ensembl_gene_id': 'ensembl_gene_id',
    'ensembl_transcript_id': 'ensembl_transcript_id',
    'ensembl_peptide_id': 'ensembl_peptide_id',
    'uniprot': 'uniprotswissprot',
    'trembl': 'uniprotsptrembl',

}

PRO_MAPPING = {
    'alzforum': 'Alzforum_mut',
    'araport': 'Araport',
    'cgnc': 'CGNC',
    'dictybase': 'dictyBase',
    'dto': 'DTO',
    'ecocyc': 'EcoCyc',
    'ecogene': 'EcoGene',
    'ensembl_pro': 'Ensembl',
    'ensembl_bacteria': 'EnsemblBacteria',
    'flybase': 'FlyBase',
    'hgnc': 'HGNC',
    'iuphar_fam': 'IUPHARfam',
    'iuphar': 'IUPHARobj',
    'mgi': 'MGI',
    'mro': 'MRO',
    'ncbi_gene': 'NCBIGene',
    'pbd': 'PDB',
    'pombase': 'PomBase',
    'interpro': 'PRO',
    'reactome': 'Reactome',
    'rgd': 'RGD',
    'sgd': 'SGD',
    'tdr': 'TDR',
    'uniprot': 'UniProtKB',
    'uniprot-var': 'UniProtKB_VAR',
    'wormbase': 'WormBase',
    'zfin': 'ZFIN',
}

ARRAY_MAPPING = {
    'affy',
    'affymetrix',
    'illumina',
    'agilent',
    'codelink',
    'phalanx',
}


class MappingInput(object):

    _resource_id_types = {}

    def __init__(
            self,
            type_,
            id_type_a,
            id_type_b,
            ncbi_tax_id = None,
            resource_id_type_a = None,
            resource_id_type_b = None,
        ):

        self.type = type_
        self.id_type_a = id_type_a
        self.id_type_b = id_type_b
        self.resource_id_type_a = resource_id_type_a
        self.resource_id_type_b = resource_id_type_b
        self.ncbi_tax_id = ncbi_tax_id or settings.get('default_organism')


    def _resource_id_type(self, side):

        return self.resource_id_type(
            getattr(self, 'id_type_%s' % side),
            override = getattr(self, 'resource_id_type_%s' % side),
        )


    @property
    def _resource_id_type_a(self):

        return self._resource_id_type(side = 'a')


    @property
    def _resource_id_type_b(self):

        return self._resource_id_type(side = 'b')


    @classmethod
    def resource_id_type(cls, id_type, override = None):

        return override or cls._resource_id_types.get(id_type, None)


class FileMapping(MappingInput):

    def __init__(
            self,
            id_type_a,
            id_type_b,
            input_,
            col_a,
            col_b,
            separator = None,
            header = 0,
            ncbi_tax_id = None,
            entity_type = 'protein',
        ):

        MappingInput.__init__(
            self,
            type_ = 'file',
            id_type_a = id_type_a,
            id_type_b = id_type_b,
            ncbi_tax_id = ncbi_tax_id,
        )

        self.input = input_
        self.col_a = col_a
        self.col_b = col_b
        self.separator = separator
        self.header = header
        self.entity_type = entity_type
        self.input_args = {'organism': self.ncbi_tax_id}


    def set_organism(self, ncbi_tax_id):

        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id

        if 'organism' in other_organism.input_args:

            other_organism.input_args['organism'] = ncbi_tax_id

        return other_organism


class UniprotMapping(MappingInput):

    _resource_id_type_b = 'id'

    def __init__(
            self,
            id_type_a,
            id_type_b = 'uniprot',
            ncbi_tax_id = 9606,
            swissprot = 'yes',
        ):
        """
        Defines an ID conversion table to retrieve from UniProt.

        id_type : str
            Type of accession numbers you would like to translate.
        target_id_type : str
            Type of accession numbers you would like to translate to.
        tax : int
            NCBI Taxonomy ID of the organism of interest.
        swissprot : str
            Look for SwissProt or Trembl.
            Passed directly to UniProt`s `reviewed` parameter. `yes` or `no`
            To fetch Trembl and SwissProt together, set value to None.
        mapping : bool
            Get the data from UniProt`s programmatic access query interface,
            (uniprot.org/uniprot) or the batch retrieval/id mapping service
            (uniprot.org/mapping). These have slightly different APIs and
            capabilities. Some IDs can be obtained from the former, some
            from the latter.
        """

        self.type = 'uniprot'

        MappingInput.__init__(
            self,
            type_ = 'uniprot',
            id_type_a = id_type_a,
            id_type_b = id_type_b,
            ncbi_tax_id = ncbi_tax_id,
        )

        self.ncbi_tax_id = int(ncbi_tax_id)
        self.typ = 'protein'
        self.swissprot = swissprot


    def set_organism(self, ncbi_tax_id):

        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        return other_organism


    @property
    def field(self):

        return AC_QUERY.get(self.id_type_a, (None,))[0]


    @property
    def subfield(self):

        return AC_QUERY.get(self.id_type_a, (None, None))[1]


    @staticmethod
    def resource_id_type(id_type, override = None):
        """
        For an ID type label used in pypath, returns the one used in the
        UniProt web service. If the label is not available in the built in
        list None is returned.

        Returns:
            (str): The ID type label used by UniProt; None if the input
                label is not known.
        """

        id_type = AC_QUERY.get(id_type, (None, None))
        id_type = '%s(%s)' % id_type if id_type[1] else id_type[0]

        return id_type


class UniprotListMapping(MappingInput):
    """
    Provides parameters for downloading mapping table from UniProt
    `Upload Lists` webservice.

    :arg str id_type_a:
        Custom name for one of the ID types.
    :arg str id_type_b:
        Custom name for the other ID type.
    :arg str uniprot_id_type_a:
        This is the symbol the UniProt webservice uses for the first
        name type. These are included in the module and set
        automatically, the argument only gives a way to override this.
    :arg str uniprot_id_type_b:
        Same as above just for the other ID type.
    :arg bool swissprot:
        DOwnload data only for SwissProt IDs.
    """

    _resource_id_types = AC_MAPPING

    def __init__(
            self,
            id_type_a,
            id_type_b,
            uniprot_id_type_a = None,
            uniprot_id_type_b = None,
            ncbi_tax_id = 9606,
            swissprot = None,
        ):

        MappingInput.__init__(
            self,
            type_ = 'uniprot_list',
            id_type_a = id_type_a,
            id_type_b = id_type_b,
            ncbi_tax_id = ncbi_tax_id,
            resource_id_type_a = uniprot_id_type_a,
            resource_id_type_b = uniprot_id_type_b,
        )

        self.swissprot = swissprot
        self.ac_mapping = AC_MAPPING

        self.uniprot_id_type_a = self._resource_id_type_a
        self.uniprot_id_type_b = self._resource_id_type_b

        self.entity_type = 'protein'


    def set_organism(self, ncbi_tax_id):

        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        return other_organism


class ProMapping(MappingInput):
    """
    Provides parameters for mapping table from the Protein Ontology
    Consortium.

    :arg str id_type_a:
        Custom name for one of the ID types.
    :arg str id_type_b:
        Custom name for the other ID type.
    :arg str pro_id_type_a:
        This is the symbol PRO uses to label the IDs.
        These are included in the module and set
        automatically, the argument only gives a way to override this.
    :arg str pro_id_type_b:
        Same as above just for the other ID type.
    """

    _resource_id_types = PRO_MAPPING

    def __init__(
            self,
            id_type_a,
            id_type_b = None,
            pro_id_type_a = None,
            pro_id_type_b = None,
            ncbi_tax_id = constants.NOT_ORGANISM_SPECIFIC,
        ):

        to_pro = id_type_a != 'pro'
        id_type = id_type_a if to_pro else id_type_b
        pro_id_type = (
            pro_id_type_a if to_pro else pro_id_type_b
        )

        MappingInput.__init__(
            self,
            type_ = 'pro',
            id_type_a = 'pro',
            id_type_b = id_type,
            ncbi_tax_id = -1,
            resource_id_type_a = pro_id_type_a,
            resource_id_type_b = pro_id_type_b,
        )
        self.to_pro = to_pro
        self.id_type = id_type

        self.pro_mapping = PRO_MAPPING

        self.pro_id_type = pro_id_type or self.pro_mapping[self.id_type_b]

        self.entity_type = 'protein'


class BiomartMapping(MappingInput):

    _resource_id_types = BIOMART_MAPPING

    def __init__(
            self,
            id_type_a,
            id_type_b = None,
            transcript = False,
            biomart_id_type_a = None,
            biomart_id_type_b = None,
            ncbi_tax_id = 9606,
        ):

        MappingInput.__init__(
            self,
            type_ = 'biomart',
            id_type_a = id_type_a,
            id_type_b = id_type_b,
            ncbi_tax_id = ncbi_tax_id,
            resource_id_type_a = biomart_id_type_a,
            resource_id_type_b = biomart_id_type_b,
        )

        self.biomart_id_type_a = self._resource_id_type_a
        self.biomart_id_type_b = self._resource_id_type_b
        self.attrs = (
            self.biomart_id_type_a,
            self.biomart_id_type_b,
        )

        self.biomart_mapping = BIOMART_MAPPING


class UnichemMapping(MappingInput):

    def __init__(
            self,
            id_type_a,
            id_type_b,
            ncbi_tax_id = constants.NOT_ORGANISM_SPECIFIC,
        ):

        MappingInput.__init__(
            self,
            type_ = 'unichem',
            id_type_a = id_type_a,
            id_type_b = id_type_b,
            ncbi_tax_id = constants.NOT_ORGANISM_SPECIFIC,
        )


class ArrayMapping(MappingInput):
    """
    Provides parameters for microarray probe mapping tables.

    :arg str id_type_a:
        Custom name for one of the ID types.
    :arg str id_type_b:
        Custom name for the other ID type.
    :arg str pro_id_type_a:
        This is the symbol PRO uses to label the IDs.
        These are included in the module and set
        automatically, the argument only gives a way to override this.
    :arg str pro_id_type_b:
        Same as above just for the other ID type.
    """

    _resource_id_types = ARRAY_MAPPING

    def __init__(
            self,
            id_type_a,
            id_type_b,
            ncbi_tax_id = 9606,
        ):

        MappingInput.__init__(
            self,
            type_ = 'array',
            id_type_a = self._get_id_type(id_type_a),
            id_type_b = self._get_id_type(id_type_b),
            ncbi_tax_id = ncbi_tax_id,
            resource_id_type_a = self._process_id_type(id_type_a),
            resource_id_type_b = self._process_id_type(id_type_b),
        )

        self.ensembl_id = (
            self.resource_id_type_a
                if self.id_type_a.startswith('ens') else
            self.resource_id_type_b
        )
        self.array_id = (
            self.resource_id_type_a
                if self.id_type_a in self._resource_id_types else
            self.resource_id_type_b
        )

        self.entity_type = 'protein'


    @staticmethod
    def _process_id_type(id_type):

        id_type = id_type.lower()
        id_type = 'affy' if id_type == 'affymetrix' else id_type
        id_type = 'ensg' if id_type == 'ensembl' else id_type

        if (
            id_type not in self._resource_id_types and
            id_type not in {'ensg', 'enst', 'ensp'}
        ):

            msg = (
                'Unknown ID type for microarray probe mapping: `%s`. '
                'Microarray ID types include `affy`, `illumina`, `agilent`, '
                '`codelink` and `phalanx`, all these can be translated to '
                'Ensembl gene, transcript or peptide IDs: `ensg`, `enst` '
                'or `ensp`. If you translate to some other ID type, do it '
                'in multiple steps.' % str(id_type)
            )
            _logger._log(msg)
            raise ValueError(msg)

        return id_type


class PickleMapping(MappingInput):


    def __init__(
            self,
            id_type_a,
            id_type_b,
            fname,
            ncbi_tax_id = None,
        ):

        MappingInput.__init__(
            self,
            type_ = 'pickle',
            id_type_a = id_type_a,
            id_type_b = id_type_b,
            ncbi_tax_id = ncbi_tax_id,
        )

        self.fname = fname


class NetworkInput:


    def __init__(
            self,
            name = "unknown",
            separator = None,
            id_col_a = 0,
            id_col_b = 1,
            id_type_a = "uniprot",
            id_type_b = "uniprot",
            entity_type_a = "protein",
            entity_type_b = "protein",
            is_directed = False,
            sign = False,
            input = None,
            references = None,
            extra_edge_attrs = None,
            extra_node_attrs_a = None,
            extra_node_attrs_b = None,
            header = False,
            taxon_a = 9606,
            taxon_b = 9606,
            ncbi_tax_id = 9606,
            interaction_type = 'post_translational',
            positive_filters = None,
            negative_filters = None,
            mark_source  =  None,
            mark_target  =  None,
            input_args = None,
            curl_args = None,
            must_have_references = True,
            huge = False,
            resource = None,
            unique_fields = None,
            expand_complexes = None,
            data_model = None,
            allow_loops = False,
            only_default_organism = False,
        ):
        """
        :param str mark_source:
            Creates a boolean vertex attribute and sets it True for the
            source vertex of directed interactions from this particular
            resource.
        :param str mark_target:
            Same as ``mark_source`` but for target vertices.
        """

        self.entity_type_a = entity_type_a
        self.entity_type_b = entity_type_b
        self.id_col_a = id_col_a
        self.id_col_b = id_col_b
        self.id_type_a = id_type_a
        self.id_type_b = id_type_b
        self.is_directed = is_directed
        self.input = input
        self.extra_edge_attrs = extra_edge_attrs or {}
        self.extra_node_attrs_a = extra_node_attrs_a or {}
        self.extra_node_attrs_b = extra_node_attrs_b or {}
        self.name = name
        self.separator = separator
        self.header = header
        self.refs = references or None
        self.sign = sign
        self.taxon_a = taxon_a
        self.taxon_b = taxon_b
        self.ncbi_tax_id = ncbi_tax_id
        self.interaction_type = interaction_type
        self.positive_filters = positive_filters or []
        self.negative_filters = negative_filters or []
        self.input_args = input_args or {}
        self.curl_args = curl_args or {}
        self.must_have_references = must_have_references and bool(references)
        self.huge = huge
        self.resource = self.name if resource is None else resource
        self.mark_source = mark_source
        self.mark_target = mark_target
        self.unique_fields = unique_fields or set()
        self.expand_complexes = expand_complexes
        self.data_model = data_model
        self.allow_loops = allow_loops
        self.only_default_organism = only_default_organism


    def _field(self, value, cls):

        return value if isinstance(value, cls) else cls(compact = value)


class ReadList:


    def __init__(
            self,
            name = 'unknown',
            separator = None,
            id_col = 0,
            id_type = 'uniprot',
            entity_type = 'protein',
            input = None,
            extra_attrs = None,
            header = False,
        ):

        self.entity_type = entity_type
        self.id_col = id_col
        self.id_type = id_type
        self.input = input
        self.extra_attrs = extra_attrs or {}
        self.name = name
        self.separator = separator
        self.header = header
