class DisgenetOld(_auth.DisgenetAuth):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def ddas_that_share_genes(
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
                ('pvalue_jaccard_genes', float),
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
        """
        Query disease-disease associations.

        Args:
            disease:
                If vocabulary is given:
                    Disease id (ICD9CM, ICD10,MeSH, OMIM, DO, EFO,
                    NCI, HPO, MONDO, or ORDO identifier) or list of disease ids up to 100.
                Otherwise:
                    Disease id (UMLS CUI) or list of disease ids up to 100.
            vocabulary:
                Disease Vocabulary.
                Available values : icd9cm, icd10, mesh, omim, do, efo,
                nci, hpo, mondo, ordo
            source:
                Source of the DDA.
                Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE,
                CGI, CLINGEN, CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND,
                GWASCAT, GWASDB, HPO, LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
            p_value:
                P value associated to the Jaccard Index based on the shared genes.
            limit:
                Number of associated diseases to retrieve, 10 by default.
        """

        return self._get_ddas(
            disease = disease,
            share = 'genes',
            vocabulary = vocabulary,
            source = source,
            p_value = p_value,
            limit = limit,
        )


    def ddas_that_share_variants(
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
                ('pvalue_jaccard_variants', float),
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
        """
        Query disease-disease associations.

        Args:
            disease:
                if vocabulary is given:
                    Disease id (ICD9CM, ICD10,MeSH, OMIM, DO, EFO,
                    NCI, HPO, MONDO, or ORDO identifier) or list of disease ids up to 100.
                else:
                    Disease id (UMLS CUI) or list of disease ids up to 100.
            vocabulary:
                Disease Vocabulary.
                Available values : icd9cm, icd10, mesh, omim, do, efo,
                nci, hpo, mondo, ordo
            source:
                Source of the DDA.
                Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE,
                CGI, CLINGEN, CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND,
                GWASCAT, GWASDB, HPO, LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
            p_value:
                p value associated to the Jaccard Index based on the shared genes.
            limit:
                Number of associated diseases to retrieve. 10 by default.
        """

        return self._get_ddas(
            disease = disease,
            share = 'variants',
            vocabulary = vocabulary,
            source = source,
            p_value = p_value,
            limit = limit,
        )


    def vdas_by_variants(
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
        """
        Variant-disease associations by variant.

        Args:
            variant:
                Variant (dbSNP Identifier) or list of variants.
            disease:
                if vocabulary is given:
                    Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                    or ORDO identifier) or list of disease ids to filter the results.
                else:
                    Disease id (UMLS CUI) or list of diseases to filter the results.
            gene:
                Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
            source:
                Source of the VDA.
                Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
            min_score:
                Min value of the variant-disease score range.
            max_score:
                Max value of the variant-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the variant.
            max_dsi:
                Max value of the DSI range for the variant.
            min_dpi:
                Min value of the DPI range for the variant.
            max_dpi:
                Max value of the DPI range for the variant.
            limit:
                Number of VDAs to retrieve.
        """

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID')
        variant = self._list_to_str(variant, 'Variant ID', limit = 100)

        return self._get_vdas(
            gene = gene,
            disease = disease,
            variant = variant,
            vocabulary = None,
            by = 'variant',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            limit = limit,
        )


    def vdas_by_genes(
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
        """
        Variant-disease associations by gene.

        Args:
            gene:
                Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes.
            variant:
                Variant (dbSNP Identifier) or list of variants to filter the results.
            disease:
                if vocabulary is given:
                    Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                    or ORDO identifier) or list of disease ids to filter the results.
                else:
                    Disease id (UMLS CUI) or list of diseases to filter the results.
            source:
                Source of the VDA.
                Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
            min_score:
                Min value of the variant-disease score range.
            max_score:
                Max value of the variant-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the variant.
            max_dsi:
                Max value of the DSI range for the variant.
            min_dpi:
                Min value of the DPI range for the variant.
            max_dpi:
                Max value of the DPI range for the variant.
            limit:
                Number of VDAs to retrieve.
        """

        gene = self._list_to_str(gene, 'Gene ID', limit = 100)
        disease = self._list_to_str(disease, 'Disease ID')
        variant = self._list_to_str(variant, 'Variant ID')

        return self._get_vdas(
            gene = gene,
            disease = disease,
            variant = variant,
            vocabulary = None,
            by = 'gene',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            limit = limit,
        )


    def vdas_by_diseases(
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
        """
        Variant-disease associations by disease.

        Args:
            disease:
                if vocabulary is given:
                    Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                    or ORDO identifier) or list of disease ids up to 100.
                else:
                    Disease id (UMLS CUI) or list of diseases up to 100.
            vocabulary:
                Disease Vocabulary.
                Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
            variant:
                Variant (dbSNP Identifier) or list of variants to filter the results.
            gene:
                Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
            source:
                Source of the VDA.
                Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
            min_score:
                Min value of the variant-disease score range.
            max_score:
                Max value of the variant-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the variant.
            max_dsi:
                Max value of the DSI range for the variant.
            min_dpi:
                Min value of the DPI range for the variant.
            max_dpi:
                Max value of the DPI range for the variant.
            limit:
                Number of VDAs to retrieve.
        """

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID', limit = 100)
        variant = self._list_to_str(variant, 'Variant ID')

        return self._get_vdas(
            gene = gene,
            disease = disease,
            variant = variant,
            vocabulary = vocabulary,
            by = 'disease',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            limit = limit,
        )

    def vdas_by_source(
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
        """
        Variant-disease associations by source.

        Args:
            source:
                Source of the VDA.
                Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
            disease:
                Disease id (UMLS CUI) or list of diseases to filter the results..
            variant:
                Variant (dbSNP Identifier) or list of variants to filter the results..
            gene:
                Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results..
            min_score:
                Min value of the variant-disease score range.
            max_score:
                Max value of the variant-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the variant.
            max_dsi:
                Max value of the DSI range for the variant.
            min_dpi:
                Min value of the DPI range for the variant.
            max_dpi:
                Max value of the DPI range for the variant.
            limit:
                Number of VDAs to retrieve.
        """

        disease = self._list_to_str(disease, 'Disease ID')
        variant = self._list_to_str(variant, 'Variant ID')
        gene = self._list_to_str(gene, 'Gene ID')

        return self._get_vdas(
            gene = gene,
            disease = disease,
            variant = variant,
            vocabulary = vocabulary,
            by = 'source',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            limit = limit,
        )


    def gdas_by_genes(
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
        """
        Gene-disease associations by gene.

        Args:
            gene:
                Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes up to 100.
            disease:
                Disease id (UMLS CUI) or list of disease ids up to 100.
            source:
                Source of the GDA.
                Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
                CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
                LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
            min_score:
                Min value of the gene-disease score range.
            max_score:
                Max value of the gene-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the gene.
            max_dsi:
                Max value of the DSI range for the gene.
            min_dpi:
                Min value of the DPI range for the gene.
            max_dpi:
                Max value of the DPI range for the gene.
            min_pli:
                Min value of the pLI range.
            max_pli:
                Max value of the pLI range.
            limit:
                Number of GDAs to retrieve.
        """

        gene = self._list_to_str(gene, 'Gene ID', limit = 100)
        disease = self._list_to_str(disease, 'Disease ID', limit = 100)

        return self._get_gdas(
            gene = gene,
            disease = disease,
            uniprot = None,
            vocabulary = None,
            by = 'gene',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            min_pli = min_pli,
            max_pli = max_pli,
            limit = limit,
        )

    def gdas_by_diseases(
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
        """
        Gene-disease associations by disease.

        Args:
            disease:
                if vocabulary is given:
                    Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                    or ORDO identifier) or list of disease ids up to 100.
                else:
                    Disease id (UMLS CUI) or list of diseases up to 100.
            vocabulary:
                Disease Vocabulary.
                Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
            gene:
                Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
            by:
                Return associations by:
                Avaliable values : genes, disease, uniprot, source
            source:
                Source of the GDA.
                Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
                CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
                LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
            min_score:
                Min value of the gene-disease score range.
            max_score:
                Max value of the gene-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the gene.
            max_dsi:
                Max value of the DSI range for the gene.
            min_dpi:
                Min value of the DPI range for the gene.
            max_dpi:
                Max value of the DPI range for the gene.
            min_pli:
                Min value of the pLI range.
            max_pli:
                Max value of the pLI range.
            limit:
                Number of GDAs to retrieve.
        """

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID', limit = 100)

        return self._get_gdas(
            gene = gene,
            disease = disease,
            uniprot = None,
            vocabulary = vocabulary,
            by = 'disease',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            min_pli = min_pli,
            max_pli = max_pli,
            limit = limit,
        )


    def gdas_by_uniprots(
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
        """
        Gene-disease associations by UniProt ID.

        Args:
            uniprot:
                UniProt identifier or list of UniProt identifiers up to 100.
            disease:
                Disease id (UMLS CUI) or list of diseases to filter the results.
            source:
                Source of the GDA.
                Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
                CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
                LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
            min_score:
                Min value of the gene-disease score range.
            max_score:
                Max value of the gene-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the gene.
            max_dsi:
                Max value of the DSI range for the gene.
            min_dpi:
                Min value of the DPI range for the gene.
            max_dpi:
                Max value of the DPI range for the gene.
            min_pli:
                Min value of the pLI range.
            max_pli:
                Max value of the pLI range.
            limit:
                Number of GDAs to retrieve.
        """

        uniprot = self._list_to_str(uniprot, 'Uniprot ID', limit = 100)
        disease = self._list_to_str(disease, 'Disease ID')

        return self._get_gdas(
            gene = None,
            disease = disease,
            uniprot = uniprot,
            vocabulary = None,
            by = 'uniprot',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            min_pli = min_pli,
            max_pli = max_pli,
            limit = limit,
        )


    def gdas_by_source(
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
        """
        Gene-disease associations by source.

        Args:
            source:
                Source of the GDA.
                Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
                CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
                LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
            gene:
                Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes to filter the results.
            disease:
                Disease id (UMLS CUI) or list of diseases to filter the results.
            min_score:
                Min value of the gene-disease score range.
            max_score:
                Max value of the gene-disease score range.
            min_ei:
                Min value of the evidence index range.
            max_ei:
                Max value of the evidence index range.
            disease_type:
                DisGeNET Disease Type.
                Available values : disease, phenotype, group
            disease_class:
                MeSH Disease Classes
                Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
                C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
                F01, F02, F03
            min_dsi:
                Min value of the DSI range for the gene.
            max_dsi:
                Max value of the DSI range for the gene.
            min_dpi:
                Min value of the DPI range for the gene.
            max_dpi:
                Max value of the DPI range for the gene.
            min_pli:
                Min value of the pLI range.
            max_pli:
                Max value of the pLI range.
            limit:
                Number of GDAs to retrieve.
        """

        gene = self._list_to_str(gene, 'Gene ID')
        disease = self._list_to_str(disease, 'Disease ID')

        return self._get_gdas(
            gene = gene,
            disease = disease,
            uniprot = None,
            vocabulary = None,
            by = 'source',
            source = source,
            min_score = min_score,
            max_score = max_score,
            min_ei = min_ei,
            max_ei = max_ei,
            disease_type = disease_type,
            disease_class = disease_class,
            min_dsi = min_dsi,
            max_dsi = max_dsi,
            min_dpi = min_dpi,
            max_dpi = max_dpi,
            min_pli = min_pli,
            max_pli = max_pli,
            limit = limit,
        )


    def vda_evidences_by_variant(
            self,
            variant,
            gene = None,
            disease = None,
            source: str = None,
            min_year: int = None,
            max_year: int = None,
            min_score: float = None,
            max_score: float = None,
            limit: int = None,
            offset: int = None,
            get_all: bool = True
        ):

        variant = self._list_to_str(variant, 'variant ID', limit = 100)

        if disease != None:

            disease = self._list_to_str(disease, 'disease ID', limit = 100)

        if gene != None:

            gene = self._list_to_str(gene, 'gene ID')

        return self._get_evidences(
            of = 'vda',
            by = 'variant',
            gene = gene,
            disease = disease,
            variant = variant,
            source = source,
            min_year = min_year,
            max_year = max_year,
            min_score = min_score,
            max_score = max_score,
            limit = limit,
            offset = offset,
            get_all = get_all
        )


    def vda_evidences_by_disease(
            self,
            disease,
            variant = None,
            gene = None,
            source: str = None,
            min_year: int = None,
            max_year: int = None,
            min_score: float = None,
            max_score: float = None,
            limit: int = None,
            offset: int = None,
            get_all: bool = True
        ):

        disease = self._list_to_str(disease, 'disease ID', limit = 100)

        if variant != None:

            variant = self._list_to_str(variant, 'variant ID', limit = 100)

        if gene != None:

            gene = self._list_to_str(gene, 'gene ID')

        return self._get_evidences(
            of = 'vda',
            by = 'disease',
            gene = gene,
            disease = disease,
            variant = variant,
            source = source,
            min_year = min_year,
            max_year = max_year,
            min_score = min_score,
            max_score = max_score,
            limit = limit,
            offset = offset,
            get_all = get_all
        )


    def gda_evidences_by_gene(
            self,
            gene,
            disease = None,
            source: str = None,
            min_year: int = None,
            max_year: int = None,
            min_score: float = None,
            max_score: float = None,
            limit: int = None,
            offset: int = None,
            get_all: bool = True
        ):

        gene = self._list_to_str(gene, 'gene ID', limit = 100)

        if disease != None:

            disease = self._list_to_str(disease, 'disease ID', limit = 100)

        return self._get_evidences(
            of = 'gda',
            by = 'gene',
            gene = gene,
            disease = disease,
            variant = None,
            source = source,
            min_year = min_year,
            max_year = max_year,
            min_score = min_score,
            max_score = max_score,
            limit = limit,
            offset = offset,
            get_all = get_all
        )


    def gda_evidences_by_disease(
            self,
            disease,
            gene = None,
            source: str = None,
            min_year: int = None,
            max_year: int = None,
            min_score: float = None,
            max_score: float = None,
            limit: int = None,
            offset: int = None,
            get_all: bool = True
        ):

        disease = self._list_to_str(disease, 'disease ID', limit = 100)

        if gene != None:

            gene = self._list_to_str(gene, 'gene ID', limit = 100)

        return self._get_evidences(
            of = 'vda',
            by = 'disease',
            gene = gene,
            disease = disease,
            variant = None,
            source = source,
            min_year = min_year,
            max_year = max_year,
            min_score = min_score,
            max_score = max_score,
            limit = limit,
            offset = offset,
            get_all = get_all
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
                ('disease1_nshare', int),
                ('disease2_nshare', int),
                ('disease1_disease_class', Tuple[str]),
                ('disease2_disease_class', Tuple[str]),
                ('disease1_disease_class_name', Tuple[str]),
                ('disease2_disease_class_name', Tuple[str]),
                ('jaccard_share', float),
                ('pvalue_jaccard_share', float),
                ('source', str),
                ('ngenes1', int),
                ('ngenes2', int),
                ('nshare', int),
                ('nvariants1', int),
                ('nvariants2', int),
                ('diseaseid1', str),
                ('diseaseid2', str),
            ],
        ):
        """
        Returns Disease-Disease Associations with given query.

        disease:
            if vocabulary is given:
                Disease id (ICD9CM, ICD10,MeSH, OMIM, DO, EFO,
                NCI, HPO, MONDO, or ORDO identifier) or list of disease ids
            else:
                Disease id (UMLS CUI) or list of disease ids
        share:
            Return associations that share:
            Avaliable values : genes, variants
        vocabulary:
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo,
            nci, hpo, mondo, ordo
        source:
            Source of the DDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE,
            CGI, CLINGEN, CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND,
            GWASCAT, GWASDB, HPO, LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        p_value:
            p value associated to the Jaccard Index based on the shared genes.
        limit:
            Number of associated diseases to retrieve.
            Default value : 10
        """

        disease = self._list_to_str(disease, 'disease ID', limit = 100)

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
                self._get_float(entry[f'pvalue_jaccard_{share}']),
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
        """
        Returns Variant-Disease Associations with given query.

        gene:
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes.
        disease:
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids.
            else:
                Disease id (UMLS CUI) or list of diseases.
        variant:
            Variant (dbSNP Identifier) or list of variants.
        vocabulary:
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
        by:
            Return associations by:
            Avaliable values : gene, disease, variant, source
        source:
            Source of the VDA.
            Available values : CURATED, BEFREE, ALL, CLINVAR, GWASCAT, GWASDB, UNIPROT
        min_score:
            Min value of the variant-disease score range.
        max_score:
            Max value of the variant-disease score range.
        min_ei:
            Min value of the evidence index range.
        max_ei:
            Max value of the evidence index range.
        disease_type:
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        disease_class:
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        min_dsi:
            Min value of the DSI range for the variant.
        max_dsi:
            Max value of the DSI range for the variant.
        min_dpi:
            Min value of the DPI range for the variant.
        max_dpi:
            Max value of the DPI range for the variant.
        limit:
            Number of VDAs to retrieve.
        """

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
            self._log('Problem in function call. Check arguments.')

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
        """
        Returns Gene-Disease Associations with given query.

        gene:
            Gene (NCBI Entrez Identifier or HGNC Symbol) or list of genes.
        disease:
            if vocabulary is given:
                Disease id (ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO,
                or ORDO identifier) or list of disease ids.
            else:
                Disease id (UMLS CUI) or list of diseases.
        uniprot:
            Disease id (UMLS CUI) or list of disease ids.
        vocabulary:
            Disease Vocabulary.
            Available values : icd9cm, icd10, mesh, omim, do, efo, nci, hpo, mondo, ordo
        by:
            Return associations by:
            Avaliable values : gene, disease, uniprot, source
        source:
            Source of the GDA.
            Available values : CURATED, INFERRED, ANIMAL_MODELS, ALL, BEFREE, CGI, CLINGEN,
            CLINVAR, CTD_human, CTD_mouse, CTD_rat, GENOMICS_ENGLAND, GWASCAT, GWASDB, HPO,
            LHGDN, MGD, ORPHANET, PSYGENET, RGD, UNIPROT
        min_score:
            Min value of the gene-disease score range.
        max_score:
            Max value of the gene-disease score range.
        min_ei:
            Min value of the evidence index range.
        max_ei:
            Max value of the evidence index range.
        disease_type:
            DisGeNET Disease Type.
            Available values : disease, phenotype, group
        disease_class:
            MeSH Disease Classes
            Available values : C01, C04, C05, C06, C07, C08, C09, C10, C11, C12,
            C13, C14, C15, C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26,
            F01, F02, F03
        min_dsi:
            Min value of the DSI range for the gene.
        max_dsi:
            Max value of the DSI range for the gene.
        min_dpi:
            Min value of the DPI range for the gene.
        max_dpi:
            Max value of the DPI range for the gene.
        min_pli:
            Min value of the pLI range.
        max_pli:
            Max value of the pLI range.
        limit:
            Number of GDAs to retrieve.
        """

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
            self._log('Problem in function call. Check arguments.')

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

    def _get_evidences(
            self,
            of: ['gda' ,'vda'],
            by: ['gene', 'disease', 'variant'],
            gene: [str, [str]] = None,
            disease: [str, [str]] = None,
            variant: [str, [str]] = None,
            source: str = None,
            min_year: int = None,
            max_year: int = None,
            min_score: float = None,
            max_score: float = None,
            limit: int = None,
            offset: int = None,
            get_all: bool = True
        ):

        url = f'{self._api_url}/{of}/evidences/{by}/'
        get_params = dict()

        if of == 'gda' and by == 'gene' and gene != None:
            url += gene

            if disease != None:
                get_params['diasease'] = disease

        elif of == 'vda' and by == 'variant' and gene != None:
            url += variant

            if gene != None:
                get_params['gene'] = gene

            if disease != None:
                get_params['disease'] = disease

        elif  by == 'disease' and disease != None:
            url += disease

            if gene != None:
                get_params['gene'] = gene

            if of == 'vda' and variant != None:
                get_params['variant'] = variant

        else:
            self._log('Problem in function call. Check arguments.')
            return None

        if source != None:
            get_params['source'] = source

        if min_year != None:
            get_params['min_year'] = str(min_year)

        if max_year != None:
            get_params['max_year'] = str(max_year)

        if min_score != None:
            get_params['min_score'] = str(min_score)

        if max_score != None:
            get_params['max_score'] = str(max_score)

        if limit != None:
            get_params['limit'] = str(max(1, limit))

        if offset != None:
            get_params['offset'] = offset

        result = []

        while True:
            data = self._retrieve_data(url, get_params)
            result.extend(data['results'])
            url = data['next']

            if not get_all or url == None:
                break

        return result

    def _list_to_str(
        self, list_obj: List[str], name: str, limit: Optional[int] = None
    ) -> List[str]:
        """
        Joins the list object and returns a string

        list_obj:
            List object to be processed
        name:
            Name of items, like Gene ID or something
        limit:
            Maximum number of list items
        """

        if isinstance(list_obj, list):

            if limit != None and len(list_obj) > limit:

                self._log(f'Maximum length of {name} are {limit}.')
                self._log(f'First {limit} items will be used.')

                return ','.join(list_obj[:limit])

            return ','.join(list_obj)

        return list_obj

    @_auth.DisgenetAuth._if_authenticated
    @_auth.DisgenetAuth._delete_cache
    def _retrieve_data(
        self, url: str, get_params: Union[List[str], Dict[str, str]]
    ) -> List[Dict[str, str]]:
        """
        Retrieves the data with given request body

        url:
            Query url
        get_params:
            GET request parameters
        """

        headers = ['accept: */*', f'Authorization: Bearer {self._api_key}']

        get_params_extend = list()

        try:
            if isinstance(get_params['disease_class'], list):
                get_params_extend = [
                    f'disease_class = {value}' for value in get_params['disease_class']
                ]

                del get_params['disease_class']

        except KeyError:
            pass

        get_params['format'] = 'json'
        get_params = [f'{key} = {value}' for key, value in get_params.items()]

        if get_params_extend:
            get_params.extend(get_params_extend)

        c = curl.Curl(url = url, get = get_params, req_headers = headers)

        if c.status == 0 or c.status == 200:
            result = c.result
            result = json.loads(result)

            return result

        self._log(f'An error occurred with the code {c.status}')

    def _get_int(self, str_obj) -> int:
        """
        Returns an int if the object is not None

        str_obj:
            String to be processed
        """

        if str_obj != None and not isinstance(str_obj, int):
            return int(str_obj)

        return str_obj

    def _get_float(self, str_obj) -> float:
        """
        Returns a float if the object is not None

        str_obj:
            String to be processed
        """

        if str_obj != None and not isinstance(str_obj, float):
            return float(str_obj)

        return str_obj

    def _get_string(self, obj) -> str:
        """
        Returns a string if the object is not None

        obj:
            Object to be processed
        """

        if obj != None and not isinstance(obj, str):
            return str(obj).strip()

        return obj

    def _get_tuple(self, str_obj: str, delim: str) -> Tuple[str]:
        """
        Returns a splitted tuple with given delimiter

        str_obj:
            String object to be processed
        delim:
            Char to split the str_obj
        """

        if str_obj != None and not isinstance(str_obj, tuple):
            return tuple([item.strip() for item in str_obj.split(delim)])

        return str_obj
