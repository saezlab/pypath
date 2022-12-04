from collections import namedtuple

import pandas as pd

import pypath.resources.urls as urls

def get_clinvar():
    '''
    Downloads Clinvar data -> variant summaries
    Returns list of namedtuple variants
    '''

    result = set() # first decleration is set to prevent recurrency. But at the end it will return as a list 
    Variant = namedtuple(
            'Variant',
            [
                'allele_id',
                'variant_type',
                'name',
                'gene_id',
                'gene_symbol',
                'clinical_significance',
                "rs",
                "phenotype_ids",
                "phenotype_list",
                "origin_simple",
                "variation_id"
            ],
            defaults = None
        )

    url = urls.urls['clinvar']['url']
    variants = pd.read_csv(url, compression='gzip', header=0, sep="\t", quotechar='"', low_memory=False)
  
    for _, sample in variants.iterrows():

        pi = sample['PhenotypeIDS'].replace("|", ";").split(";")
        pl = sample['PhenotypeList'].replace("|", ";").split(";")

        zipped = list(set(zip(pl, pi)))
        res = list(zip(*zipped))

        # phenotype list and phenotype ids are matched by their index
        phenotype_list_adj = res[0]
        phenotype_ids_adj = res[1]

        variant = Variant(
            allele_id = sample['#AlleleID'],
            variant_type = sample['Type'],
            name = sample['Name'],
            gene_id = sample['GeneID'],
            gene_symbol = sample['GeneSymbol'],
            clinical_significance = sample['ClinicalSignificance'],
            rs = sample['RS# (dbSNP)'],
            phenotype_ids = phenotype_ids_adj,
            phenotype_list = phenotype_list_adj,
            origin_simple = sample['OriginSimple'],
            variation_id = sample['VariationID']
        )
        result.add(variant)

    return list(result)