from pypath.share import curl
from pypath.resources.urls import urls


def ctdbase_general(url):

    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding="utf-8",
        default_mode="r",
        compressed=True,
        compr="gz",
    )

    relations = list()
    fieldnames = None

    for line in c.result:

        if line.startswith("#"):

            line = line.strip(" #\n")
            line = line.split("\t")

            if len(line) > 1:
                fieldnames = line
                fieldnames = [fieldname for fieldname in fieldnames if fieldname != '']

            continue

        data = line.split("\t")

        if data[-1] == "\n":
            del data[-1]

        for i, v in enumerate(data):

            is_list = "|" in v
            has_sublist = "^" in v

            if is_list:
                v = v.split("|")
            
                if has_sublist:
                    v = [element.split("^") for element in v]

            elif has_sublist:
                v = [v.split("^")]

            data[i] = v

        data = {
            fieldname: element if element != "" else None
            for (fieldname, element) in zip(fieldnames, data)
        }

        relations.append(data)

    return relations


def chemical_gene():
    """
    Fields:

    ChemicalName
    ChemicalID (MeSH identifier)
    CasRN (CAS Registry Number, if available)
    GeneSymbol
    GeneID (NCBI Gene identifier)
    GeneForms (list)
    Organism (scientific name)
    OrganismID (NCBI Taxonomy identifier)
    Interaction
    InteractionActions (list)
    PubMedIDs (list)
    """

    url = urls['ctdbase']['url'] % "CTD_chem_gene_ixns.tsv.gz"
    return ctdbase_general(url)


def chemical_disease():
    """
    Fields:

    ChemicalName
    ChemicalID (MeSH identifier)
    CasRN (CAS Registry Number, if available)
    DiseaseName
    DiseaseID (MeSH or OMIM identifier)
    DirectEvidence (list)
    InferenceGeneSymbol
    InferenceScore
    OmimIDs (list)
    PubMedIDs (list)
    """

    url = urls['ctdbase']['url'] % "CTD_chemicals_diseases.tsv.gz"
    return ctdbase_general(url)


def disease_pathway():
    """
    Fields:

    DiseaseName
    DiseaseID (MeSH or OMIM identifier)
    PathwayName
    PathwayID (KEGG or REACTOME identifier)
    InferenceGeneSymbol (a gene via which the association is inferred)
    """

    url = urls['ctdbase']['url'] % "CTD_diseases_pathways.tsv.gz"
    return ctdbase_general(url)


def chemical_phenotype():
    """
    Fields:

    ChemicalName
    ChemicalID (MeSH identifier)
    CASRN (CAS Registry Number, if available)
    PhenotypeName
    PhenotypeID (GO identifier)
    CoMentionedTerms (list ) entries formatted as {Name, Id, Source}
    Organism (scientific name)
    OrganismID (NCBI Taxonomy identifier)
    Interaction
    InteractionActions (list) {Interaction, Action}
    AnatomyTerms (MeSH term; list) entries formatted as {SequenceOrder, Name, Id}
    InferenceGeneSymbols (list) entries formatted as {Name, Id}
    PubMedIDs (list)
    """

    url = urls['ctdbase']['url'] % "CTD_pheno_term_ixns.tsv.gz"
    result = ctdbase_general(url)

    result = modify_dict(result,
        ('comentionedterms',        ['name', 'id', 'source']),
        ('anatomyterms',            ['sequenceorder', 'name', 'id']),
        ('inferencegenesymbols',    ['name', 'id']),
        ('interactionactions',       ['interaction', 'action']),
    )

    return result


def gene_disease():
    """
    Fields:

    GeneSymbol
    GeneID (NCBI Gene identifier)
    DiseaseName
    DiseaseID (MeSH or OMIM identifier)
    DirectEvidence (list)
    InferenceChemicalName
    InferenceScore
    OmimIDs (list)
    PubMedIDs (list)
    """

    url = urls['ctdbase']['url'] % "CTD_genes_diseases.tsv.gz"
    return ctdbase_general(url)


def chemical_vocabulary():
    """
    Fields:

    ChemicalName
    ChemicalID (MeSH identifier)
    CasRN (CAS Registry Number, if available)
    Definition
    ParentIDs (identifiers of the parent terms; list)
    TreeNumbers (identifiers of the chemical's nodes; list)
    ParentTreeNumbers (identifiers of the parent nodes; list)
    Synonyms (list)
    """

    url = urls['ctdbase']['url'] % "CTD_chemicals.tsv.gz"
    return ctdbase_general(url)


def gene_vocabulary():
    """
    Fields:

    GeneSymbol
    GeneName
    GeneID (NCBI Gene identifier)
    AltGeneIDs (alternative NCBI Gene identifiers; list)
    Synonyms (list)
    BioGRIDIDs (list)
    PharmGKBIDs (list)
    UniprotIDs (list)
    """

    url = urls['ctdbase']['url'] % "CTD_genes.tsv.gz"
    return ctdbase_general(url)


def disease_vocabulary():
    """
    Fields (non-OBO):

    DiseaseName
    DiseaseID (MeSH or OMIM identifier)
    Definition
    AltDiseaseIDs (alternative identifiers; list)
    ParentIDs (identifiers of the parent terms; list)
    TreeNumbers (identifiers of the disease's nodes; list)
    ParentTreeNumbers (identifiers of the parent nodes; list)
    Synonyms (list)
    SlimMappings (MEDIC-Slim mappings; list)
    """

    url = urls['ctdbase']['url'] % "CTD_diseases.tsv.gz"
    return ctdbase_general(url)


def pathway_vocabulary():
    """
    Fields:

    PathwayName
    PathwayID (KEGG or REACTOME identifier)
    """

    url = urls['ctdbase']['url'] % "CTD_pathways.tsv.gz"
    return ctdbase_general(url)


def anatomy_vocabulary():
    """
    Fields:

    AnatomyName
    AnatomyID (MeSH identifier)
    Definition
    AltAnatomyIDs (alternative identifiers; list)
    ParentIDs (identifiers of the parent terms; list)
    TreeNumbers (identifiers of the anatomical term's nodes; list)
    ParentTreeNumbers (identifiers of the parent nodes; list)
    Synonyms (list)
    ExternalSynonyms (list)
    """

    url = urls['ctdbase']['url'] % "CTD_anatomy.tsv.gz"
    return ctdbase_general(url)

def modify_dict(dict_list, *entry_pairs):

    for i, element in enumerate(dict_list):

        for key, new_keys in entry_pairs:

            element[key] = map_keys(
                new_keys,
                element[key]
            )

        dict_list[i] = element
    
    return dict_list

def map_keys(keys, entry):

    if entry == None:
        return None
    
    result = list()

    for values in entry:

        temp_dict = dict()

        for key, value in zip(keys, values):
            temp_dict[key] = value
        
        result.append(temp_dict)

    return result
