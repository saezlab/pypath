SELECT
    term.name AS go_name,
    term.term_type AS go_aspect,
    term.acc AS go_id,
    species.ncbi_taxa_id AS ncbi_tax_id,
    gene_product.symbol AS genesymbol,
    dbxref.xref_key AS uniprot
FROM gene_product
INNER JOIN
    species
    ON gene_product.species_id = species.id
INNER JOIN
    gene_product_synonym
    ON gene_product.id = gene_product_synonym.gene_product_id
INNER JOIN
    dbxref
    ON gene_product.dbxref_id = dbxref.id
INNER JOIN
    association
    ON gene_product.id = association.gene_product_id
INNER JOIN
    term
    ON association.term_id = term.id
WHERE
    species.ncbi_taxa_id = %u AND
    %s
    %s
    dbxref.xref_dbname = 'UniProtKB'
GROUP BY
    xref_key,
    term.acc;
