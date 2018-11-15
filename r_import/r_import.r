#!/usr/bin/env Rscript

# Dénes Türei RWTH Aachen 2017
# turei.denes@gmail.com

require(igraph)
require(dplyr)

url <- paste0(
    'http://omnipathdb.org/interactions?',
    'fields=sources,references&genesymbols=1'
)

download_omnipath <- function(){
    
    read.table(url, sep = '\t', header = TRUE)
    
}

get_genesymbols <- function(raw){
    
    return(rbind(raw %>% select(uniprot = source,
                                label = source_genesymbol),
                 raw %>% select(uniprot = target,
                                label = target_genesymbol)) %>%
           distinct())
    
}

omnipath_data_frame <- function(raw = NULL, directed = FALSE){
    
    if (is.null(raw)) raw <- download_omnipath()
    
    if(directed){
        
        raw <- raw %>% filter(is_directed == 1)
        
    }
    
    up_gs <- get_genesymbols(raw)
    
    raw <- raw %>% select(- c(source_genesymbol, target_genesymbol))
    
    return(list(edges = raw, vertices = up_gs))
    
}

omnipath_graph <- function(raw = NULL, directed = FALSE){

    op_dfs <- omnipath_data_frame(raw = raw, directed = directed)
    
    op_g   <- graph_from_data_frame(op_dfs$edges,
                                    directed = directed,
                                    vertices = op_dfs$vertices)
    
    E(op_g)$sources    <- strsplit(E(op_g)$sources,    ';')
    E(op_g)$references <- strsplit(E(op_g)$references, ';')
    
    return(op_g)
    
}
