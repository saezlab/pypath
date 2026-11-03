# PyPath Module Groups for Parallel Analysis

Total modules: 190 (divided into 10 groups of ~20 modules each)

## Group 1 (20 modules)
- transmir
- lncdisease
- icellnet
- ipi
- ca1
- wang
- ctdbase
- hpmr
- laudanna
- compleat
- signor
- dorothea
- uniprot
- phosphatome
- progeny
- expasy
- biogrid
- celltypist
- oma
- protmapper

## Group 2 (20 modules)
- reactome_old
- uniprot_idmapping
- pazar
- trrust
- acsn
- panglaodb
- proteinatlas
- netbiol
- lambert2018
- kegg
- ebi
- comppi
- encode
- homologene
- clinvar
- humancellmap
- pubmed
- cellinker
- cpad
- phosphoelm

## Group 3 (20 modules)
- graphviz
- tfcensus
- mcam
- proteins
- guide2pharma
- credentials
- sider
- hippie
- iptmnet
- humsavar
- kea
- kegg_api
- cellphonedb
- opm
- matrisome
- intogen
- cellcellinteractions
- matrixdb
- mirdeathdb
- cell

## Group 4 (20 modules)
- almen2009
- biogps
- uniprot_db
- pubchem
- corum
- deathdomain
- integrins
- adrecs
- msigdb
- lmpid
- opentargets
- kirouac2010
- drugcentral
- locate
- cytosig
- mppi
- trip
- collectri
- ncrdeathdb
- mirtarbase

## Group 5 (20 modules)
- mitab
- cosmic
- intact
- abs
- cancercellmap
- phobius
- zhong2015
- gpcrdb
- pharos
- hgnc
- cancerdrugsdb
- surfaceome
- embrace
- diseases
- pathophenodb
- topdb
- pdb
- compath
- italk
- phosphonetworks

## Group 6 (20 modules)
- twosides
- membranome
- mirecords
- spike
- pathwaycommons
- huri
- cellchatdb
- oreganno
- havugimana
- gutmgene
- cellcall
- embopress
- ensembl
- connectomedb
- common
- li2012
- scconnect
- hprd
- ddinter
- lincs

## Group 7 (20 modules)
- ielm
- pepcyber
- domino
- switches_elm
- lrdb
- signalink
- mimp
- cancersea
- ontology
- i3d
- instruct
- humap
- pfam
- eutils
- htri
- go
- science
- offsides
- interpro
- talklr

## Group 8 (20 modules)
- csa
- mir2disease
- macrophage
- dgidb
- threedid
- lncrnadb
- pro
- unichem
- string
- baccin2019
- stitch
- wojtowicz2020
- genecards
- mirbase
- phosphopoint
- elm
- cspa
- adhesome
- cpdb
- exocarta

## Group 9 (20 modules)
- ramilowski2015
- complexportal
- depod
- netpath
- negatome
- pisa
- kinasedotcom
- phosphosite
- pdzbase
- drugbank
- biomodels
- celltalkdb
- threedcomplex
- reaction
- imweb
- dip
- hpo
- dbptm
- biomart
- innatedb

## Group 10 (10 modules)
- lipidmaps/structures
- hmdb/metabolites
- hmdb/proteins
- hmdb/xml
- hmdb/common
- hmdb/visual
- hmdb/structures
- new_stitch/actions_test
- disgenet/_api/simple
- disgenet/_api/schema

## Instructions for Each Agent

Each agent should:
1. Read the categorization instructions from `categorization_instructions.md`
2. For each module in their assigned group:
   - Locate the module file in `/pypath/inputs/`
   - Read and analyze the source code
   - Look up the module's functions in `modules.json`
   - Categorize according to the YAML schema
3. Output results as a YAML file named `group_X_results.yaml` where X is the group number

## Expected Output Format

Each agent should produce a YAML file containing an array of module categorizations:

```yaml
modules:
  - module_name: example_module
    functions:
      - example_function1
      - example_function2
    maintenance_status:
      category: frequent
      last_known_update: 2024-01
      notes: "Actively maintained database"
    # ... rest of categorization
  
  - module_name: another_module
    # ... categorization
```