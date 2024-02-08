# Script to obtain Enrichment of set genes
# Given a known targets it recovers the entrezID
#Author: Marco Espina
# Marco Espina NH project 
#Date: 21 sept 2021

# Libraries----
#BiocManager::install("biomaRt")
library(tidyverse)
library(readxl)
library(biomaRt)
#library(AnnotationDbi) # if biomaRt is not working
#library(org.Hs.eg.db) # if biomaRt is not working
library(multiMiR) # la versión más actual debe ser la 1.12.0, pero no influye en las bases de datos

# Read excel file ----
setdw("~/MarcoEspina/EnrichmentwithR")
my_data <- read_excel("Diferencial de genes comunes para scriptR.xlsx") 
my_data

# Obtaining entrez ID of the targets ----
#   Using bioMart ----
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#listFilters(ensembl)
listAttributes(ensembl)

#Genes DEG comunes NH vs FPI ----
Genes_DEG_comunes_NHvsFPI <- my_data$`Genes DEG comunes NHvsFPI`
entrezIDs_Genes_DEG_comunes_NHvsFPI <- getBM(attributes = c('entrezgene_id', 'hgnc_id', 'hgnc_symbol', 'ensembl_gene_id'), 
                             filters = 'hgnc_symbol', 
                             values = Genes_DEG_comunes_NHvsFPI, 
                             mart = ensembl)
entrezIDs_Genes_DEG_comunes_NHvsFPI


#Genes DEG comunes NH vs Ctrls ----
Genes_DEG_comunes_NHvsCtrls <- my_data$`Genes DEG comunes NH vs Ctrl`
entrezIDs_Genes_DEG_comunes_NHvsCtrls <- getBM(attributes = c('entrezgene_id', 'hgnc_id', 'hgnc_symbol', 'ensembl_gene_id'), 
                                             filters = 'hgnc_symbol', 
                                             values = Genes_DEG_comunes_NHvsCtrls, 
                                             mart = ensembl)
entrezIDs_Genes_DEG_comunes_NHvsCtrls

#Genes DEG comunes FPI vs Ctrls ----
Genes_DEG_comunes_FPIvsCtrls <- my_data$`Genes DEG comunes FPIvsCtrl`
entrezIDs_Genes_DEG_comunes_FPIvsCtrls <- getBM(attributes = c('entrezgene_id', 'hgnc_id', 'hgnc_symbol', 'ensembl_gene_id'), 
                                             filters = 'hgnc_symbol', 
                                             values = Genes_DEG_comunes_FPIvsCtrls, 
                                             mart = ensembl)
entrezIDs_Genes_DEG_comunes_FPIvsCtrls

# Enrichment with enrichR ----
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021", "Reactome_2016", 
         "WikiPathways_2019_Human", "MSigDB_Hallmark_2020")

# comunes NH vs FPI ---- 

enr_Genes_DEG_comunes_NHvsFPI <- enrichr(Genes_DEG_comunes_NHvsFPI, dbs)
write.csv(enr_Genes_DEG_comunes_NHvsFPI$KEGG_2021_Human, "DEG comunes/enr_comunesNHvsFPI_KEGG.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsFPI$GO_Biological_Process_2021, "DEG comunes/enr_comunesNHvsFPI_GOBP.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsFPI$Reactome_2016, "DEG comunes/enr_comunesNHvsFPI_Reactome.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsFPI$WikiPathways_2019_Human, "DEG comunes/enr_comunesNHvsFPI_Wiki.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsFPI$MSigDB_Hallmark_2020, "DEG comunes/enr_comunesNHvsFPI_Hallmark.csv", row.names=T)

# comunes NH vs Ctrls ---- 

enr_Genes_DEG_comunes_NHvsCtrls <- enrichr(Genes_DEG_comunes_NHvsCtrls, dbs)
write.csv(enr_Genes_DEG_comunes_NHvsCtrls$KEGG_2021_Human, "DEG comunes/enr_comunesNHvsCtrl_KEGG.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsCtrls$GO_Biological_Process_2021, "DEG comunes/enr_comunesNHvsCtrl_GOBP.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsCtrls$Reactome_2016, "DEG comunes/enr_comunesNHvsCtrl_Reactome.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsCtrls$WikiPathways_2019_Human, "DEG comunes/enr_comunesNHvsCtrl_Wiki.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_NHvsCtrls$MSigDB_Hallmark_2020, "DEG comunes/enr_comunesNHvsCtrl_Hallmark.csv", row.names=T)

# comunes FPI vs Ctrls ---- 

enr_Genes_DEG_comunes_FPIvsCtrls <- enrichr(Genes_DEG_comunes_FPIvsCtrls, dbs)
write.csv(enr_Genes_DEG_comunes_FPIvsCtrls$KEGG_2021_Human, "DEG comunes/enr_comunesFPIvsCtrl_KEGG.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_FPIvsCtrls$GO_Biological_Process_2021, "DEG comunes/enr_comunesFPIvsCtrl_GOBP.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_FPIvsCtrls$Reactome_2016, "DEG comunes/enr_comunesFPIvsCtrl_Reactome.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_FPIvsCtrls$WikiPathways_2019_Human, "DEG comunes/enr_comunesFPIvsCtrl_Wiki.csv", row.names=T)
write.csv(enr_Genes_DEG_comunes_FPIvsCtrls$MSigDB_Hallmark_2020, "DEG comunes/enr_comunesFPIvsCtrl_Hallmark.csv", row.names=T)