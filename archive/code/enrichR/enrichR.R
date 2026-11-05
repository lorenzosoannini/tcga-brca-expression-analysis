# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html

library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)

rm(list=ls())

setwd("/Users/lorenzo/Desktop/Bioinf/")
################################################
source("enrichR/getEnrichmentPlot.R")
################################################
dataset <- "BRCA"

dirRes <- "results/"

if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory",dirRes,"already exists"))
}

dirDataset <- paste0(dirRes,dataset,"/")

if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory",dirDataset,"already exists"))
}

dirEnrich <- paste0(dirDataset,"Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}
################################################
top_term <- 10
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirDataset,"DEG.txt")

dbs <- listEnrichrDbs() #lista di tutti i db

dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "DisGeNET")

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")

# keep the rows with desired direction only
input_list <- input_list[input_list$direction == "DOWN",]

input_list <- input_list$genes

df <- enrichr(input_list, dbs)

BP <- df$GO_Biological_Process_2021
MF <- df$GO_Molecular_Function_2021

KEGG <- df$KEGG_2021_Human

DISGENET <- df$DisGeNET
  
getEnrichmentPlot(BP,"GO_BP",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(MF,"GO_MF",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(KEGG,"KEGG",top_term,thr_pval,dirEnrich)
getEnrichmentPlot(DISGENET,"DISGENET",top_term,thr_pval,dirEnrich)


