
#load config file
library(rjson)
json_data<-fromJSON(file="Tool.config")

gene_length_file = json_data$config$transcriptome$gene_length_file
ROSMAP_count_file = json_data$config$transcriptome$readcounts



library(tidyverse)
library(dplyr)
#custom setting parameter

counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}


############### Transform from counts to TPMs ##############
#gene_length_file = "C:/Data/gencode.v24.annotation.gtf.gene_length_unionexion.txt" # union exon length
gene_length =read.table(gene_length_file, header=F, sep = "\t", stringsAsFactors=F)
names(gene_length) = c("ENSGID", "length")
# remove gene name version
gene_length$ENSGID= do.call(rbind, strsplit(gene_length$ENSGID, '.', fixed=T))[,1]

#ROSMAP_count_file = "C:/Data/ROSMAP_all_counts_matrix.txt"
ROSMAP_count = read.table(ROSMAP_count_file, header=T, sep = "\t", stringsAsFactors=F, check.names=F)
ROSMAP_count = ROSMAP_count[-c(1:4),]
ROSMAP_count$feature = do.call(rbind, strsplit(ROSMAP_count$feature, '.', fixed=T))[,1]
comm = intersect(gene_length$ENSGID, ROSMAP_count$feature)

rownames(ROSMAP_count) = ROSMAP_count$feature; ROSMAP_count = ROSMAP_count[,-1]
rownames(gene_length) = gene_length$ENSGID; 
ROSMAP_count = ROSMAP_count[comm,]
gene_length = gene_length[comm,]
ROSMAP_tpm = counts_to_tpm(ROSMAP_count, gene_length$length)




write.csv(ROSMAP_tpm,"TPM_matrix.csv",header = T)

