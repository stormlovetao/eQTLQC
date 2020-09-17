###########################################################
### This R script transform read count file to TPM metrics
### Input: read count file (first column is gene ID, following columns 
###        are samples' IDs; separated by tab), gene length file (two columns, separated by tab)
### Output: TPM metrics file
###########################################################
#load config file
library(rjson)
json_data<-fromJSON(file="Tool.config")

gene_length_file = json_data$config$transcriptome$gene_length_file
Read_count_file = json_data$config$transcriptome$readcounts



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

gene_length =read.table(gene_length_file, header=F, sep = "\t", stringsAsFactors=F)
names(gene_length)[1:2] = c("ENSGID", "length")
# remove gene version
gene_length$ENSGID= do.call(rbind, strsplit(gene_length$ENSGID, '.', fixed=T))[,1]

#Read_count_file load
Read_count = read.table(Read_count_file, header=T, sep = "\t", stringsAsFactors=F, check.names=F)
names(Read_count)[1] = "feature"
# remove gene version
Read_count$feature = do.call(rbind, strsplit(Read_count$feature, '.', fixed=T))[,1]
comm = intersect(gene_length$ENSGID, Read_count$feature)

rownames(Read_count) = Read_count$feature; Read_count = Read_count[,-1]
rownames(gene_length) = gene_length$ENSGID; 
Read_count = Read_count[comm,]
gene_length = gene_length[comm,]
tpm = counts_to_tpm(Read_count, gene_length$length)

write.csv(tpm,"TPM_matrix.csv",header = T)

