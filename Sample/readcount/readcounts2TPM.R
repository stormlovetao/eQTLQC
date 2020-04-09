
#load config file
library(rjson)
json_data<-fromJSON(file="Tool.config")

gene_length_file = json_data$config$transcriptome$gene_length_file
Count_file = json_data$config$transcriptome$readcounts

#Count_file = "./GD660.GeneQuantCount.txt"
#gene_length_file = "C:/Data/gencode.v24.annotation.gtf.gene_length_unionexion.txt"

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


Count = read.table(Count_file, header=T, sep = "\t", stringsAsFactors=F, check.names=F)
#user should pre-format the read count and delete extra rows and cols
#Count = Count[-c(1:4),]
#Count = Count[,-c(2:4)]
Count$TargetID = do.call(rbind, strsplit(Count$TargetID, '.', fixed=T))[,1]
comm = intersect(gene_length$ENSGID, Count$TargetID)

rownames(Count) = Count$TargetID; Count = Count[,-1]
rownames(gene_length) = gene_length$ENSGID; 
Count = Count[comm,]
gene_length = gene_length[comm,]
tpm_data = counts_to_tpm(Count, gene_length$length)




write.csv(tpm_data,"TPM_matrix.csv")

