wget https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD660.GeneQuantCount.txt.gz

gunzip GD660.GeneQuantCount.txt.gz

# Format the readcount file:
# Tab seperator, has header name, the first column is gene ID, and following columns are samples
echo 'data = read.table("GD660.GeneQuantCount.txt", header=T, sep = "\t", stringsAsFactors=F, check.names=F)' > format_input.R
echo 'names(data)[1] = "feature"' >> format_input.R
echo 'data = data[,-c(2:4)] # remove not to use columns' >> format_input.R
echo 'write.table(data, file="GD660.GeneQuantCount.txt", row.names=F, col.names=T, sep = "\t", quote=F)' >> format_input.R
Rscript format_input.R
rm format_input.R

#python3 Tool.py







