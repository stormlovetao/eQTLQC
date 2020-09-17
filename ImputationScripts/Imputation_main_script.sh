
# Steps for preparation of imputation based on Michigan Imputation Server, and post-imputation QC
# Author: Tao Wang

# Step 1, run data preparation script, HRC-1000G-check-bim.pl
# 
perl /home/tw83/twang/AMP/tools/HRC-1000G-check-bim.pl -b DATA.QC.bim -f DATA.QC.frq -r /home/tw83/twang/AMP/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

##########################
# VCFQC.sh
#!/bin/bash
# remove SNPs with imputation R2 < 0.3
# remove SNPs with MAF < 0.05
# remove SNPs with more than 2 alleles
[[ ! -d QC ]] && mkdir QC
for (( i = 2; i <= 22; i++ )); do
        zcat chr${i}.info.gz | awk 'NR!= 1 && $7>=0.3 {print $1}' > QC/chr${i}.R2.txt # extract snps with r2 >= 0.3      
        vcftools --gzvcf chr${i}.dose.vcf.gz --snps QC/chr${i}.R2.txt --maf 0.05 --min-alleles 2 --max-alleles 2 --recode  --stdout| gzip -c > QC/chr${i}.dos.postQC.vcf.gz
done