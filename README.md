
# Readme
A Tool.config file is necessary to run this tool!

This tool is mostly wrote in R3.5 and Python3.7, please make sure it can run on your enviroment.

Before use this QC tool, run a script to install dependence by running:
```
Rscript dependence.R
```


The Tool.config file is actually a *json* with one object “config", which include two object "transcriptome" and "genotyping". In both transcriptome and genotyping, the first key is **"usage"** which represent if this part is used or not. For example, if the usage in transcriptome is "FALSE" and in genotyping is "TRUE" then the tool will only run the genotyping part, all parameter in transcriptome is ignored.

The details of each parameters are listed below.
## Transcriptome

### Create TPM
If **usage** in transcriptome is TRUE *(or true, T, t)*, you need to select input data type including **fastq, bam, readcounts and TPM(FPKM/RPKM) matrix** by assigning the corresponding values. Note that if the input is a readcounts file, a gene length file is also needed.

If input files are fastq or bam, one(or two if fastq is paried-end)file(s) only represent one individual, you have to put all files in one directory and set the parameter *fastq* or  *bam* in the config file as the directory.

*e.g. "fastq":"/home/Document/experiment/" means all fastq files are in this directory. So does the "bam"*

After determind fastq or bam, the tool will count TPM(Transcripts per million reads) with <a href="https://github.com/deweylab/RSEM" target="_blank">**RSEM**</a> tool, you may need to install it first and confirm the value of each key in object "RSRM", the details of them are listed:

| parameter for RSEM| necessity(Y/N)|  remark  |
| --------   | -----:  | :----:    |
|rsem_path|  Y  |The directory of RSEM|
|pnum|N|number of process for RSEM|
|paired-end|Y|TRUE or FALSE|
|reffile|necessary|The name of the reference used|
|mapping software|necessary|Choose a tool to align reads including bowtie/bowtie2/star|
|bowtie software path|necessary|Directory of mapping software chosen|
|output_genome_bam|optional|alignments in genomic coordinates instead in transcript|
|estimate-rspd|optional|Set this option if you want to estimate the read start position distribution (RSPD) from data. Otherwise, RSEM will use a uniform RSPD|
|append_names|optional|tells RSEM to append gene_name/transcript_name to the result files|

After RSEM, each individuals will have a .result file as result, which contain a colume "TPM", then the tool will read all the .result files and create a "TPM_matrix.csv" for later analysis.

If input is readcount and gene length file, the tool will call readcounts2TPM.R and create a "TPM_matrix.csv".

If input is a TPM matrix, or we have already create it from other format, we can do the next step.

### transcriptomeQC
The following process is base on a TPM file.
In this part, we will do quality control and some normalizatation work base on TPM data.

#### Three statistics to identify outliers

##### RLE

The RLE (Relative Log Expression) analysis is based on the hypothesis that in a sample, the expression of most genes should be at similar level, and only a small portion of genes show very high or low expression levels. Suppose the gene expression matrix *G* with genes on the rows and samples on the columns, RLE analysis ﬁrst calculate the median expression value across samples for every gene. Then eachexpression valueinthe matrix issubtractedbythe median expression value of the corresponding gene. For each sample (column), the distribution of residuals of expression values should be centered close to zero, and has small variations. The RLE plot aligns the expression distributions of all samples, in the form of boxplots, side by side in increasing order of the interquartile range (IQR). The samples locate at the right most side are more likely to be outliers


##### Hierarchical clustering

We use(1 - spearman's correlation) as distances between samples and assume samples should be homogeneous, so all samples should have short distances between others and clustered homogeneously. So samples far away others may be outliers. So we first use distances which is 1-spearmen correlation to hierarchically clustered, than use top 100 gene sort by variance to calculate Mahalanobis Distance, then a chi2 p-value will be calculated based on mahalanobis distance. Then clusters with ≥ 60\% samples with Bonferroni-corrected p-values < 0.05 were marked as outliers.

##### D-statistic

The D-statistic of each sample is the median Spearman’s correlations with other samples, which reﬂects the distance to other samples. And the outlier is far from the peak of D-statistics distribution

After this three steps, we will do quantile normalization, if users have provided coveriance information, we can do gender test with sex information. We adopt combat function implemented in in SVA package to adjust the batch effects and the fsva function also in SVA package to adjust other known covariants and latent covariants.

| parameter        | necessity(Y/N)   |  remark  |
| --------   | -----:  | :----:  |
|RLEFilterPercent|Y|sample with higher RLE may be outliers|
|DSFilterPercent|Y|D-statistic|
|pvalues.cut|Y|cutoff p-value|
|cluster_level|Y|The max level of hcluster when finding outlier|
|coveriance|Y|TRUE or FALSE|
|coveriance_file|N|coveriance information file, necessary when coveriance is TRUE|

## Genotyping

When "usage":"TRUE", the genotyping part is functional. The tools is use plink to process genotype data, if your data is in VCF fromat you can also use plink to change it into plink format. However, only plink 1.9+ can change VCF to plink so make sure you have install plink1.9+ or you may have to change VCF to plink yourself. The tool also use *smartpca* to run PCA analysis, please install it first.

| parameter        | necessity(Y/N)   |  remark  |
| --------   | -----:  | :----:  |
|plink|N|Rawdata in plink format|
|VCF|N|Rawdata in VCF format(plink and VCF should have one and only one value)|
|workdir|N|Set the work directory|
|genotyping_rate|Y|the cutoff of genotyping rate|
|call_rate|Y|the cutoff of call rate|
|HWE|Y|HWE Test's cutoff p-value|
|MAF|Y|Minor Allele Frequency|
|F_outlier_n_sd|Y|determind F outlier with fmean $\pm$F_outlier_n_sd\*fsd|
|population|N|Population stratification| 

#### PCA
| parameter for smartpca| remark  |
| ------ | :----:  |
|m|default(5)maximum number of outlier removal iterations.To turn off outlier removal, set -m 0.|
|k|default(10)number of principal components to output|
|t|default(10)number of principal components along which to remove outliers during each outlier removal iteration|

#### Imputation
If the usage of Imputation is TRUE, the tool will use <a href="https://www.well.ox.ac.uk/~wrayner/tools/" target="_blank">**McCarthy's tool**</a> to do prepreparing checking. The tool can use three kinds of reference panel: HRC, 1000G and CAAPA, use the parameter in config file to select and determind the reference panel.

## Run Demo data
To show how this QCTool works, we use the E-GEUV-1's data as a demo and run the QC pipeline on it. In /Sample/fastq/ there is an *run_demo.sh*, which can automatically download three sample's RNA-seq data in fastq format in to /fastq/ and download reference to /ref/, than create reference files and use RSEM to calculate the expression data for each sample. The tool will create a */TranscriptomeWorkplace" directory and place the result including mapped bam file, which can also be used to create a TPM matrix, and each sample's expression data. After calculate expression for each sample, the results will be distilled and generate a TPM_matrix.csv file under /Sample/fastq/ and read for downstream analysis.

However, three samples are too less to run the downstream pipeline. To show the whole function of the tool, we use summary of alignment (read count) as another demo which is placed in /Sample/readcount/, by running *run_demo.sh*, a read count file will be download and calculate expression levels in tpm with a gene length file. it will also create a TPM_matrix.csv file and can be used for downstream analysis. The result will be stored in *expression.postQC.xls* and generate a *Report.html* to show the details of preprocessing under the working directory.

As for genotype data, we also create a demo in /Sample/VCF/, by running *run_demo.sh* it will download a <a href="https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GEUVADIS.chr1.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz" target="_blank">genotype data</a> in VCF format and preprocess it with QC tool. It will generate PLINK format data for each step of preprocessing and create a table in markdown named *GeneReport.md* containing details of each step. Here is an example:


|step|SNPs removed|SNP pass|sample removed|sample pass|
| --------   | :-----:  | :----:    | :-----:| :---: |
|genotyping rate|84133|409726|0|465|
|call rate|0|409726|0|465|
|sex check|0|409726|0|465|
|HWE test|5000|404726|0|465|
|Mishap|0|404726|0|465|
|MAF|239983|164743|0|465|
|F-outlier and IBD|0|164743|140|325|
|PCA outliers|0|164743|0|325|



# eQTLQC
