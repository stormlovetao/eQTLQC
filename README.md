# Introduction
Advances in the expression quantitative trait loci (eQTL) studies have provided valuable insights into the mechanism of diseases and traits-associated genetic variants. However, it remains computationally challenging to evaluate and control the data quality required in eQTL analysis for researcher with a limited background in programming. There is an urgent need to develop a powerful and user-friendly tool to automatically process the raw dataset and perform the eQTL mapping afterwards. In this work, we present a pipeline for eQTL analysis, termed as eQTLQC, with automated data preprocessing of both genotype data and gene expression data, especially RNA-seq based. Our pipeline provides comprehensive quality control and normalization approaches. Multiple automated techniques are used to concatenate the pipeline to reduce manual interventions. And a configuration file is provided to adjust parameters and control the processing logic. The source code, demo data and instructions are freely available <a href="https://github.com/stormlovetao/eQTLQC" target="_blank">**here**</a>.


# Preliminaries

This tool is written in R 3.5 and Python 3.7, please double check your computing environment.

The "Tool.config" file is in JSON format, and required to configurate parameters and control the processes.

Before using this tool, run the following R script to install dependencies:
```
Rscript dependence.R
```
# Configuration

The configuration file, named "Tool.config" and in JSON format, is used to setup parameters and control the pipeline.
JSON is a human-readable text format to store "attribute-value" pairs. More details introducing JOSN could be found at <a href="https://en.wikipedia.org/wiki/JSON" target="_blank">**wikipedia**</a>. 
The Tool.config file consists of two main sub-objects: "transcriptome" and "genotyping". The key named **"usage"** is used as switch to turn on or off the processing sections of gene expression data ("transcriptome") and genotype data ("genotyping"). For example, if the usage in "transcriptome" is set to "FALSE" and in "genotyping" is set to "TRUE", then the tool will only run the genotyping processing part, all parameters in transcriptome will be ignored.

The details of parameters are introduced as follows.

## Transcriptome data preprocessing

If **usage** in transcriptome is TRUE *(or true, T, t)*, one and only one of the following formats is required: **Fastq**, **BAM**, **Readcounts** or **normalized metrics(TPM/FPKM/RPKM)**.  Note that if the input is a readcounts file, additional gene length file is also needed.

### FASTQ/BAM inputs
If input sources are in FASTQ or BAM format, FASTQ or BAM files should be under same folder path. For example:
"fastq":"/home/rawdata/fastq/" or "bam":"/home/rawdata/bam/".
For pair-ended reads, filenames of FASTQ should be ended with "_1.fastq" and "_2.fastq".

eQTLQC will process the raw inputs into the TPM (Transcripts per million reads) metrics using <a href="https://github.com/deweylab/RSEM" target="_blank">**RSEM**</a>. Parameters of RSEM could be setup in configuration file by following keys:

| Keys| Mandatory/optional| Defaults| Remarks  |
| --------   | -----:     | --------| :----    |
|rsem_path   |  Mandatory  | NA |string, the directory of RSEM installation|
|pnum        |  optional   |  8 |string, number of threads to use|
|paired-end  |  Mandatory  | TRUE| boolean, TRUE or FALSE|
|reffile     |  Mandatory  | NA | string, the name of the reference used|
|mapping_software|Mandatory|bowtie2| string, choose alignment software: "bowtie","bowtie2" or "star"|
|mapping_software_path|optional| NA|string, Path of mapping software chosen. If not given, the path to the executables is assumed to be in the user's PATH environment variable|
|output_genome_bam|optional|TRUE|boolean, alignments in genomic coordinates instead in transcript|
|estimate-rspd|optional|TRUE|boolean, set this option if you want to estimate the read start position distribution (RSPD) from data. Otherwise, RSEM will use a uniform RSPD|
|append_names|optional|boolean, append gene_name/transcript_name to the result files|
|add_parameters|optional||


### Readcounts/Normalized metrics inputs
If readcounts table is provided for key "readcounts", eQTLQC will scale it into TPM values. A gene length file is also requested from user with key "gene_length_file". See gene length demo file format: Sample/readcount/gencode.v24.annotation.gtf.gene_length_unionexon. 

If users do not want to transform read counts into TPM, users could directly assign readcounts file path to the key "metrics".

If normalized metrics are provided for the key "metrics", such as TPM/RPKM/FPKM or even readcounts, eQTLQC will directly turn into downstream preprocessing procedures.

ReadCount/Metrics file format requirements: The input file should has sample ID as column name; the first column is gene ID (If ENGSID is used, the version behind '.' will be ignored), and the following columns are samples; Fields in each row should be separated by tab separator.

### Quality control on gene expression profiles
Quality control procedures will be applied on the "metrics" user provided or eQTLQC generated. 

The following three methods will be used to identify outliers with problematic expression profiles.

#### RLE

The RLE (Relative Log Expression) analysis is based on the hypothesis that in a sample, the expression of most genes should be at similar level, and only a small portion of genes show very high or low expression levels. Suppose the gene expression matrix *G* with genes on the rows and samples on the columns, RLE analysis ﬁrst calculate the median expression value across samples for every gene. Then each expression value in the matrix is subtracted by the median expression value of the corresponding gene. For each sample (column), the distribution of residuals of expression values should be centered close to zero, and has small variations. The RLE plot aligns the expression distributions of all samples, in the form of boxplots, side by side in increasing order of the interquartile range (IQR). The samples locate at the right most side are more likely to be outliers. In default, we set the right most 5% samples as candidate outliers.


#### Hierarchical clustering

Hierarchical clustering is also a widely used approach to exclude sample outliers. Similarity between each pair of samples is first measured by metrics such as Pearson's correlation coefficient or Spearman's correlation coefficient. The sample distance matrix, obtained by one minus similarity scores, is used to perform hierarchical clustering. Usually, samples with problematic expression profiles will be far from other samples with normal expression profiles in the clustering dendrogram.

We measure the Mahalanobis distances between each sample and the distribution of all samples. The chi-squared P value is calculated for each sample in each cluster. Clusters with <= *cluster_percent* (60% in default) samples with Bonferroni-corrected P-values <= *pvalues_cut* (0.05 in default) are marked as outlying clusters, and all samples included will be marked as candidate outliers.

#### D-statistic

The D-statistic of each sample is defined as the median Spearman's correlation coefficients with other samples, which reﬂects the distance to other samples. The outliers will be likely far from the peak of D-statistics distribution.  We consider the left most 5\% samples as candidate outliers.

And conducting the analysis of above analysis, eQTLQC considers the intersection of candidate outliers reported by the three methods as final set of sample outliers. 

### Normalization and covariates adjustment
In eQTLQC, we use SVA to adjust the known and hidden covariates. To be specified, we use the *combat* function in SVA to adjust the batch effects, and use the *fsva* function to adjust other known covariates and latent covariates.

Covariates file format requirements: columns should be sperated by tab. The first column should be sample's ID which should be consistant with sample ID in expression data. "sex" and "batch" are two reserved covariate column names, which will be used in gender-mismatch check and adjusting batch effects respectively. Other covariates could be named arbitrarily. Sex, batch and other covariates in string will be treated as factors, and other covariates in numbers will be treated as numeric.


| Keys| Mandatory/optional| Defaults|  Remarks  |
| ----| -----              | -----   |  ---      |
|RLEFilterPercent|Mandatory| 0.05  |numeric, RLE distribution, right-most percentage|
|DSFilterPercent |Mandatory| 0.05  |numeric, D-statistic distribution, left-most percentage|
|pvalues_cut     |Mandatory| 0.05  |numeric, cutoff p-value|
|cluster_level   |Mandatory| 5     |numeric, The maximum level of hcluster when detecting outlying clusters|
|covariates      |Mandatory|FALSE  |boolean, covariates available or not |
|covariates_file |Optional  |NA     |string, covariates file path, necessary when covariates is TRUE|


## Genotype data preprocessing

If the key "usage" is set to "TRUE", the genotype preprocessing part will be activited. eQTLQC accepts genotype data in PLINK and VCF formats, and VCF files will be converted to PLINK formats for further processing. eQTLQC adopts PLINK for genotype QC, and <a href="https://www.cog-genomics.org/plink/" target="_blank">**Plink1.9+**</a> is required. To specify population stratification, eQTLQC adopts SmartPCA in <a href="https://www.hsph.harvard.edu/alkes-price/software/" target="_blank">**EIGENSOFT**</a>, which is also required.

| Keys| Mandatory/optional|  Remarks  |
| --------   | -----:  | :----:  |
|PLINK          |Optional|Raw data in PLINK format, mandatory if VCF is null|
|VCF            |Optional|Raw data in VCF format, mandatory if PLINK is null|
|work_dir       |Optional|Set the work directory|
|genotyping_rate|Mandatory|the cutoff of genotyping rate|
|call_rate      |Mandatory|the cutoff of call rate|
|HWE            |Mandatory|HWE Test's cutoff p-value|
|MAF            |Mandatory|Minor Allele Frequency|
|F_outlier_n_sd |Mandatory|determind F outlier with fmean $\pm$F_outlier_n_sd\*fsd|
|population     |Optional|Remove population outlier| 

#### PCA
| parameter for smartpca| remark  |
| ------ | :----:  |
|m|default(5)maximum number of outlier removal iterations.To turn off outlier removal, set -m 0.|
|k|default(10)number of principal components to output|
|t|default(10)number of principal components along which to remove outliers during each outlier removal iteration|

#### Imputation
If the usage of Imputation is TRUE, the tool will use <a href="https://www.well.ox.ac.uk/~wrayner/tools/" target="_blank">**McCarthy's tool**</a> to do prepreparing checking. The tool can use three kinds of reference panel: HRC, 1000G and CAAPA, use the parameter in config file to select and determind the reference panel.

## Run Demo data
To show how eQTLQC works, we use the E-GEUV-1's data as a demo and run the QC pipeline on it. In /Sample/fastq/ there is an *run_demo.sh*, which can automatically download three sample's RNA-seq data in fastq format in to /fastq/ and download reference to /ref/, than create reference files and use RSEM to calculate the expression data for each sample. The tool will create a */TranscriptomeWorkplace" directory and place the result including mapped bam file, which can also be used to create a TPM matrix, and each sample's expression data. After calculate expression for each sample, the results will be distilled and generate a TPM_matrix.csv file under /Sample/fastq/ and read for downstream analysis.

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
