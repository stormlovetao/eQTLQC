# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import pandas as pd
import sys
import json
import re


config_file = 'Tool.config'
with open('Tool.config','r') as f:
    pardata = json.load(f)

wkdir = os.getcwd()
wkdir = os.path.join(wkdir,'TranscriptomeWorkplace/')
isExist = os.path.exists(wkdir)
if not isExist:
    os.mkdir(wkdir)

if (pardata['config']['transcriptome']['usage'].lower() in ['true','t']):
    #fastq
    
    
    def createTPMfile(filenames):
        firstSample = True
        for file in filenames:
            sampleTable = pd.read_csv(wkdir+file+'.genes.results',sep='\t')
            if firstSample:
                tpm_table = sampleTable[{'gene_id','TPM'}]
                tpm_table['gene_id'] = tpm_table['gene_id'].map(lambda x:x.split('_')[0])
                tpm_table.rename(columns={'TPM':file},inplace = True)
                firstSample = False
            else:
                tpm_table[file] = sampleTable['TPM']
        tpm_table.to_csv('TPM_matrix.csv',index = False)
    
    
    
    if (pardata['config']['transcriptome'].get('fastq','')!='' ):
        rsem_path = pardata['config']['transcriptome']['RSEM']['rsem_path']
        pnum = pardata['config']['transcriptome']['RSEM']['pnum']
        reffile = pardata['config']['transcriptome']['RSEM']['reffile']
        mapping_software = pardata['config']['transcriptome']['RSEM']['mapping software']
        mapping_software_path = pardata['config']['transcriptome']['RSEM']['mapping software path']
        est = ' --estimate-rspd' if (pardata['config']['transcriptome']['RSEM']['estimate-rspd'].lower() in ['true','t']) else ''
        outgenomebam = ' --output-genome-bam' if (pardata['config']['transcriptome']['RSEM']['output_genome_bam'].lower() in ['true','t']) else ''
        appendnames = ' --append-names' if (pardata['config']['transcriptome']['RSEM']['append_names'].lower() in ['true','t']) else ''
        pairedend = ' --paired-end' if (pardata['config']['transcriptome']['RSEM']['paired-end'].lower() in ['true','t']) else ''
        
        Options = '-p '+pnum + pairedend +' --'+mapping_software+ ' --' + mapping_software + '-path '+ mapping_software_path + est + appendnames + outgenomebam
        print('Options ='+Options)

        inputdir = pardata['config']['transcriptome']['fastq']
        filenames = os.listdir(inputdir)
        filelist = list()
        if pairedend == ' --paired-end':
            for file in filenames:
                x = file.split('_')
                if x[0] not in filelist:
                    filelist.append(x[0])
            filenames = filelist
            for file in filenames:
                command = rsem_path + '/rsem-calculate-expression ' + Options +' --paired-end '+inputdir+file+'_1.fastq '+inputdir+file+'_2.fastq '+reffile+' '+wkdir+file
                print('command:') 
                print(command)
                os.system(command)
        else:
            for file in filenames:
                command = rsem_path + '/rsem-calculate-expression ' + Options +' '+inputdir+file+'_1.fastq '+reffile+' '+wkdir+file
                print('command:') 
                print(command)
                os.system(command)
        
        createTPMfile(filenames)
                
    #bam
    elif (pardata['config']['transcriptome'].get('bam','') != ''):
        rsem_path = pardata['config']['transcriptome']['RSEM']['rsem_path']
        pnum = pardata['config']['transcriptome']['RSEM']['pnum']
        reffile = pardata['config']['transcriptome']['RSEM']['reffile']
        est = ' --estimate-rspd' if (pardata['config']['transcriptome']['RSEM']['estimate-rspd'].lower() in ['true','t']) else ''
        outgenomebam = ' --output-genome-bam' if (pardata['config']['transcriptome']['RSEM']['output_genome_bam'].lower() in ['true','t']) else ''
        appendnames = ' --append-names' if (pardata['config']['transcriptome']['RSEM']['append_names'].lower() in ['true','t']) else ''
        pairedend = ' --paired-end' if (pardata['config']['transcriptome']['RSEM']['paired-end'].lower() in ['true','t']) else ''
        
        Options = '-p '+pnum + ' --bam'  + pairedend + est + appendnames + outgenomebam
        print('Options ='+Options)
        
        inputdir = pardata['config']['transcriptome']['bam']
        filenames = os.listdir(inputdir)
        for file in filenames:
            command = rsem_path + '/rsem-calculate-expression ' + Options + ' '+inputdir+file +' '+ reffile+' '+wkdir+file.replace('.bam','')
            print('command:') 
            print(command)
            os.system(command)
            
        createTPMfile(filenames)
            
    #readcounts
    elif (pardata['config']['transcriptome'].get('readcounts','') != ''):
        countsfile = pardata['config']['transcriptome'].get('readcounts','')
        genelengthfile = pardata['config']['transcriptome']['gene_length_file']
        command = 'Rscript readcounts2TPM.R'
        os.system(command)
        
    #TPM
    elif (pardata['config']['transcriptome'].get('TPM','') != ''):
        command = 'cp '+pardata['config']['transcriptome'].get('TPM','')+' TPM_matrix.csv'
        os.system(command)
        
    else:
        print('Error: did not send input')
            

    command='R -e \"rmarkdown::render(\'report.Rmd\')\"'
    os.system(command)

if(pardata['config']['genotyping']['usage'].lower() in ['true','t']):
    if(pardata['config']['genotyping'].get('plink','')!=''):
        Rawdata = pardata['config']['genotyping'].get('plink','')
    elif(pardata['config']['genotyping'].get('VCF','')!=''):
        command = "plink --vcf "+pardata['config']['genotyping']['VCF']+" --allow-extra-chr --make-bed --out VCF2plink"
        os.system(command)
        Rawdata = "VCF2plink"
    workdir = pardata['config']['genotyping']['workdir']
    PCA_m = pardata['config']['genotyping']['smartpca']['m']
    PCA_k = pardata['config']['genotyping']['smartpca']['k']
    PCA_t = pardata['config']['genotyping']['smartpca']['t']
    initdir = os.getcwd()
    os.chdir(workdir)
    GeneReport = open("GeneReport.md",'wt')
    #step1 genotype rate
    print("step1")
    command = "plink --bfile "+Rawdata+" --geno "+pardata['config']['genotyping']['genotyping_rate']+" --make-bed --out DATA.geno > step1.txt"
    os.system(command)
    with open('step1.txt','r') as f:
        logfile = f.read()
    matchObj = re.compile(r'\d+'+' markers to be included from.+')
    match = matchObj.search(logfile).group()
    SNPs_pass = match.split()[0]
    matchObj = re.compile(r'\d+'+' individuals read from.+')
    match = matchObj.search(logfile).group()
    sample_pass = match.split()[0]
    matchObj = re.compile(r'\d+'+' SNPs failed missingness test.+')
    match = matchObj.search(logfile).group()
    SNPs_removed = match.split()[0]
    SNPs_pass = str(int(SNPs_pass) - int(SNPs_removed))
    print("|step|SNPs removed|SNP pass|sample removed|sample pass|" ,file = GeneReport)
    print("| --------   | :-----:  | :----:    | :-----:| :---: |" ,file = GeneReport)
    print("|genotyping rate|%s|%s|%s|%s|"%(SNPs_removed,SNPs_pass,'0',sample_pass) ,file = GeneReport)

    #step2 call rate
    print("step2")
    command = "plink --bfile DATA.geno --mind "+pardata['config']['genotyping']['call_rate']+" --make-bed --out DATA.filterSubject > step2.txt"
    os.system(command)
    with open('step2.txt','r') as f:
        logfile = f.read()
    matchObj = re.compile('.+'+' individuals removed for low genotyping.+')
    match = matchObj.search(logfile).group()
    tmp = match.split()
    sample_removed = tmp[0]
    sample_pass = str(int(sample_pass) - int(sample_removed))
    print("|call rate|%s|%s|%s|%s|"%('0',SNPs_pass,sample_removed,sample_pass) ,file = GeneReport)

    #step3 sex check
    print("step3")
    command = "plin --bfile DATA.filterSubject --chr X --make-bed --out DATA.chrX"
    os.system(command)
    ChrXSNP = open('DATA.chrX.bim')
    ChrXSNPnum = len(ChrXSNP.readlines())
    if ChrXSNPnum >= 10000:
        command = "plink --bfile DATA.filterSubject --check-sex --out DATA.filterSubject.Sex > step3.txt"
        os.system(command)
        command = "awk \'BEGIN{OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6}\' DATA.filterSubject.Sex.sexcheck > DATA.filterSubject.Sex.txt"
        os.system(command)
    ## Check DATA.filterSubject.Sex.txt and see if any sex mismatch (problem), and generated a file DATA.Sex.mismatch.txt
        command = "awk \'$5 != \"OK\"\' DATA.filterSubject.Sex.txt | awk \'NR!=1{print $1,$2}\' > DATA.Sex.mismatch.txt"
        os.system(command)
        command = "plink --bfile DATA.filterSubject --remove DATA.Sex.mismatch.txt --make-bed --out DATA.filterSubject.filterSexMis >> step3.txt"
        os.system(command)
        with open('step3.txt','r') as f:
            logfile = f.read()
        matchObj = re.compile('.+'+' individuals removed')
        match = matchObj.search(logfile).group()
        sample_removed = match.split()[0]
        sample_pass = str(int(sample_pass) - int(sample_removed))
        print("|sex check|%s|%s|%s|%s|"%('0',SNPs_pass,sample_removed,sample_pass) ,file = GeneReport)
    #step4 heterozygous hapliod
        print("step4")
        command = "plink --bfile DATA.filterSubject.filterSexMis --set-hh-missing --make-bed --out DATA.filterSubject.filterhh > step4.txt"
        os.system(command)
    else:
        print("dont have enough SNPs on chr X, cannot do check sex")
        print("|sex check|%s|%s|%s|%s|"%('0',SNPs_pass,'0',sample_pass) ,file = GeneReport)
        print("step4")
        command = "plink --bfile DATA.filterSubject --set-hh-missing --make-bed --out DATA.filterSubject.filterhh > step4.txt"
        os.system(command)

    #step5 HWETest
    print("step5")
    command = "plink --bfile DATA.filterSubject.filterhh --hwe "+pardata['config']['genotyping']['HWE'] + " --out DATA.Hardy --make-bed > step5.txt"
    os.system(command)
    with open('step5.txt','r') as f:
        logfile = f.read()
    matchObj = re.compile('.+'+'markers to be excluded based on HWE test.+')
    match = matchObj.search(logfile).group()
    SNPs_removed = match.split()[0]
    SNPs_pass = str(int(SNPs_pass) - int(SNPs_removed))
    print("|HWE test|%s|%s|%s|%s|"%(SNPs_removed,SNPs_pass,'0',sample_pass) ,file = GeneReport)

    #step6  Remove Test-mishap SNPs
    print("step6")
    command = "plink --bfile DATA.Hardy --test-mishap --out DATA.Mis"
    os.system(command)
    command = "awk \'BEGIN{OFS=\"\t\"}{if($8 < 1e-9){print $1,$2,$3,$4,$5,$6,$7,$8,$9 }}\' DATA.Mis.missing.hap |  grep \'HETERO\' | cut -f1 | sort | uniq > DATA.mishap.txt"
    os.system(command)
    mishap = open('DATA.mishap.txt')
    SNPs_removed = str(len(mishap.readlines()))
    mishap.close()
    SNPs_pass = str(int(SNPs_pass) - int(SNPs_removed))
    print("|Mishap|%s|%s|%s|%s|"%(SNPs_removed,SNPs_pass,'0',sample_pass) ,file = GeneReport)

    command = "plink --bfile DATA.Hardy --exclude DATA.mishap.txt --make-bed --out DATA.Mis > step6.txt"
    os.system(command)

    #step7 MAF(Minor Allele Frequency)
    print("step7")
    command = "plink --bfile DATA.Mis --maf "+pardata['config']['genotyping']['MAF']+" --make-bed --out DATA.MAF > step7.txt"
    os.system(command)
    with open('step7.txt','r') as f:
        logfile = f.read()
    matchObj = re.compile('.+'+'SNPs failed frequency test.+')
    match = matchObj.search(logfile).group()
    SNPs_removed = match.split()[0]
    SNPs_pass = str(int(SNPs_pass) - int(SNPs_removed))
    print("|MAF|%s|%s|%s|%s|"%(SNPs_removed,SNPs_pass,'0',sample_pass) ,file = GeneReport)

    #step8,step9
    command = "plink --bfile DATA.MAF --indep 50 5 2 --out DATA.MAF"
    os.system(command)
    command = "plink --bfile DATA.MAF --extract DATA.MAF.prune.in --make-bed --out DATA.MAF.subset"
    os.system(command)
    #step 8 remove F outliers
    print("step8")
    command = "plink --bfile DATA.MAF.subset --het --out DATA.Het"
    os.system(command)
    command = "awk \'BEGIN{OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6}' DATA.Het.het > DATA.Het.het1"
    os.system(command)

    data = pd.read_csv('DATA.Het.het1', sep='\t')
    sd = data['F'].std()
    Fmean = data['F'].mean()
    n_sd = pardata['config']['genotyping']['F_outlier_n_sd']
    n_sd = int(n_sd)
    outlier = data[data['F'] < Fmean - n_sd*sd]
    outlier = outlier.append(data[data['F'] > Fmean + n_sd*sd])
    outlier[['FID','IID']].to_csv('DATA.QC.outlier.txt',sep=' ',index = False,header = False)

    # Step9 IBS/IBD Filtering subjects
    print("step9")
    command = "plink --bfile DATA.MAF.subset --genome --out DATA.MAF.subset.relatedness"
    os.system(command)
    command = "awk \'{if($10>0.1875){print $1,$2}}\' DATA.MAF.subset.relatedness.genome | sed 1d >> DATA.QC.outlier.txt"
    os.system(command)
    command = "awk \'{if($10>0.98){print $3,$4}}\' DATA.MAF.subset.relatedness.genome | sed 1d >> DATA.QC.outlier.txt"
    os.system(command)
    command = "plink --bfile DATA.MAF --remove DATA.QC.outlier.txt --make-bed --out DATA.QC > step9.txt"
    os.system(command)
    with open('step9.txt','r') as f:
        logfile = f.read()
    matchObj = re.compile('.+'+'individuals removed.+')
    match = matchObj.search(logfile).group()
    sample_removed = match.split()[0]
    sample_pass = str(int(sample_pass) - int(sample_removed))
    print("|F-outlier and IBD|%s|%s|%s|%s|"%('0',SNPs_pass,sample_removed,sample_pass) ,file = GeneReport)

    #PCA
    print("doing PCA analysis")
    command = "plink --bfile DATA.QC --indep 50 5 2 --out DATA.QC.prune"
    os.system(command)
    command = "plink --bfile DATA.QC --extract DATA.QC.prune.prune.in --make-bed --out DATA.QC.pruned"
    os.system(command)
    command = "cp DATA.QC.pruned.bim DATA.QC.pruned.pedsnp"
    os.system(command)
    command = "cp DATA.QC.pruned.fam DATA.QC.pruned.pedind"
    os.system(command)
    command = "smartpca -i DATA.QC.pruned.bed -a DATA.QC.pruned.pedsnp -b DATA.QC.pruned.pedind -o DATA.QC.pruned.pca -p DATA.QC.pruned.plot -e DATA.QC.pruned.eval -l DATA.QC.smartpca.log -m "+PCA_m+" -k "+PCA_k+" -t "+PCA_t
    os.system(command)
    command = "awk \'/REMOVED/ {print $3}\' DATA.QC.smartpca.log | awk \'BEGIN{FS=\":\";OFS=\"\t\"}{print $1,$2}\' > DATA.QC.Pop.outlier.txt"
    os.system(command)
    command = "plink --bfile DATA.QC --remove DATA.QC.Pop.outlier.txt --make-bed --out DATA.QC.Pass"
    os.system(command)
    PCAoutlier = open('DATA.QC.Pop.outlier.txt')
    sample_removed = len(PCAoutlier.readlines())
    PCAoutlier.close()
    sample_pass = str(int(sample_pass) - int(sample_removed))
    print("|PCA outliers|%s|%s|%s|%s|"%('0',SNPs_pass,sample_removed,sample_pass) ,file = GeneReport)

    #More
    if 'population' in pardata['config']['genotyping'].keys():
        print("merge data with hapmap")
        command = "awk \'BEGIN{OFS=\"\t\"} NR == FNR {array[$2\"\t\"$5\"\t\"$6]=1} NR > FNR {if(array[$2\"\t\"$5\"\t\"$6]==1){print $2}}\' "+pardata['config']['genotyping']['population']+".bim DATA.QC.bim > ShareSNP.txt"
        os.system(command)
        command = "plink --bfile DATA.QC --extract ShareSNP.txt --make-bed --out DATA.QC.Share"
        os.system(command)
        command = "plink --bfile "+pardata['config']['genotyping']['population']+" --extract ShareSNP.txt --make-bed --out Hapmap3.QC.Share"
        os.system(command)
        command = "plink --bfile DATA.QC.Share --bmerge Hapmap3.QC.Share.bed Hapmap3.QC.Share.bim Hapmap3.QC.Share.fam --make-bed --out DATA.Hapmap3"
        os.system(command)
        command = "plink --bfile DATA.Hapmap3 --indep 50 5 2 --out DATA.Hapmap3.prune"
        os.system(command)
        command = "plink --bfile DATA.Hapmap3 --extract DATA.Hapmap3.prune.prune.in --make-bed --out DATA.Hapmap3.pruned"
        os.system(command)
        command = "cp DATA.Hapmap3.pruned.bim DATA.Hapmap3.pruned.pedsnp"
        os.system(command)
        command = "cp DATA.Hapmap3.pruned.fam DATA.Hapmap3.pruned.pedind"
        os.system(command)
        command = "awk \'{print $1,$2,$3,$4,$5,-9}\' DATA.Hapmap3.pruned.pedind > DATA.Hapmap3.pruned.pedind9"
        os.system(command)
        command = "smartpca -i DATA.Hapmap3.pruned.bed -a DATA.Hapmap3.pruned.pedsnp -b DATA.Hapmap3.pruned.pedind9 -o DATA.Hapmap3.pruned.pca -p DATA.Hapmap3.pruned.plot -e DATA.Hapmap3.pruned.eval -l DATA.Hapmap3.smartpca.log -m "+PCA_m+" -k "+PCA_k+" -t "+PCA_t
        os.system(command)
        command = "Rscript "+initdir+"/ggplot.population.R"
        os.system(command)
        
    #Imputation
    if (pardata['config']['genotyping']['Imputation']['usage'].lower() in ['true','t']):
        command = "plink --bfile DATA.QC.Pass --freq --out DATA.QC.Pass"
        os.system(command)
        if(pardata['config']['genotyping']['Imputation'].get('HRC','') != ''):
            command = "perl HRC-1000G-check-bim.pl -b DATA.QC.Pass.bim -f DATA.QC.Pass.frq -r "+pardata['config']['genotyping']['Imputation'].get('HRC','')+" -h"
            os.system(command)
        elif(pardata['config']['genotyping']['Imputation'].get('1000G','') != ''):
            command = "perl HRC-1000G-check-bim.pl -b DATA.QC.Pass.bim -f DATA.QC.Pass.frq -r "+pardata['config']['genotyping']['Imputation'].get('1000G','')+" -g"
            os.system(command)
        elif(pardata['config']['genotyping']['Imputation'].get('CAAPA','') != ''):
            command = "perl HRC-1000G-check-bim.pl -b DATA.QC.Pass.bim -f DATA.QC.Pass.frq -r "+pardata['config']['genotyping']['Imputation'].get('CAAPA','')+" -h"
            os.system(command)
        else:pass
        os.system("./Run-plink.sh")
