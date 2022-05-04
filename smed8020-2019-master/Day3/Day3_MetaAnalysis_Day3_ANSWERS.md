# ANSWERS
# DAY 3 - Meta-analysis for quantitative trait using METAL

## FOCUS and LEARNING GOALS

> The aim for this session is to get familiar with running a genome wide association study meta-analysis which  enables researchers to  gather data from many studies and analyse them together. 

There are several motivations for meta-analysis. One is the ability to increase power to detect small effect sizes or rare variant effects by increasing the study sample size. Many methods for meta analysis rely on using summary statistics therefore rendering the need to share individual level data unnecessary. This makes it easier to share data for meta analysis as the summary statistics are not deemed sensitive information and are typically made publicly available when papers are published in peer-reviewed journals. Finally, meta-analysis across genetic ancestries is the most statistically robust approach rather than pooling all ancestries together in one GWAS as it results in little or no loss of efficiency (as compared to analysis of combined data-sets) and reduces population stratification.
Some [slides](MetaAnalysis.pdf) with extra informaiton may be helpful.  

**Suggested reading:**

* [Willer, C. J., Li, Y. & Abecasis, G. R. METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics 26, 2190–2191 (2010).](https://academic.oup.com/bioinformatics/article/26/17/2190/198154)  
* [Nielsen, J. et al. Loss-of-function genomic variants with impact on liver-related blood traits highlight potential therapeutic targets for cardiovascular disease. Biorxiv. (2019) ](https://www.biorxiv.org/content/10.1101/597377v1)  

[METAL](http://csg.sph.umich.edu/abecasis/metal/) was developed at the University of Michigan as a tool for meta-analysis of  genome-wide association analysis. For running METAL you can use either test statistics and standard error or p-values. For more info on METAL see the web-links in this document and the suggested readings paper. 

## TASK: Running a trans-ethnic meta-analysis using METAL  

Today you will run a meta-analysis to combine three studies using METAL. METAL has been pre-installed on our lab. Because of time restraints we have made a small data-set that will run within reasonable time. A separate set of files will therefore be used for plotting results. 
    
The phenotype in today's practical is low density lipoprotein (LDL) cholesterol and we will be using data from three large studies: HUNT, Biobank Japan and Global Lipids Genetics Consortium (GLGC).  

### Task outline  

1. Log into your lab
- SSH into the smed8020‐2019 lab: `ssh smed8020‐2019`  
- Once inside the lab, ssh to your own machine: `ssh ubuntu@smed8020‐2019‐nodeX` (password: ubuntu)
1. Gather summary statistics from GWAS for low density lipoprotein (LDL) cholesterol in three separate studies (HUNT, Global Lipids Genetics Consortium, and UK Biobank) and check details for the files.
2. Run a meta-analysis using METAL
3. Have a coffee or a biobreak or ask questions to the lecturers. Questions for consideration are in ****bold****. 
4. View meta-analysis results


### Important points to consider  

When running a meta analysis there are many issues that need to be addressed.
* the availability of summary statistics
* phenotype related questions such as:
  * what phenotypes are available?
  * what are the phenotypes based on (self-reported, Electronic Health Records, physician curated)?
  * how are the phenotypes constructed?
  * are they comparable to your defined phenotype?
* genotype related questions such as: 
  * differences in genotyping and imputation, some markers will be study specific
  * genome build
  * flipped markers
  * population stratification

### Instructions  
1. Organizing summary statistics  

Usually you would download publically available summary statistics from the internet to your local machine. For convience for this practical, the data can be downloaded from [here](https://ntnu.box.com/s/0veanic525k4m8m45x9rs4onh7snwpfw)

* The original summary statistics from Biobank Japan (BBJ) of LDL cholesterol in N=72,866 can be found [here](https://humandbs.biosciencedbc.jp/files/hum0014/hum0014_README_QTL_GWAS.html)  
`/mnt/scratch/day3/data/BBJ-LDL-preMeta.txt`  
The columns are CHR     POS38   SNPID   Allele1 Allele2 AC_Allele2      AF_Allele2      N       BETA    SE      p.value log10P

* The original summary statistics of joint analysis of metabochip and GWAS data for LDL cholesterol in N=89,138 from the Global Lipids Genetics Consortium (GLGC) can be found [here](http://csg.sph.umich.edu/willer/public/lipids2013/)  
`/mnt/scratch/day3/data/GLGC-LDL-preMeta.txt`  
The columns are SNP_hg18        SNP_hg19        rsid    A1      A2      beta    se      N       P-value Freq.A1.1000G.EUR

* The summary statistics of LDL cholesterol from the HUNT study in N=67,429.   
`/mnt/scratch/day3/data/HUNT-LDL-preMeta.txt`  
The columns are CHR     POS38   SNPID   Allele1 Allele2 AC_Allele2      AF_Allele2      N       BETA    SE      Tstat   p.value varT    varTstar
        log10P  AC_Allele1

2. Check your summary statistics to make sure they're ready for meta-analysis.

The human reference genome has been updated over the years and variants are given different coordinates in different versions. 
The latest human reference genome GRCh38 was released from the Genome Reference Consortium on 17 December 2013.  
The previous human reference genome (GRCh37) was the nineteenth version (hg19).  
The version before this was NCBI Build 36.1	released March 2006	(hg18). 
You can see more [here](https://genome.ucsc.edu/FAQ/FAQreleases.html#release1). hg19 is still widely used and people are slowly converting to hg38.  
****From the summary statistic headers, can you tell  what reference genome versions are used for each study?****  

It looks like BBJ and HUNT have SNP coordinates from hg38, but GLGC has summary statistics from hg18 and hg19. 
We must use [UCSC listOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert the hg19 coordinates to hg38 before meta-analysis. 

Create a .bed file file from GLGC-LDL-preMeta.txt using linux tools `awk` and `sed` in the following the command:
`awk ' NR > 1 {print $2"\t"$3"\t"$4"\t"$5}' GLGC-LDL-preMeta.txt | sed 's/:/\t/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$4}' > GLCG.hg19.bed`

Use the .bed file you made in liftOver. The liftover command requires 4 parameters in this order: 1) oldFile 2) map.chain 3) newFile 4) unMapped. Execute this command.
`liftOver GLCG.hg19.bed hg19ToHg38.over.chain GLGC.h38.bed GLGC.hg38.unmapped`

Look in GLGC.hg38.unmapped. ****Were there some markers that did not get converted from build 37 to build 38? Why do you think that is?****  

All of the regions were not converted because they were "Deleted in new". This means that there was no region in hg19 that aligned well to your region in hg18. 

Other error messages that can arise:
-"Partially deleted in new" means that there was a piece of an alignment that matched your hg18 region, but not enough to convert the region to hg19. 
-"Split in new", which means that your hg18 region was split up into different parts of hg19.

To avoid extensive file manipulation on your part, we already used this .bed file to make a new version of the GLGC results: `GLGC-LDL-hg38-preMeta.txt`. 
This file also has a header that is consistent with the other two files. You will use this in the meta-analysis.

The code to create the file is here:
`join -1 4 -2 2 <(sort -k 4 GLGC.h38.bed) <(sort -k 2 GLGC-LDL-preMeta.txt) | awk -v OFS='\t' '{$5=toupper($5);$9=toupper($9)}1' | awk '{print $0"\t"substr($2, 4)"\t"$1":"$9":"$5}'  | sed  '1i\CHRPOS\tchr\tstart\tPOS38\tAllele2\tCHRPOS37\trsid\ta2\tAllele1\tBETA\tSE\tN\tp.value\tAF_Allele2\tCHR\tSNPID' > GLGC-LDL-hg38-preMeta.txt`

****Are your SNPIDs across the files formatted in the same way?****
Yes. You can check this using the grep command by checking at least one SNPID in each set of results
`grep 'chr18:4440419:G:T' BBJ-LDL-preMeta.txt`
chr18:4440419:G:T #BBJ-LDL-preMeta.txt
chr18:4440419:G:T #GLGC-LDL-hg38-preMeta.txt
chr18:4440419:G:T #HUNT-LDL-preMeta.txt

****How many variants are in each of the files?****  
```
wc -l BBJ-LDL-preMeta.txt
wc -l HUNT-LDL-preMeta.txt
wc -l GLGC-LDL-hg38-preMeta.txt
```

6108162 BBJ-LDL-preMeta.txt    
25879010 HUNT-LDL-preMeta.txt    
2436641 GLGC-LDL-hg38-preMeta.txt    

The HUNT summary statistics are so large because imputation was done with the TOPMed imputation panel,  which allows for higher resolution imputation due to the large amount of sequencing samples which make up the reference panel. 
****What imputation panel was used for GLGC?**** HINT: Check the methods of the [paper](https://www.nature.com/articles/ng.2797).

GWAS data were imputed to the International HapMap project (https://www.nature.com/articles/nature09270?page=55#Sec2).

****How many genome wide significant results are in each of the input files?****  
```
awk '$12 < 5e-8 {print 0}' HUNT-LDL-preMeta.txt | wc -l   
#4591
awk '$11 < 5e-8 {print 0}' BBJ-LDL-preMeta.txt | wc -l   
#1998
awk '$13 < 5e-8 {print 0}' GLGC-LDL-hg38-preMeta.txt | wc -l
#3070
```

3. Running METAL   

The [Wiki page for METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation#Brief_Description)  may be useful.

Input files: 
* A text file for each study with results, summarized as a table. NB column separators must be specified. 
* A column with marker name, which should be consistent across studies 
* A column indicating the tested allele 
* A column indicating the other allele 
* If you are carrying out a sample size weighted analysis (based on p-values), you will also need: 
  * A column indicating the direction of effect for the tested allele 
  * A column indicating the corresponding p-value 
  * An optional column indicating the sample size (if the sample size varies by marker) 
* If you are carrying out a meta-analysis based on standard errors, you will need: 
  * A column indicating the estimated effect size for each marker 
  * A column indicating the standard error of this effect size estimate 
  * The header for each of these columns must be specified so that METAL knows how to interpret the data. 
 
A meta-analysis script, `LDL_metal.sh` has been created for you. You can run it with the following commands: 
1. Create a config file with the bash script `LDL_METAL.sh` by filling in the appropriate arguments instead of "file1",  "file2",  "file3" and using "LDL_METAL" as your output prefix.
`bash LDL_METAL.sh  file1 file2 file3 LDL_METAL > LDL_METAL.conf`    
i.e. `bash LDL_metal.sh HUNT-LDL-preMeta.txt GLGC-LDL-hg38-preMeta.txt BBJ-LDL-preMeta.txt LDL_METAL_META > LDL_METAL.conf`
2. Run metal with the config file (this will take less than 20 minutes)
`metal LDL_METAL.conf > LDL_METAL.log`  
If you would like to time your analysis you can use the time program.  
`/usr/bin/time -o test_time -v metal LDL_METAL.conf > LDL_METAL.log`

****What type of meta-analysis did you run (fixed or random effects? sample size or inverse variance based?) What is the difference?****  
Fixed effects using p-value and sample size. By default, METAL combines p-values across studies taking into account a study specific weight (typically, the sample size) and direction of effect. This behavior can be requested explicitly with the `SCHEME SAMPLESIZE` command.    
    
An alternative can be requested with the `SCHEME STDERR` command and weights effect size estimates using the inverse of the corresponding standard errors. To enable this option, you will also need to specify which of your input columns contains standard error information using the `STDERRLABEL` command (or `STDERR` for short). While standard error based weights are more common in the biostatistical literature, if you decide to use this approach, it is very important to ensure that effect size estimates (beta coefficients) and standard errors use the same units in all studies (i.e. make sure that the exact same trait was examined in each study and that the same transformations were applied). Inconsistent use of measurement units across studies is the most common cause of discrepancies between these two analysis strategies.

The difference between the fixed effects and random effects models is that fixed effects meta-analysis assumes that the genetic effects are the same across the different studies. Fixed effects models provide narrower confidence intervals and significantly lower P-values for the variants than random effects models.   
 
The random effects model assumes that the mean effect (of each SNP) in each study is different, with those means usually assumed to be chosen from a Gaussian distribution. The variance of that Gaussian distribution, and thus the amount of between-study heterogeneity, is estimated by the model.   

****Did you use genomic control? In what situations is it useful to use genomic control****  
No. METAL has the ability to apply a genomic control correction to all input files. METAL will estimate the inflation of the test statistic by comparing the median test statistic to that expected by chance, and then apply the genomic control correction to the p-values (for SAMPLESIZE weighted meta-analysis) or the standard error (for STDERR weighted meta-analysis). This should only be applied to files with whole genome data (i.e. should not be used for settings where results are only available for a candidate locus or a small number of SNPs selected for follow-up of GWAS results). Genomic control settings can be customized for each input file. We recommend applying genomic control correction to all input files that include genomewide data and, in addition, to the meta-analysis results. To apply genomic control to the meta-analysis results, just perform an initial meta-analysis and then load the initial set of results into METAL to get final, genomic control adjusted results.   

****What does it mean to set the minimum weight to 10,000?****   
METAL does not require that all input files report a result for every marker. Any available data is used. To restrict the output to only markers that have at least a specific number of individuals analysed (or weight), use a command like the following:   
`MINWEIGHT 10000`   
This will restrict the output to show only Markers with a total sample size of at least 10,000 individuals.   

****What is the difference between "ANALYZE" and "ANALYZE HETEROGENEITY"?****  
The METAL `ANALYZE HETEROGENEITY` requires a second pass of analysis to decide whether observed effect sizes (or test statistics) are homogeneous across samples.    

****How might you create the config file if your summary statistics files had different header labels?****   
The header for each of these columns must be specified so that METAL knows how to interpret the data. An example is given below:   

```
# === DESCRIBE AND PROCESS THE FOURTH INPUT FILE ===
MARKER MARKERNAME
ALLELE EFFECTALLELE NON_EFFECT_ALLELE
EFFECT EFFECT1
PVALUE PVALUE
WEIGHT NONMISS
PROCESS newfile4.txt
``` 

4. View the meta-analysis results

Some informative outptut was printed to stdout as METAL was running. ****What was the smallest p-value and how many markers was the meta-analysis completed for?****
There will be a .tbl and .tbl.info file created from the meta-analysis. You can use `less` to view the files.   
# Smallest p-value is 1.24e-652 at marker 'chr19:45412079:T:C'   

****Will we use the same genome wide significance threshold as in step 2? Why or why not?****  
We will use the same genome wide significance threshold as the number of independent markers across the genome has not changed.   

****How many genome wide significant results are there now?**** HINT: Use code like in #2 but replace `$10` with the column number with the p-value and use the file name for your meta-analysis results.

```
awk '$8 < 5e-8 {print 0}' LDL_METAL_MultiStudy.txt | wc -l   
#8293
```

5. Subset the results to markers in >1 study

METAL will perform a meta-analysis even on markers which are only present in one of the sub-studies. We are only interested in markers present in more than one study. 
The column labelled direction shows '?', '+', or '-' to indicate missingness, positive direction of effect, or negative direction of effect, respectively.  
One can use the `subset_meta_analysis.r` Rscript to exclude markers with two or more '?'.  
Execute the following command to subset the results. This will take < 5 minutes.  
`Rscript subset_meta_analysis.r --input LDL_METAL_META1.tbl --output LDL_METAL_MultiStudy.txt`

****How many genome wide significant results are there now?**** 

```
awk '$8 < 5e-8 {print 0}' LDL_METAL_MultiStudy.txt | wc -l   
#3292
```

6. Plot the meta-analysis results

To visually inspect your results for significant findings you can make a QQ-plot like Day 2's practical. 
****How does the inflation appear to you?****
Despiste not using `GENOIMIC CONTROL` the lambda is below 1.1, indicating that the results are not inflated.
   
`Rscript QQplot.r --input LDL_METAL_MultiStudy.txt --pvalue P-value --af Freq1 --prefix LDL_METAL_MultiStudy --break.top 120`  
Download the `*_QQ.png` file.  
From a terminal logged into `smed8020-2019-home`:   
`scp ubuntu@smed8020-2019-nodeXXX:day3/LDL_METAL_MultiStudy_QQ.png ~/`  
From a new terminal on your local machine:  
`scp smed8020-2019:LDL_METAL_MultiStudy_QQ.png .`  
The file should exist in whatever the default directory your terminal opens to. You can find this with `pwd`.  

****What is the lambda value for the smallest minor allele frequency (MAF) bin?****  
`cat *_lambda.txt`. 

```
lambda	frequency_bin
1.06591779149991	(0.05,0.5]
1.07819408540481	(0.005,0.05]
1.1080877366439	(0.001,0.005]
0.761649180425455	[0,0.001]
```

