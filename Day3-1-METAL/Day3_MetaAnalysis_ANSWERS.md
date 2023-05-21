#### 2. Check your summary statistics to make sure they're ready for meta-analysis.

****From the summary statistic headers, can you tell what reference genome versions are used for each study?****  

It looks like BBJ and HUNT have SNP coordinates from hg38, but GLGC has summary statistics from hg18 and hg19.    

****Were there some markers that did not get converted from hg19 to hg38? Why do you think that is?****       

Use less to view GLGC.hg38.unmapped:
In terminal:    
```less GLGC.hg38.unmapped```.   

All of the regions were not converted because they were "Deleted in new". This means that there was no region in hg19 that aligned well to your region in hg18. 

Other error messages that can arise:
-"Partially deleted in new" means that there was a piece of an alignment that matched your hg18 region, but not enough to convert the region to hg19. 
-"Split in new", which means that your hg18 region was split up into different parts of hg19. 

### 2.2 Check the file formats and headers 

****What is the header of each file? What does the `-n 1` parameter do in `head`?****
In terminal:    
```
# This prints the first n lines of the file 
head -n 1 file
```

```
#Use head to check the headers
head -n 1 BBJ-LDL-preMeta-U.txt
head -n 1 HUNT-LDL-preMeta.txt
head -n 1 GLGC-LDL-hg38-preMeta-U.txt
```
The columns are:    
CHR: Chromosome    
POS38: Position in hg38    
SNPID: Marker ID    
Allele1: Allele1    
Allele2: Allele2    
AC_Allele2: Allele Count    
AF_Allele2: Allele 2 Ferquency    
N: Sample size    
BETA: Effect size    
SE: Standard error    
p.value: P-value    

****Are your SNPIDs across the files formatted in the same way?****     
```
#Use head to check SNPID formatting
head -n 2 BBJ-LDL-preMeta-U.txt | cut -f 3 
head -n 2 GLGC-LDL-hg38-preMeta-U.txt | cut -f 3
head -n 2 HUNT-LDL-preMeta-U.txt | cut -f 3 
```
Yes, we need the SNPID to be consistent across files. The header could be called something different, but the software will match the markers across studies based on the column. 

Yes. You can check this using the grep command by checking at least one SNPID in each set of results
`grep 'chr18:4440419:G:T' BBJ-LDL-preMeta-U.txt`
chr18:4440419:G:T #BBJ-LDL-preMeta-U.txt
chr18:4440419:G:T #GLGC-LDL-hg38-preMeta.txt
chr18:4440419:G:T #HUNT-LDL-preMeta.txt

2.3 How many variants will we be meta-analyzing?     
****How many variants are in each of the files?****       
```
#use the wc function to count the original line numbers in the files
wc -l BBJ-LDL-preMeta.txt
wc -l HUNT-LDL-preMeta.txt
wc -l GLGC-LDL-hg38-preMeta.txt
```
6108162 BBJ-LDL-preMeta.txt
5788779 HUNT-LDL-preMeta.txt
2436641 GLGC-LDL-hg38-preMeta.txt      

The HUNT summary statistics originally had millions of variants because imputation was done with the TOPMed imputation panel, which allows for higher resolution imputation due to the large amount of sequencing samples which make up the reference panel. We have subsetted the input files to only include variants seen in all 3 studies. We only want to perform meta-analysis on variants tested in 2 or more studies.    

```
#use the wc function to count the formatted and subset line numbers in the files
wc -l BBJ-LDL-preMeta-U.txt
wc -l HUNT-LDL-preMeta-U.txt
wc -l GLGC-LDL-hg38-preMeta-U.txt
```
2007882 BBJ-LDL-preMeta-U.txt   
2007882 HUNT-LDL-preMeta-U.txt   
2007882 GLGC-LDL-hg38-preMeta-U.txt   

****What imputation panel was used for GLGC?**** HINT: Check the methods of the [paper](https://www.nature.com/articles/ng.2797).    

GWAS data were imputed to the International HapMap project (https://www.nature.com/articles/nature09270?page=55#Sec2).

****How many genome wide significant results are in each of the input files?****    

In terminal:
```
#use awk to identify the rows that have a p-value < 5E-8
awk '$11 < 5e-8 {print 0}' HUNT-LDL-preMeta-U.txt | wc -l
awk '$11 < 5e-8 {print 0}' BBJ-LDL-preMeta-U.txt | wc -l
awk '$10 < 5e-8 {print 0}' GLGC-LDL-hg38-preMeta-U.txt | wc -l
```   
1350
957
2397

#### 3. Running METAL   

****What type of meta-analysis did you run (fixed or random effects? sample size or inverse variance based?) What is the difference?****     

Fixed effects using inverse variance based weights.   
    
We ran the `SCHEME STDERR` command which weights the effect size estimates using the inverse of the corresponding standard errors. To enable this option, we specifed which of the input columns contained standard error information using the `STDERRLABEL` command (or `STDERR` for short). While standard error based weights are more common in the biostatistical literature, if you decide to use this approach, it is very important to ensure that effect size estimates (beta coefficients) and standard errors use the same units in all studies (i.e. make sure that the exact same trait was examined in each study and that the same transformations were applied). Inconsistent use of measurement units across studies is the most common cause of discrepancies between these two analysis strategies.

The default - METAL combines p-values across studies taking into account a study specific weight (typically, the sample size) and direction of effect. This behavior can be requested explicitly with the `SCHEME SAMPLESIZE` command. 

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

#### 4. View the meta-analysis results

****What was the smallest p-value and how many markers was the meta-analysis completed for?****     

There will be a .tbl and .tbl.info file created from the meta-analysis. You can use `less` to view the files.

In terminal:
  
```
#where we saved the stdout when we ran METAL
less LDL_METAL.log
```   
Completed meta-analysis for 2007881 markers!     
Smallest p-value is 6.40e-1688 at marker '19:44908822:C:T'    

****Do you think we will we use the same genome-wide significance threshold (5xE-8) for the meta-analysis as we used for the GWAS? Why or why not?****  
We will use the same genome wide significance threshold as the number of independent markers across the genome has not changed.   

****How many genome wide significant results are there now?****    
HINT: Use code like in *2.3* but replace `$8` with the column number that has the p-value and use the file name for your meta-analysis results.


```
awk '$8 < 5e-8 {print 0}' LDL_METAL_META1.tbl | wc -l   
#4066
```

****How does the inflation appear to you?****  

The plot looks slightly inflated. However LDL is a highly polygenic trait and therefore it could be due to polygenicity as opposed to inflation. Testing for inflation with LD score regression to help confirm this.    

****What is the lambda value for the smallest minor allele frequency (MAF) bin?****  
In terminal:
```
#use cat to view the file of lambda values
cat *_lambda.txt
```

```
lambda  frequency_bin
1.16451105532861        (0.05,0.5]
1.07376256781323        (0.005,0.05]
1.33349430385881        [0.0017,0.005]
```
