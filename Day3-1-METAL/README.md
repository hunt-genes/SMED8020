# DAY 3 - Meta-analysis for quantitative traits using METAL

Please fill the following form during the exercise today: https://nettskjema.no/a/527809
Thank you!


## FOCUS and LEARNING GOALS

>  The aim for this session is to get familiar with running a genome wide association study meta-analysis which  enables researchers to  gather data from many studies and analyse them together. 

There are several motivations for meta-analysis. One is the ability to increase power to detect small effect sizes or rare variant effects by increasing the study sample size. Many methods for meta analysis rely on using summary statistics therefore rendering the need to share individual level data unnecessary. This makes it easier to share data for meta analysis as the summary statistics are not deemed sensitive information and are typically made publicly available when papers are published in peer-reviewed journals. Finally, meta-analysis across genetic ancestries is the most statistically robust approach rather than pooling all ancestries together in one GWAS as it results in little or no loss of efficiency (as compared to analysis of combined data-sets) and reduces population stratification.
Some [slides](SMED8020_2022_MetaAnalysis.pdf) with extra informaiton may be helpful.  

**Suggested reading:**

* [Willer, C. J., Li, Y. & Abecasis, G. R. METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics 26, 2190–2191 (2010).](https://academic.oup.com/bioinformatics/article/26/17/2190/198154)  
* [Nielsen, J. et al. Loss-of-function genomic variants with impact on liver-related blood traits highlight potential therapeutic targets for cardiovascular disease. Biorxiv. (2019) ](https://www.biorxiv.org/content/10.1101/597377v1)  
* [Evangelos Evangelou & John P. A. Ioannidis. Meta-analysis methods for genome-wide association studies and beyond. Nature Reviews Genetics. 2013](https://www.nature.com/articles/nrg3472)

[METAL](http://csg.sph.umich.edu/abecasis/metal/) was developed at the University of Michigan as a tool for meta-analysis of  genome-wide association analysis. For running METAL you can use either test statistics and standard error or p-values. For more info on METAL see the web-links in this document and the suggested readings paper. 

## TASK: Running a multi-ancestry meta-analysis using METAL  

Today you will run a meta-analysis to combine three studies using METAL. METAL has been pre-installed on our lab. Because of time restraints we have made a small data-set that will run within reasonable time. The data will not generate significant results. A separate set of files will therefore be used for plotting results. 
    
The phenotype in today's practical is low density lipoprotein (LDL) cholesterol and we will be using data from three large studies: HUNT, Biobank Japan and Global Lipids Genetics Consortium (GLGC).  

### Task outline  

0. Install METAL
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

This practical makes use of many Unix/Linux commands in the terminal like `head`, `awk`, `cut`. You can read about these [here].

To know more about the steps we have done for you (e.g. liftOver) [read this!](https://github.com/hunt-genes/SMED8020/blob/main/Day3-1-METAL/ForReference.md)

## Instructions  

#### 0. Install METAL
   
For today's practical we have alreading installed the correct version here:    
`/mnt/scratch/software/METAL-2020-05-05/build/bin/metal`    

#### 1. Gather summary statistics  

Usually you would download publically available summary statistics from the internet to your local machine. For convience for this practical, we have already downloaded summary statistics from 3 studies, BBJ, HUNT, and GLGC.

# Required data

You already copied the required data (See the main [README](https://github.com/hunt-genes/SMED8020/blob/main/README.md) for reference)

Navigate to the directory with the data, remember to replace $USER with your username.
```
cd /mnt/work/workbench/$USER/SMED8020/Day3-1-METAL/
```

* The original summary statistics from Biobank Japan (BBJ) of LDL cholesterol in N=72,866 can be found [here](https://humandbs.biosciencedbc.jp/files/hum0014/hum0014_README_QTL_GWAS.html)  
```
head BBJ-LDL-preMeta.txt
```
The columns are CHR     POS38   SNPID   Allele1 Allele2 AC_Allele2      AF_Allele2      N       BETA    SE      p.value log10P

* The original summary statistics of joint analysis of metabochip and GWAS data for LDL cholesterol in N=89,138 from the Global Lipids Genetics Consortium (GLGC) can be found [here](http://csg.sph.umich.edu/willer/public/lipids2013/)  
```
head GLGC-LDL-preMeta.txt
``` 
The columns are SNP_hg18        SNP_hg19        rsid    A1      A2      beta    se      N       P-value Freq.A1.1000G.EUR

* The summary statistics of LDL cholesterol from the HUNT study in N=67,429.   
```
head HUNT-LDL-preMeta.txt
``` 
The columns are CHR     POS38   SNPID   Allele1 Allele2 AC_Allele2      AF_Allele2      N       BETA    SE  p.value 

#### 2. Check your summary statistics to make sure they're ready for meta-analysis.

### 2.1 Which genome build is used?

The human reference genome has been updated over the years and variants are given different coordinates in different versions. 
The latest human reference genome GRCh38 was released from the Genome Reference Consortium on 17 December 2013.  
The previous human reference genome (GRCh37) was the nineteenth version (hg19).  
The version before this was NCBI Build 36.1	released March 2006	(hg18). 
You can see more [here](https://genome.ucsc.edu/FAQ/FAQreleases.html#release1). hg19 is still widely used and people are slowly converting to hg38.       

****From the summary statistic headers, can you tell what reference genome versions are used for each study?****  

It looks like BBJ and HUNT have SNP coordinates from hg38, but GLGC has summary statistics from hg18 and hg19. 
We can use [UCSC listOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert the hg19 coordinates to hg38 before meta-analysis. We can use liftOver on the command line or via the web. To avoid extensive file manipulation on your part, we already used this .bed file to make a new version of the GLGC results: `GLGC-LDL-hg38-preMeta.txt`. This file also has a header that is consistent with the other two files. 

```
head GLGC-LDL-hg38-preMeta.txt
```

Look in `GLGC.hg38.unmapped`. ****Were there some markers that did not get converted from hg19 to hg38? Why do you think that is?****       

You would then use R or another tool to merge GLGC.hg38.bed and GLGC-LDL-preMeta.txt (the hg19 position should be shared between them).  To simplify the practical we have done this for all files:    

* BBJ-LDL-preMeta-U.txt
* HUNT-LDL-preMeta-U.txt
* GLGC-LDL-hg38-preMeta-U.txt
  

### 2.2 Check the file formats and headers   

****What is the header of each file? What does the `-n 1` parameter do in `head`?****
In terminal, use the `head` command with the flag `-n 1` to check the headers of the files:    
```
#Use head to check the headers
head -n 1 BBJ-LDL-preMeta-U.txt
head -n 1 HUNT-LDL-preMeta-U.txt
head -n 1 GLGC-LDL-hg38-preMeta-U.txt
```

****Are your SNPIDs across the files formatted in the same way?****     
```
#Use head to check SNPID formatting
head -n 2 BBJ-LDL-preMeta-U.txt | cut -f 3 
head -n 2 GLGC-LDL-hg38-preMeta-U.txt | cut -f 3
head -n 2 HUNT-LDL-preMeta-U.txt | cut -f 3 
```
Yes, we need the SNPID to be consistent across files. The header could be called something different, but the software will match the markers across studies based on the column. 

### 2.3 How many variants will we be meta-analyzing?     
****How many variants are in each of the files?****  

In terminal:
```
#use the wc function to count the line numbers in the files
wc -l BBJ-LDL-preMeta-U.txt
wc -l HUNT-LDL-preMeta-U.txt
wc -l GLGC-LDL-hg38-preMeta-U.txt
```

```
#check the original files before formatting and subsetting 
wc -l BBJ-LDL-preMeta.txt
wc -l HUNT-LDL-preMeta.txt
wc -l GLGC-LDL-hg38-preMeta.txt
```

The HUNT and BBJ summary statistics originally had millions of variants because imputation was done with the TOPMed imputation panel, which allows for higher resolution imputation due to the large amount of sequencing samples which make up the reference panel. We have subsetted the input files to only include variants seen in all 3 studies. We only want to perform meta-analysis on variants tested in 2 or more studies.    

****What imputation panel was used for GLGC?**** HINT: Check the methods of the [paper](https://www.nature.com/articles/ng.2797). This is sort of a trick question, there wasn't a clear imputation method described in the 2013 paper. See the answers for an explanation.   

****How many genome wide significant results are in each of the input files?****    

In terminal:
```
#use awk to identify the rows that have a p-value < 5E-8
awk '$11 < 5e-8 {print 0}' HUNT-LDL-preMeta-U.txt | wc -l
awk '$11 < 5e-8 {print 0}' BBJ-LDL-preMeta-U.txt | wc -l
awk '$10 < 5e-8 {print 0}' GLGC-LDL-hg38-preMeta-U.txt | wc -l
```   

#### 3. Running METAL   

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
 
### 3.1. Create config file
A shell wrapper script will be used to create the config file needed to run METAL. This script, `LDL_metal.sh`, has been created for you. Create a config file with the bash script `LDL_METAL.sh` by filling in the appropriate arguments instead of "file1",  "file2",  "file3" and using "LDL_METAL" as your output prefix.

In terminal:
```
bash LDL_metal.sh HUNT-LDL-preMeta-U.txt GLGC-LDL-hg38-preMeta-U.txt BBJ-LDL-preMeta-U.txt LDL_METAL_META > LDL_METAL.conf
```   

In terminal:
```
# You can view the config file with less
less LDL_METAL.conf

#press "q" to exit
```   

### 3.2. Run metal 
We will run METAL with the config file we just made. This should take less than 20 minutes.

Remember, we have already installed METAL here:    
`/mnt/scratch/software/METAL-2020-05-05/build/bin/metal`   

In terminal:
```
#Run METAL using the config file
/mnt/scratch/software/METAL-2020-05-05/build/bin/metal LDL_METAL.conf > LDL_METAL.log
```     

Note: If you would like to time your analysis you can use the time program.  
`/usr/bin/time -o test_time -v /home/benb/scratch/software/METAL-2020-05-05/build/bin/metal LDL_METAL.conf`    

While the meta-analysis runs, consider the following questions:  

****What type of meta-analysis did you run (fixed or random effects? sample size or inverse variance based?) What is the difference?****     

****Did you use genomic control? In what situations is it useful to use genomic control****       

****What does it mean to set the minimum weight to 10,000?****        

****What is the difference between "ANALYZE" and "ANALYZE HETEROGENEITY"?****       

****How might you create the config file if your summary statistics files had different header labels?****       
Example of a 4th study that you might want to include (newfile4.txt):    
The columns are:    
CHR     POS38   MARKERNAME   NON_EFFECT_ALLELE EFFECTALLELE      AF_EFFECTALLELE      NONMISS       EFFECT    SE  PVALUE     

And here is the first line of the file:    
18      32811534        18:32811534:C:G C       G       0.45      87041   -0.0015 0.0075  0.7029    


#### 4. View the meta-analysis results

Some informative output was printing to "standard output" as METAL was running. We saved it in a file named `LDL_METAL.log`. Check out the information there. Just so you know, "standard output" is called stdout, and in the terminal, stdout defaults to the user's screen.    

****What was the smallest p-value and how many markers was the meta-analysis completed for?****     

There will be a .tbl and .tbl.info file created from the meta-analysis. You can use `less` to view the files.

In terminal:
```
#output file
less LDL_METAL_META1.tbl
```   
```
#info about the output file format
less LDL_METAL_META1.tbl.info
```   
```
#where we saved the stdout when we ran METAL
less LDL_METAL.log
```   

****Do you think we will we use the same genome-wide significance threshold (5xE-8) for the meta-analysis as we used for the GWAS? Why or why not?****  

****How many genome wide significant results are there now?****    
HINT: Use code like in *2.3* but replace `$8` with the column number that has the p-value and use the file name for your meta-analysis results.

## 5. Subset to markers in more than 1 study
Note: We pre-processed the files so you don't have to subset the results to markers in >1 study, but you might need this information in the future if you have not pre-processed your input files.

METAL will perform a meta-analysis even on markers which are only present in one of the sub-studies. We are only interested in markers present in more than one study. 
The column labelled "direction" shows '?', '+', or '-' to indicate missingness, positive direction of effect, or negative direction of effect, respectively.  
One can use the `subset_meta_analysis.r` Rscript to exclude markers with more than one '?'. This is R code like we use in RStudio, but it's packaged in a script so we can call it from the command line and pass it parameters, like the input file.

First install the R packages we need:
In terminal:
```
#Type R  and press enter to open R
install.packages("data.table")
install.packages('stringr')
install.packages("optparse")
#Type q() and press "n" to exit
```

```
#subset the results to variants with more than 1 study, may take 5 minutes
Rscript subset_meta_analysis.r --input LDL_METAL_META1.tbl --output LDL_METAL_MultiStudy.txt
#if this doesn't work in the terminal, open subset_metal_analysis_manual.r. Add in the file names and parameters to run it.
```

#### 6. Plot the meta-analysis results

#### 6.1 QQ-plot
To visually inspect your results for significant findings you can make a QQ-plot. We have a script `QQplot.R` which creates an image file with the plot and a text file with lambda values. You don't need RStudio for this.

First install the R packages we need:
In terminal:
```
#Type R  and press enter to open R
install.packages("plotrix")
install.packages("data.table")
install.packages("RColorBrewer")
install.packages("optparse")
#Type q() and press "n" to exit
```

In terminal:
```
#make a QQplot using an Rscript
Rscript QQplot.r --input LDL_METAL_META1.tbl --pvalue P-value --af Freq1 --prefix LDL_METAL_MultiStudy --break.top 120
#if this doesn't work, open QQplot_manual.r in RStudio. Add in the file names and parameters to run it.
```
The image file should exist in whatever the default directory your R is writing into, which should be your current working directory. You can find this with `pwd`. Open the file to inspect the QQ-plot.

****How does the inflation appear to you?****  

****What is the lambda value for the smallest minor allele frequency (MAF) bin?****  
In terminal:
```
#use cat to view the file of lambda values
cat *_lambda.txt
```

#### 6.2 Forest plot
Another useful comparison of input studies and the meta-analysis is a Forest plot. 
You can read more about R code to make this plot [here](https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html).

We will make a Forest plot for a lead SNP in APOE.

In RStudio:
Open the ForestPlot.R script. 
Run the code line by line to generate a forest plot. 

#### 6.3 Manhattan plot (Currently broken, sorry!)

Copy the script into your folder for meta-analysis.
```
cd /mnt/work/workbench/$USER/SMED8020/Day3-1-METAL/
mv /mnt/scratch/data/ManhattanPlot.r .
```
You may need to specify some of the parameters. Use `Rscript ManhattanPlot.r --help` to figure out what.
In terminal:
```
Rscript ManhattanPlot.r --input LDL_METAL_META1.tbl --pvalue P-value --af Freq1 --prefix LDL_METAL_MultiStudy
```

#### 6.4 Open Targets
Check out the APOE region on the [Open Targets platform](https://genetics.opentargets.org/variant/19_44886339_G_A). This platform integrates a lot of data for interrogating genetic variants as drug targets.

I would also recommend [this example](https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html) of a meta-analysis by Matti Pirinen at the University of Helsinki using R.

You can review [answers](Day3_MetaAnalysis_ANSWERS.md) to the questions in this practical.
