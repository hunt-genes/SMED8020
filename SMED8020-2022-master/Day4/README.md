# DAY 4 - Practical Exercise – Polygenic Risk Scores
## Overview
In this practical, we will generate at polygenic risk score (PRS) for height and see how much variation in height it explains. This practical is a reduced set of analyses provided in the tutorial by Shing Wan Choi and colleagues, which can be found at the following link[PRS Guide](https://choishingwan.github.io/PRS-Tutorial/). Please review the tutorial for more detailed information, and calculation of the PRS using the additional methods that were discussed in the lecture.

The practical is separated into four main sections, corresponding to the guide in the Nature Protocols paper discussed in the lecture [PRS Guide paper](https://www.nature.com/articles/s41596-020-0353-1):
1.	Quality control (QC) of the base data
2.	QC of the target data
3.	Calculating and analysing the PRS
4.	Visualising the PRS results
Please refer to the paper as you work through the practical. 

Any command line code will be presented in red text, R code in green text and questions to answer in blue text.

# QC of Base Data
The first step in Polygenic Risk Score (PRS) analyses is to generate or obtain the base data (GWAS summary statistics). Ideally these will correspond to the most powerful GWAS results available on the phenotype under study. In this example, we will use GWAS summary statistics from simulated height data (these will be provided in the folder with this document).

## Reading the base data file
**Height.gwas.txt.gz** is compressed. To read its content, you can type:

```bash
gunzip -c Height.gwas.txt.gz | head
```

which will display the first 10 lines of the file. The **Height.gwas.txt.gz** file contains the following columns:

- **CHR**: The chromosome in which the SNP resides
- **BP**: Chromosomal co-ordinate of the SNP
- **SNP**: SNP ID, usually in the form of rs-ID
- **A1**: The effect allele of the SNP
- **A2**: The non-effect allele of the SNP
- **N**: Number of samples used to obtain the effect size estimate
- **SE**: The standard error (SE) of the effect size esimate
- **P**: The P-value of association between the SNP genotypes and the base phenotype
- **OR**: The effect size estimate of the SNP, if the outcome is binary/case-control. If the outcome is continuous or treated as continuous then this will usually be BETA
- **INFO**: The imputation information score
- **MAF**: The minor allele frequency (MAF) of the SNP

## QC checklist: Base data

## \# Heritability check
We recommend that PRS analyses are performed on base data with a chip-heritability estimate $h_{snp}^{2} > 0.05$. 
The chip-heritability of a GWAS can be estimated using e.g. LD Score Regression (LDSC). 
Our height GWAS data are simulated to have a chip-heritability much greater than 0.05 and so we can move on to the next QC step. 

## \# Effect allele
It is important to know which allele is the effect allele and which is the non-effect allele for PRS association results to be in the correct direction.

!!! Important
    Some GWAS results files do not make clear which allele is the effect allele and which is the non-effect allele.
    If the incorrect assumption is made in computing the PRS, then the effect of the PRS in the target data will be in the wrong direction.

    To avoid misleading conclusions the effect allele from the base (GWAS) data must be known.


## \# Genome build
The height summary statistic are on the same genome build as the target data that we will be using. 
You must check that your base and target data are on the same genome build, and if they are not then use a tool such as [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to make the builds consistent across the data sets.

## \# Standard GWAS QC
As described in the paper, both the base and target data should be subjected to the standard stringent QC steps performed in GWAS. 
If the base data have been obtained as summary statistics from a public source, then the typical QC steps that you will be able to perform on them are to filter the SNPs according to INFO score and MAF. 
SNPs with low minor allele frequency (MAF) or imputation information score (INFO) are more likely to generate false positive results due to their lower statistical power (and higher probability of genotyping errors in the case of low MAF). 
Therefore, SNPs with low MAF and INFO are typically removed before performing downstream analyses.
We recommend removing SNPs with MAF < 1% and INFO < 0.8 (with very large base sample sizes these thresholds could be reduced if sensitivity checks indicate reliable results).
These SNP filters can be achieved using the following code:

=== "Using bash"
    ```bash 
    gunzip -c Height.gwas.txt.gz |\
    awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
    gzip  > Height.gz
    ```
The bash code above does the following:

1. Decompresses and reads the **Height.gwas.txt.gz** file
2. Prints the header line (`NR==1`)
3. Prints any line with MAF above 0.01 (`$11` because the eleventh column of the file contains the MAF information)
4. Prints any line with INFO above 0.8 (`$10` because the tenth column of the file contains the INFO information)
5. Compresses and writes the results to **Height.gz**

## \# Mismatching SNPs
SNPs that have mismatching alleles reported in the base and target data are either resolvable by "strand-flipping" the alleles to their complementary alleles in e.g. the target data, such as for a SNP with A/C in the base data and G/T in the target, or non-resolvable, such as for a SNP with C/G in the base and C/T in the target. 
Most polygenic score software perform strand-flipping automatically for SNPs that are resolvable, and remove non-resolvable mismatching SNPs.


Since we need the target data to know which SNPs have mismatching alleles, we will perform this strand-flipping in the target data.

## \# Duplicate SNPs
If an error has occurred in the generation of the base data then there may be duplicated SNPs in the base data file.
Most PRS software do not allow duplicated SNPs in the base data input and thus they should be removed, using a command such as the one below: 

```bash
gunzip -c Height.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip - > Height.nodup.gz
```

The above command does the following:

1. Decompresses and reads the **Height.gz** file
2. Count number of time SNP ID was observed, assuming the third column contian the SNP ID (`seen[$3]++`). If this is the first time seeing this SNP ID, print it. 3. 
Compresses and writes the results to **Height.nodup.gz**

??? note "How many duplicated SNPs are there?"
    There are a total of `2` duplicated SNPs

## \# Ambiguous SNPs
If the base and target data were generated using different genotyping chips and the chromosome strand (+/-) that was used for either is unknown, then it is not possible to pair-up the alleles of ambiguous SNPs (i.e. those with complementary alleles, either C/G or A/T SNPs) across the data sets, because it will be unknown whether the base and target data are referring to the same allele or not. While allele frequencies could be used to infer which alleles are on the same strand, the accuracy of this could be low for SNPs with MAF close to 50% or when the base and target data are from different populations. Therefore, we recommend removing all ambiguous SNPs to avoid introducing this potential source of systematic error.

Ambiguous SNPs can be removed in the base data and then there will be no such SNPs in the subsequent analyses, since analyses are performed only on SNPs that overlap between the base and target data.

Nonambiguous SNPs can be retained using the following:
```bash
gunzip -c Height.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > Height.QC.gz
```

??? note "How many non-ambiguous SNPs were there?"
    There are `499,617` non-ambiguous SNPs



## \# Sample overlap and relatedness
Since the target data were simulated there are no overlapping samples between the base and target data here (see the relevant section of the paper for discussion of the importance of avoiding sample overlap).
Closely related individuals within and between the base and the target data may lead to overfitted results, limiting the generalizability of the results (see the relevant sections of the paper). Relatedness within the target data is tested in the Target Data section.

The Height.QC.gz base data are now ready for using in downstream analyses.

# QC of Target Data
Target data consist of individual-level genotype-phenotype data, usually generated within your lab/department/collaboration. For this tutorial, we have simulated some genotype-phenotype data using the 1000 Genomes Project European samples. You can download the data here

Unzip the data as follow:

```bash
unzip EUR.zip
```
## QC checklist: Target data
Below are the QC steps that comprise the QC checklist for the target data.

## \# Sample size
We recommend that users only perform PRS analyses on target data of at least 100 individuals. The sample size of our target data here is 503 individuals. 

## \# Genome build
As stated in the base data section, the genome build for our base and target data is the same, as it should be.

## \# Standard GWAS QC
The target data must be quality controlled to at least the standards 
implemented in GWAS studies, e.g. removing SNPs with low genotyping rate, 
low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing
individuals with low genotyping rate 
(see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).

The following `plink` command applies some of these QC metrics to the target data:


```bash
plink \
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC
```

Each of the parameters corresponds to the following

| Parameter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| maf | 0.01 | Removes all SNPs with minor allele frequency less than 0.01. Genotyping errors typically have a larger influence on SNPs with low MAF. Studies with large sample sizes could apply a lower MAF threshold|
| hwe | 1e-6 | Removes SNPs with low P-value from the Hardy-Weinberg Equilibrium Fisher's exact or chi-squared test. SNPs with significant P-values from the HWE test are more likely affected by genotyping error or the effects of natural selection. Filtering should be performed on the control samples to avoid filtering SNPs that are causal (under selection in cases). When phenotype information is included, plink will automatically perform the filtering in the controls. |
| geno | 0.01 | Excludes SNPs that are missing in a high fraction of subjects. A two-stage filtering process is usually performed (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).|
| mind | 0.01 | Excludes individuals who have a high rate of genotype missingness, since this may indicate problems in the DNA sample or processing. (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/) for more details).|
| make-just-fam | - | Informs `plink` to only generate the QC'ed sample name to avoid generating the .bed file.  |
| write-snplist | - | Informs `plink` to only generate the QC'ed SNP list to avoid generating the .bed file. |
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |

??? note "How many SNPs and samples were filtered?"
    - `14` samples were removed due to a high rate of genotype missingness
    - `5,353` SNP were removed due to missing genotype data
    -  `944` SNPs were removed due to being out of Hardy-Weinberg Equilibrium
    - `5,061` SNPs were removed due to low minor allele frequency

!!! note
    Normally, we can generate a new genotype file using the new sample list.
    However,  this will use up a lot of storage space. Using `plink`'s
    `--extract`, `--exclude`, `--keep`, `--remove`, `--make-just-fam` and `--write-snplist` functions, we can work 
    solely on the list of samples and SNPs without duplicating the 
    genotype file, reducing the storage space usage.  
    
Very high or low heterozygosity rates in individuals could be due to DNA contamination or to high levels of inbreeding. Therefore, samples with extreme heterozygosity are typically removed prior to downstream analyses. 

First, we perform pruning to remove highly correlated SNPs:

```bash
plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC
```

Each of the parameters corresponds to the following

| Parameter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| keep | EUR.QC.fam | Informs `plink` that we only want to use samples in `EUR.QC.fam` in the analysis |
| extract | EUR.QC.snplist | Informs `plink` that we only want to use SNPs in `EUR.QC.snplist` in the analysis |
|indep-pairwise| 200 50 0.25 | Informs `plink` that we wish to perform pruning with a window size of 200 variants, sliding across the genome with step size of 50 variants at a time, and filter out any SNPs with LD $r^2$ higher than 0.25|
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |


This will generate two files 1) **EUR.QC.prune.in** and 2) **EUR.QC.prune.out**. All SNPs within **EUR.QC.prune.in** have a pairwise $r^2 < 0.25$. 


Heterozygosity rates can then be computed using `plink`:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC
```

This will generate the **EUR.QC.het** file, which contains F coefficient estimates for assessing heterozygosity.
We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean, which can be performed using the following `R` command (assuming that you have R downloaded, then you can open an `R` session by typing `R` in your terminal):

```R
    dat <- read.table("EUR.QC.het", header=T) # Read in the EUR.het file, specify it has header
    m <- mean(dat$F) # Calculate the mean  
    s <- sd(dat$F) # Calculate the SD
    valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
    write.table(valid[,c(1,2)], "EUR.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
    q() # exit R
    ```
??? note "How many samples were excluded due to high heterozygosity rate?"
    - `2` samples were excluded

## \# Ambiguous SNPs
These were removed during the base data QC.

## \# Mismatching SNPs
SNPs that have mismatching alleles reported in the base and target data may be resolvable by strand-flipping the alleles to their complementary alleles in e.g. the target data, such as for a SNP with A/C in the base data and G/T in the target. Most PRS software will perform strand-flipping automatically. Check your software does this before calculating your PRS; if it does not, there is some R code in the online tutorial for performing strand-flipping (https://choishingwan.github.io/PRS-Tutorial/target/).

```bash
Rscript Ambiguous_SNPs.R
```

## \# Duplicate SNPs
Make sure to remove any duplicate SNPs in your target data (these target data were simulated and so include no duplicated SNPs).

## \# Sex chromosomes 
Sometimes sample mislabelling can occur, which may lead to invalid results. One indication of a mislabelled sample is a difference between reported sex and that indicated by the sex chromosomes. While this may be due to a difference in sex and gender identity, it could also reflect mislabeling of samples or misreporting and, thus, individuals in which there is a mismatch between biological and reported sex are typically removed. A sex check can be performed in PLINK, in which individuals are called as females if their X chromosome homozygosity estimate (F statistic) is < 0.2 and as males if the estimate is > 0.8.

Before performing a sex check, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
A sex check can then easily be conducted using `plink`
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC
```

This will generate a file called **EUR.QC.sexcheck** containing the F-statistics for each individual. Individuals are typically called as being biologically male if the F-statistic is > 0.8 and biologically female if F < 0.2.

```R
    # Read in file
    valid <- read.table("EUR.valid.sample", header=T)
    dat <- read.table("EUR.QC.sexcheck", header=T)
    valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
    write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
    q() # exit R
    ```
??? note "How many samples were excluded due mismatched Sex information?"
    - `4` samples were excluded

## \# Sample overlap
Since the target data were simulated there are no overlapping samples between the base and target data here (see the relevant section of [the paper](https://www.nature.com/articles/s41596-020-0353-1) for discussion of the importance of avoiding sample overlap). 

## \# Relatedness
Closely related individuals in the target data may lead to overfitted results, limiting the generalisability of the results. 

Before calculating the relatedness, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
Individuals that have a first or second degree relative in the sample ($\hat{\pi} > 0.125$) can be removed with the following command:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC
```

??? note "How many related samples were excluded?"
    - `0` samples were excluded

## Generate final QC'ed target data file
After performing the full analysis, you can generate a QC'ed data set with the following command:

```bash
plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch \
    --a1-allele EUR.a1
```

Each of the parameters corresponds to the following

| Parameter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| keep | EUR.QC.rel.id | Informs `plink` that we only want to keep samples in `EUR.QC.rel.id` |
| extract | EUR.QC.snplist | Informs `plink` that we only want to use SNPs in `EUR.QC.snplist` in the analysis |
| exclude | EUR.mismatch | Informs `plink` that we wish to remove any SNPs in `EUR.mismatch`|
| a1-allele |  EUR.a1 | Fix all A1 alleles to those specified in `EUR.a1` |
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |