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
2. Count number of time SNP ID was observed, assuming the third column contian the SNP ID (`seen[$3]++`). If this it the first time seeing this SNP ID, print it. 3. Compresses and writes the results to **Height.nodup.gz**

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

