
## Set up 
#set working directory, here exemplified with a folder on the Desktop of the user tno010
# Example: setwd("C://Users/tno010/Desktop/SMED8020/2021/Day2/") or in workbench e.g.
setwd("/mnt/work/..../SMED8020/Day2/")

## Then investigate the results from the logistic regressions before QC is done
# Reading in the results from the logistic regressions performed in Plink
assoc <- read.table("output/results.assoc.logistic", header=T, as.is=T)

# You need to have the packages "ggplot2" and "scales" installed, do that by:
# if you have them, no need to run line 10
#install.packages("ggplot2", "scales")

# load the R packages (making them active)
library("ggplot2", "scales")

# make a function to reverse and log transform the p-values
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf))
}

# make the plot
assoc_plot <- ggplot(assoc, aes(BP, P)) +
  geom_point() +  scale_y_continuous(trans=reverselog_trans(base=10),
                                     labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="-log10 p-values", x ="position", y = "-log10 p-value") +
  theme_bw()

# look at the plot
assoc_plot

#save the plot
ggsave("output/assoc_plot.pdf", assoc_plot)

## Perform and evaluate Quality Control (QC) alongside Plink

### Sample Filters
# (a)  elevated missing data rates or outlying heterozygosity rate**
### Evaluating QC results - missingness and heterozygosity - from Plink in R 

# Reading the files:
het <- read.table("output/het.het", header=T)
mis <- read.table("output/mis.imiss", header=T)

#First, calculate the observed heterozygosity rate per individual using the formula (N(NM) - O(Hom))/N(NM)
#Then, plot the proportion of missing genotypes and the heterozygosity rate
# Create new data frame with only relevant info for plotting
mishet <- data.frame(FID=het$FID, IID=het$IID, het.rate=(het$N.NM - het$O.HOM)/het$N.NM, mis.rate=mis$F_MISS)

# make plot of missing genotypes and the heterozygosity rate
mishet_plot <- ggplot(mishet, aes(mis.rate, het.rate)) +
  geom_point() +
  labs(title="", x ="Proportion of missing genotypes", y = "Heterozygosity rate") +
  theme_bw()

# investigate plot
mishet_plot

# save plot
ggsave("output/mishet.pdf", mishet_plot)


## Evaluating which individuals to exclude based on set cut off values
# by making plot
mishet_excl_plot <- ggplot(mishet, aes(mis.rate, het.rate)) +
  geom_point() +
  labs(title="", x ="Proportion of missing genotypes", y = "Heterozygosity rate") +
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dashed") +
  geom_hline(yintercept=0.23, linetype="dashed") +
  geom_hline(yintercept=0.4, linetype="dashed")

# investigate plot
mishet_excl_plot

# save plot
ggsave("output/mishet_excl_1.pdf", mishet_excl_plot)

#How many samples failed each of the filters? And what is the total number of failures?
# individuals with mis.rate > 0.01
fail_mis_qc <- mishet[mishet$mis.rate > 0.01,]

# individuals with het.rate <0.23 and individuals with  het.rate >0.4
fail_het_qc <- mishet[mishet$het.rate < 0.23 | mishet$het.rate > 0.4,]

# save details for individuals that fail these criteria:
write.table(fail_mis_qc, "output/fail_mis_qc.txt", row.names=F, col.names=T, quote=F)
write.table(fail_het_qc, "output/fail_het_qc.txt", row.names=F, col.names=T, quote=F)


## Identify duplicated or related individuals using the identity-by-descent (IBD) report from Plink
# The column [PI_HAT] is a mean IBD per individual, i.e. (0\*genome$Z0 + 1\*genome$Z1 + 2\*genome$Z2)/2
#PI_HAT > 0.1875 corresponds to a half way between second and third degree relatives.

#read in day2pruned.genome file
genome <- read.table("output/pruned.genome", 
                     header=T, as.is=T)
genome <- genome[genome$PI_HAT > 0.1875,]

# Identify duplicated individuals by plotting the standard error of the IBD sharing vs. the mean IBD sharing:
# compute Mean(IBD)
mean.ibd <- 0*genome$Z0 + 1*genome$Z1 + 2*genome$Z2

# compute Var(IBD)
var.ibd <- ((0 -mean.ibd)^2)*genome$Z0 +
  ((1 -mean.ibd)^2)*genome$Z1 +
  ((2 -mean.ibd)^2)*genome$Z2

# compute SE(IBD)
se.ibd <- sqrt(var.ibd)

# make dataframe with relevant info for plotting
ibd.stat <- data.frame(mean.ibd, se.ibd)

# make the plot
ibd.stat_plot <- ggplot(ibd.stat, aes(mean.ibd, se.ibd)) +
  geom_point() +
  labs(title="", x ="Mean IBD", y = "SE IBD") +
  theme_bw()

# investigate plot
ibd.stat_plot

# save plot
ggsave("output/ibd.stat.pdf", ibd.stat_plot)

# identify two pairs of identical individuals (genome$Z2=1, i.e. they share all their alleles IBD).
duplicate <- genome[mean.ibd == 2,]
duplicate

# save details for one in each pair into a file
fail_ibd_qc <- data.frame(FID=duplicate$FID2, IID=duplicate$IID2)
write.table(fail_ibd_qc, "output/fail_ibd_qc.txt", row.names=F, col.names=T, quote=F)


## Remove all individuals failing QC
#First, concatenate all the files listing individuals failing the previous QC steps into a single file
fail_mis_qc <- read.table("output/fail_mis_qc.txt",header=T,as.is=T)
fail_het_qc <- read.table("output/fail_het_qc.txt", header=T, as.is=T)
fail_ibd_qc <- read.table("output/fail_ibd_qc.txt", header=T, as.is=T)

# combine these
fail_qc <- data.frame(FID=c(fail_mis_qc$FID, fail_het_qc$FID,  fail_ibd_qc$FID), IID=c(fail_mis_qc$IID, fail_het_qc$IID,  fail_ibd_qc$IID))

# find unique individuals, i.e. exclude any duplicates
fail_qc <- unique(fail_qc)

# save details for individuals that fail any of the QC criteria
write.table(fail_qc, "output/fail_qc.txt", row.names=F, col.names=F, quote=F)

#####  SNP Filters
# Similar steps to those above

# Excessive missing data rate
# read in day2.qc.ind.mis.lmiss file
lmis <- read.table("output/qc.ind.mis.lmiss",header=T)

# plot histogram of the fraction of missing genotypes:
mis1_genot_plot <- ggplot(lmis, aes(F_MISS)) +
  geom_histogram(bins=9) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="Fraction of missing data", x ="Fraction of missing genotypes", y = "Number of SNPs") +
  theme_bw()

# make plot, ignore the warning
mis1_genot_plot

# save plot
ggsave("output/qc_ind_lmis_1.pdf", mis1_genot_plot)

# make plot with cut off indicated
mis_genot_plot_2 <- ggplot(lmis, aes(F_MISS)) +
  geom_histogram(bins=9) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="Fraction of missing data", x ="Fraction of missing genotypes", y = "Number of SNPs") +
  theme_bw() +
  geom_vline(xintercept=0.04, linetype="dashed")

# make plot, ignore the warning
mis_genot_plot_2

# save plot
ggsave("output/qc_ind_lmis_2.pdf", mis_genot_plot_2)


##  Examine the frequency of the minor alleles by computing them in 
# read in frequencies file
freq <- read.table("output/qc.ind.freq.frq", header=T)

# plot minor alleles
maf_plot <- ggplot(freq, aes(MAF)) +
  geom_histogram(bins=10) +
  labs(title="Minor allele frequencies", x ="MAF", y = "Number of SNPs") +
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dashed")

# make plot
maf_plot

# save plot
ggsave("output/maf.pdf", maf_plot)

## Different genotype call rates between cases and controls
# read in qc.ind.call.rate.missing file
diffmiss <- read.table("output/qc.ind.call.rate.missing", header=T, as.is=T)

# save SNPs with p-value < 0.000001 into a file
diff.miss <- diffmiss[diffmiss$P<0.000001,]

# save details for those who fail
write.table(diff.miss$SNP, "output/fail_diffmiss_qc.txt", row.names=F, col.names=F, quote=F)


## Investigate the results from the logistic regressions before QC and after is done

# After filtering of individuals and samples lets evaluate results after re-testing associations 
# read in the association test results (before QC)
assoc <- read.table("output/results.assoc.logistic", header=T, as.is=T)

# read in the association test results (after QC)
assoc.qc <- read.table("output/day2.qc.log.assoc.logistic", header=T, as.is=T)

# "Manhattan plots
# make plot of results before QC
QC_before_plot <- ggplot(assoc, aes(BP, P)) +
  geom_point() +  scale_y_continuous(trans=reverselog_trans(base=10),
                                     labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="-log10 p-values (before QC)", x ="position", y = "-log10 p-value") +
  theme_bw()
# investigate plot
QC_before_plot

# make plot of results after QC
QC_after_plot <- ggplot(assoc.qc, aes(BP, P)) +
  geom_point() +  scale_y_continuous(trans=reverselog_trans(base=10),
                                     labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="-log10 p-values (after QC)", x ="position", y = "-log10 p-value") +
  theme_bw()
# investigate plot
QC_after_plot

# Saving plots
ggsave("output/comparison_beforeQC.pdf", QC_before_plot)
ggsave("output/comparison_afterQC.pdf", QC_after_plot)

# Q-Q plots of p-values
# plotting the observed p-values vs. the expected p-values

# observed -log10 p-values (before QC)
p.obs <- -log10(sort(assoc$P,decreasing=F))

# expected -log10 p-values (before QC)
p.exp <- -log10( 1:length(p.obs)/length(p.obs) )

qq.stat <- data.frame(p.obs, p.exp)

qq_plot_beforeQC <- ggplot(qq.stat, aes(p.exp, p.obs)) +
  geom_point() +
  labs(title="QQ plot (before QC)", x ="expected -log10 p-value", y = "observed -log10 p-value") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")

qq_plot_beforeQC

#save plot
ggsave("output/qq_plot_beforeQC.pdf", qq_plot_beforeQC)

#observed -log10 p-values (after QC)
p.obs <- -log10(sort(assoc.qc$P,decreasing=F))

#expected -log10 p-values (after QC)
p.exp <- -log10( 1:length(p.obs)/length(p.obs) )

# make dataframe with relevant info for plotting
qq.stat <- data.frame(p.obs, p.exp)

# make plot
qq_plot_afterQC <- ggplot(qq.stat, aes(p.exp, p.obs)) +
  geom_point() +
  labs(title="QQ plot (after QC)", x ="expected -log10 p-value", y = "observed -log10 p-value") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")

# investigate plot
qq_plot_afterQC

# save plot
ggsave("output/qq_plot_after.pdf", qq_plot_afterQC)

## Identifying regions of association
# identify a set of SNPs for follow-up replication according to p-value threshold -log10 (p-value) > 4 
assoc.qc[-log10(assoc.qc$P)>4,]
