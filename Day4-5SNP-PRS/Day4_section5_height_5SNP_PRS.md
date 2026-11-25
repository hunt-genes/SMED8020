# 5. Bonus: A PRS you can calculate by hand (5 SNPs for height)

The genome-wide PRS above is built from tens of thousands of SNPs, so it is hard to *see* what PLINK is actually doing. In this section we build a tiny PRS for height from just **5 SNPs**. Because there are so few, we can calculate the score by hand in R, check that PLINK gives the same answer, and look at the genotypes behind each person's score.

A PRS is just a sum. For person *j*:

$$PRS_j = \sum_{i=1}^{k} \beta_i \times G_{ij}$$

where $k$ is the number of SNPs, $\beta_i$ is the GWAS effect size of SNP *i* (from the base data) and $G_{ij}$ is the number of copies (0, 1 or 2) of the effect allele that person *j* carries at SNP *i* (from the target data).

We will calculate the score in five ways:

| Alternative | Score | Range with 5 SNPs |
|:-:|:-|:-:|
| 1 | Unweighted: count the height-increasing alleles | 0 to 10 |
| 2 | Weighted: sum of $\beta \times$ allele count (by hand in R) | 0 to $2\sum\beta$ |
| 3 | Weighted allele count: Alternative 2 rescaled back to the 0–10 "allele" scale | 0 to 10 |
| 4 | Weighted score from PLINK (`--score ... sum`) | same as 2 |
| 5 | Weighted allele count from PLINK's default average score | same as 3 |

> This section uses files made earlier in the practical: **EUR.clumped**, **Height.QC.Transformed**, the **EUR.QC** genotype files, and the `pheno` object and `null.r2` from Section 3. Run Sections 1–3 first.

## 5.1 Pick 5 independent SNPs

We take the 5 SNPs with the smallest P-values after clumping. Clumping has already made sure they are not in LD with each other (r² < 0.1), so each SNP adds independent information.

To make the unweighted count meaningful, every SNP should be counted in the same direction. We therefore **flip SNPs with a negative BETA**, so the effect allele is always the *height-increasing* allele. Swapping the effect allele changes the sign of BETA, but not its size.

###### RStudio session:
```R
n_snps <- 5

# Index SNPs after clumping, sorted by P-value
clumped <- read.table("EUR.clumped", header=T)
clumped <- clumped[order(clumped$P),]
top <- head(clumped$SNP, n_snps)

# Get alleles and effect sizes from the QC'ed base data
base <- read.table("Height.QC.Transformed", header=T, stringsAsFactors=F)
effects <- base[base$SNP %in% top, c("CHR","BP","SNP","A1","A2","BETA","P")]

# Make the effect allele the height-increasing allele
flip <- effects$BETA < 0
effects$EA    <- ifelse(flip, effects$A2, effects$A1)   # effect allele
effects$OA    <- ifelse(flip, effects$A1, effects$A2)   # other allele
effects$BETA  <- abs(effects$BETA)
effects <- effects[order(effects$P),]
effects

# Files for PLINK
write.table(effects$SNP, "height_5snps.snplist", quote=F, row.names=F, col.names=F)
write.table(effects[,c("SNP","EA")], "height_5snps.allele", quote=F, row.names=F, col.names=F)
write.table(effects[,c("SNP","EA","BETA")], "height_5snps.score", quote=F, row.names=F, col.names=F)
```

**??? note "Which chromosomes are the 5 SNPs on? Which SNP has the largest effect?"**

> The height data in this practical are simulated, so these SNPs are not real height genes. In a real analysis, you would usually choose the lead SNPs from a large published GWAS, for example GIANT. Two examples are SNPs near *HMGA2* and *GDF5*.

## 5.2 Get the genotypes of the 5 SNPs

`--recode A` writes the genotypes as allele counts (0/1/2), one column per SNP. `--recode-allele` tells PLINK to count our **effect allele** rather than its default allele.

###### In terminal:
```bash
plink \
    --bfile EUR.QC \
    --extract height_5snps.snplist \
    --recode A \
    --recode-allele height_5snps.allele \
    --out height_5snps
```

This creates **height_5snps.raw**. Its columns are FID, IID, PAT, MAT, SEX and PHENOTYPE, followed by one column per SNP, named `SNP_allele`. For example, `rs123_A` is the number of A alleles.

###### RStudio session:
```R
raw <- read.table("height_5snps.raw", header=T)
geno <- raw[, 7:ncol(raw)]                  # only the genotype columns
colnames(geno) <- sub("_.*", "", colnames(geno))   # keep just the rsID
geno <- geno[, effects$SNP]                 # same SNP order as 'effects'
head(cbind(raw[,c("FID","IID")], geno))
```

**??? note "Look at the first few people. What do the numbers 0, 1 and 2 mean?"**

## 5.3 Alternative 1 – Unweighted score (allele count)

The simplest score counts how many height-increasing alleles each person carries. Each SNP counts equally.

###### RStudio session:
```R
raw$unw_score <- rowSums(geno)
table(raw$unw_score)
```

**??? note "What is the lowest and highest possible score with 5 SNPs? What range do you see?"**

## 5.4 Alternative 2 – Weighted score by hand

Now each allele is multiplied by its effect size, so SNPs with larger effects count more. Let us first calculate it for **one person**, so every step is visible:

###### RStudio session:
```R
person <- 1
worked <- data.frame(SNP        = effects$SNP,
                     EA         = effects$EA,
                     BETA       = effects$BETA,
                     n_alleles  = as.numeric(geno[person,]))
worked$contribution <- worked$BETA * worked$n_alleles
worked
sum(worked$contribution)   # this is person 1's weighted PRS
```

Then calculate the score for everyone at once. This is a matrix multiplication of the genotype table by the vector of betas:

```R
raw$w_score <- as.vector(as.matrix(geno) %*% effects$BETA)
raw$w_score[person]   # same as the worked example above
```

## 5.5 Alternative 3 – Weighted allele count

The weighted score is in units of BETA, which can be hard to interpret. We can rescale it back to the 0–10 scale of Alternative 1, so that it reads as the number of height-increasing alleles. Each allele is weighted by its relative effect size. To do this, we divide by the average beta:

$$\text{weighted allele count} = PRS_{weighted} \times \frac{k}{\sum \beta}$$

```R
raw$w_allele_count <- raw$w_score * n_snps / sum(effects$BETA)
summary(raw[,c("unw_score","w_score","w_allele_count")])
cor(raw[,c("unw_score","w_score","w_allele_count")])
```

**??? note "Why are w_score and w_allele_count perfectly correlated, but unw_score is not?"**

## 5.6 Alternatives 4 and 5 – Let PLINK do it

Now we calculate the same score with PLINK, using the score file we made in 5.1 (column 1 = SNP, 2 = effect allele, 3 = BETA).

###### In terminal:
```bash
# Alternative 4: sum of BETA x allele count
plink \
    --bfile EUR.QC \
    --score height_5snps.score 1 2 3 sum \
    --out height_5snps_sum

# Alternative 5: PLINK's default score (the average per allele)
plink \
    --bfile EUR.QC \
    --score height_5snps.score 1 2 3 \
    --out height_5snps_avg
```

With `sum`, PLINK reports `SCORESUM`, which should equal our hand-calculated `w_score`. By default, PLINK reports `SCORE`: the sum divided by the number of alleles counted (`CNT` = 2 × number of SNPs = 10). Multiplying it by 2 × n_snps gives back the sum, and then we rescale as in Alternative 3:

###### RStudio session:
```R
# Alternative 4
plink_sum <- read.table("height_5snps_sum.profile", header=T)
# Alternative 5 - total weighted number of alleles from PLINK's average score
plink_avg <- read.table("height_5snps_avg.profile", header=T)
plink_avg$w_t_score2 <- ((plink_avg$SCORE*(n_snps*2))/sum(effects$BETA))*n_snps

check <- merge(raw[,c("FID","IID","w_score","w_allele_count")],
               plink_sum[,c("FID","IID","SCORESUM")], by=c("FID","IID"))
check <- merge(check, plink_avg[,c("FID","IID","w_t_score2")], by=c("FID","IID"))
head(check)
all.equal(check$w_score, check$SCORESUM, tolerance=1e-4)          # Alternative 2 vs 4
all.equal(check$w_allele_count, check$w_t_score2, tolerance=1e-4) # Alternative 3 vs 5
```

**??? note "Do your hand-calculated scores match PLINK's?"**

> If a genotype is missing, PLINK fills in the expected allele count (2 × allele frequency), so the scores may differ slightly for those people. Our QC removed SNPs with more than 1% missingness, so this should rarely happen here.

## 5.7 How much does a 5-SNP PRS tell us about height?

We use the same model as in Section 3: height on the PRS plus sex and 6 PCs. We compare the R² with the null model (`null.r2`).

###### RStudio session:
```R
pheno.5 <- merge(pheno, raw[,c("FID","IID","unw_score","w_score","w_allele_count")],
                 by=c("FID","IID"))
covs <- c("Sex", paste0("PC",1:6))   # check the covariate names with colnames(pheno)

res5 <- NULL
for(s in c("unw_score","w_score","w_allele_count")){
    f <- as.formula(paste("Height ~", s, "+", paste(covs, collapse=" + ")))
    m <- lm(f, data=pheno.5)
    co <- summary(m)$coeff[s,]
    res5 <- rbind(res5, data.frame(Score=s,
                                   BETA=co[1], SE=co[2], P=co[4],
                                   R2=summary(m)$r.squared - null.r2,
                                   row.names=NULL))
}
res5
# Compare with the genome-wide "best-fit" PRS from Section 3
prs.result[which.max(prs.result$R2),]
```

The coefficient for `unw_score` and `w_allele_count` has a simple interpretation: **the difference in height per extra height-increasing allele**.

Finally, plot height by the number of height-increasing alleles:

```R
boxplot(Height ~ unw_score, data=pheno.5,
        xlab="Number of height-increasing alleles (out of 10)",
        ylab="Height", col="lightblue", border="grey40")
means <- tapply(pheno.5$Height, pheno.5$unw_score, mean)
lines(seq_along(means), means, type="b", pch=19, col="firebrick", lwd=2)
table(pheno.5$unw_score)   # number of people in each group
```

**??? note "Is there a trend in height across the allele-count groups? How large is the spread of height within each group compared with the difference between groups?"**

**??? note "How much variance in height does the 5-SNP PRS explain compared with the genome-wide best-fit PRS? Why?"**

> **Take-home message:** a PRS is just a weighted count of alleles. Every PRS method does this; methods differ only in *which* SNPs they include and *what weights* they use (C+T, LDpred, PRS-CS, etc.). For a highly polygenic trait like height, a few top SNPs explain very little variation. Thousands of small effects together explain much more, which is why genome-wide scores perform better.
