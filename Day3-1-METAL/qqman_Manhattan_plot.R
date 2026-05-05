#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressMessages(library(data.table))
suppressMessages(library(qqman))

setwd("/mnt/work/workbench/nikhila/SMED8020/Day3-1-METAL")

# -------------------------
# Input / Output
# -------------------------
infile  <- "LDL_METAL_META1.tbl"
outfile <- "LDL_Manhattan.png"

# -------------------------
# Read data
# -------------------------
dt <- fread(infile)

# -------------------------
# Parse MarkerName
# Format: CHR:POS:REF:ALT
# -------------------------
dt[, c("CHR", "BP", "REF", "ALT") := tstrsplit(MarkerName, ":", fixed=TRUE)]

dt[, CHR := as.integer(CHR)]
dt[, BP  := as.integer(BP)]

# -------------------------
# Rename columns to match qqman
# -------------------------
dt[, SNP := MarkerName]
dt[, P := as.numeric(`P-value`)]

# Remove missing
dt <- dt[!is.na(P)]
dt[P == 0, P := 1e-300]

# -------------------------
# Manhattan plot
# -------------------------
png(outfile, width = 1400, height = 700)

manhattan(dt,
          chr = "CHR",
          bp  = "BP",
          snp = "SNP",
          p   = "P",
          main = "LDL_METAL_Manhattan",
#         cex = 0.6, cex.axis = 0.9,
          col = c("blue4", "orange3"),
          ylim = c(0, 15),
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5),
          chrlabs = as.character(c(1:22)))

dev.off()