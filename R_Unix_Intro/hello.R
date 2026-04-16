#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript hello.R <name> <output_file>")
}

name <- args[1]
outfile <- args[2]

message <- paste("Hello", name, "from R!")

cat(message, file = outfile)
cat(message, "
")
