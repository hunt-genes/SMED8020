echo "#THIS SCRIPT EXECUTES AN ANALYSIS OF THREE STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES HUNT-LDL-preMeta.txt GLGC-LDL-preMeta.txt BBJ-LDL-preMeta.txt

#LOAD THE THREE INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# UNCOMMENT THE NEXT LINE TO ENABLE frequencies to be averaged
AVERAGEFREQ ON

# === DESCRIBE AND PROCESS THE INPUT FILES ===
MARKER SNPID
ALLELE Allele2 Allele1
EFFECT BETA
PVALUE p.value 
WEIGHT N
FREQLABEL AF_Allele2

PROCESS ${1}
PROCESS ${2}
PROCESS ${3}

OUTFILE ${4} .tbl
MINWEIGHT 10000
SCHEME STDERR
ANALYZE 

QUIT"
