
[Precourse materials](SMED_8020-pre-course_preparations_2023.pdf)

# SMED8020 Practicals

[Day 1 - Pen and paper exercises](Day1)

[Day 2 - QC, data preparation and GWAS](Day2)

[Day 3 - GWAS meta-analysis](Day3)

[Day 4 - PRS calculation and analysis](Day4)

[Day 5 - Mendelian randomization](Day5)

# Required data
cp -R /mnt/scratch/benb/data/Day2 /mnt/scratch/user/smed8020/Day2     

/mnt/scratch/benb/data/Day2    
/mnt/scratch/benb/data/Day3    
/mnt/scratch/benb/data/Day4    

[Download here](https://ntnu.box.com/s/d74fob6vo86834tuvtbesrt3hjqih0sh)

# Required software
All software is installed    

/mnt/scratch/benb/software/     

1. [R](https://www.r-project.org/) (**version 4.0.0+**)
2. PLINK 1.9 From zip archive for Practical Day 2 here or directly from: https://www.cog-genomics.org/plink2  
3. For Windows users: Ubuntu from Microsoft Store on your computer
4. [METAL](http://csg.sph.umich.edu/abecasis/metal/download/) from University of Michigan Center for Statistical genetics 

# Linux cheet sheet
| Command | Meaning | Description|
|:-:|:-:|:-|
| pwd | Print working directory | Prints path to current directory (folder) |
| . | Here | Exchangeable with the path printed by ‘pwd’ |
| .. | One level up | Gives the path until outside the current folder |
| geno | 0.01 | Excludes SNPs that are missing in a high fraction of subjects. A two-stage filtering process is usually performed (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).|
| mind | 0.01 | Excludes individuals who have a high rate of genotype missingness, since this may indicate problems in the DNA sample or processing. (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/) for more details).|
| make-just-fam | - | Informs `plink` to only generate the QC'ed sample name to avoid generating the .bed file.  |
| write-snplist | - | Informs `plink` to only generate the QC'ed SNP list to avoid generating the .bed file. |
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |
