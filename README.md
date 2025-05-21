
[Precourse materials](SMED_8020-pre-course_preparations_2023.pdf)

# SMED8020 Practicals

[Day 1 - Pen and paper exercises](Day1)

[Day 2 - QC, data preparation and GWAS](Day2)

[Day 3 - PART 1: GWAS meta-analysis](Day3-1-METAL)

[Day 3 - PART 2: FUMA](Day3-2-FUMA)

[Day 4 - PRS calculation and analysis](Day4)

[Day 5 - Mendelian randomization](Day5)

# Required data
Go to your home directory

#Note: $USER is a variable that is equal to your username
```   
cd /mnt/work/workbench/$USER
```

Clone the repository
```
git clone https://github.com/hunt-genes/SMED8020.git
```

Copy the data to your directory:    

#Note: $USER is a variable that is equal to your username
```
cp -R /mnt/scratch/data/Day3/* /mnt/work/workbench/$USER/SMED8020/Day3-1-METAL/
```

#Note: $USER is a variable that is equal to your username 
```
cp -R /mnt/scratch/data/Day4/* /mnt/work/workbench/$USER/SMED8020/Day4/
```

# Required software
All software is installed in the directory:    
`/mnt/scratch/software/`

###### In terminal:
We need to tell the terminal where plink is:    
``` 
#Note: $USER is a variable that is equal to your username
cd /mnt/work/workbench/$USER
```  

```  
cp -v /mnt/scratch/software/plink/plink /mnt/work/workbench/$USER/.local/bin/plink
```  

###### In RStudio Server session:
```   
install.packages("AER")
install.packages('stringi')
install.packages("devtools")
library(devtools)
install_github("MRCIEU/TwoSampleMR")
install.packages("plyr")
install.packages("ggplot2")
```  

# Linux cheet sheet
| Command | Meaning | Description|
|:-:|:-:|:-|
| pwd | Print working directory | Prints path to current directory (folder) |
| . | Here | Exchangeable with the path printed by ‘pwd’ |
| .. | One level up | Gives the path until outside the current folder |
| cd <path> | Change directory | Moves you to the specified path |
| TAB button | Complete or list options | Completes path/filename that you have started typing, or gives you available options if multiple matches |
| ls | list | Lists file in the current directory |
| head <file> | - | Prints first 10 lines of the specified file |
| tail <file> | - | Prints last 10 lines of the specified file |
| less <file> | - | Opens file for you to read/scroll. (exit with ‘q’) |
| wc –l <file> | - | Using the –l flag, it counts lines in the specified file |
  
## Lecture notes
https://ntnu.box.com/s/b8ayvovfccwejemcsaqpdvf6o3j02a0m
