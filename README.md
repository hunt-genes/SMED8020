
[Precourse materials](SMED_8020-pre-course_preparations_2023.pdf)

# SMED8020 Practicals

[Day 1 - Pen and paper exercises](Day1)

[Day 2 - QC, data preparation and GWAS](Day2)

[Day 3 - PART 1: GWAS meta-analysis](Day3-1-METAL)

[Day 3 - PART 2: FUMA](Day3-2-FUMA)

[Day 4 - PRS calculation and analysis](Day4)

[Day 5 - Mendelian randomization](Day5)

# Required data
`/mnt/scratch/benb/data` 

Copy the data to your directory:    
```
cp -R /mnt/scratch/benb/data/Day3/* /mnt/scratch/{user_name}/smed8020/Day3-1-METAL/    
cp -R /mnt/scratch/benb/data/Day4/* /mnt/scratch/{user_name}/smed8020/Day4/   
```

Note: Change user_name

# Required software
All software is installed in the directory:    
`/mnt/scratch/software/`      

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
