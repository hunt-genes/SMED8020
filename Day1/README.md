# DAY 1 - Gitlab orientation and laptop setup for practicals

### FOCUS and LEARNING GOALS
> The aim for this session is for you to 1) become familar with gitlab and how to navigate the practicals, 2) 
> download and install the required programs, and 3) test the scripts used in the practicals.

### Required downloads   

R version 4.0.0 https://cran.uib.no/   
RStudio https://www.rstudio.com/products/rstudio/download/   

If you are using windows, please install bash to run the suggested commands in the terminal:
https://itsfoss.com/install-bash-on-windows/   

You can use Plink from the zip archive here: https://github.com/hunt-genes/SMED8020/tree/main/Day2 
or from the originating site for PLINK: (PLINK 1.9) https://www.cog-genomics.org/plink/ 

### Specific installations and tests for practicals

#### Practical 2 - Setup for PLINK QC: 

Download the zip archive in https://github.com/hunt-genes/SMED8020/tree/main/Day2

Store the files on your local computer and direct the terminal to them by writing ```cd filepath```
This is how you make the terminal know where your files are/what directory to work from
and you can copy this and paste it into your terminal.
For example:
```cd /mnt/c/Users/name00/Desktop/SMED8020/2021/Day2/```

The suggested structure would be a folder that has Plink, the three data day2 files and a subfolder called output.
Note that the Plink.exe file needs to have the same file path as the data

If you are new to linux commands, you could make yourself familiar with basic commands in sites like: https://www.hostinger.com/tutorials/linux-commands

###### In R session:
Libraries needed for R are ggplot2 and scales   

```
install.packages("ggplot2", "scales")
```

#### Practical 5 - Mendelian randomization: https://github.com/hunt-genes/SMED8020/tree/main/SMED8020-2022-master

Test your setup   
-Download the files for practical 5 from github: https://github.com/hunt-genes/SMED8020/tree/main/SMED8020-2022-master 
-Open Rstudio or R   
-Open  Mendelian_Randomization_Practical_1_-_Single_Sample_-_Rscript.R   
-Edit input directories   
-Run   
-Open  Mendelian_Randomization_Practical_2_-_MR_base_-_Rscript.R   
-Edit input directories   
-Run   
