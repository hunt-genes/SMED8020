Meta-analysis steps in detail for your ference.

#### 0. We have [installed METAL](http://csg.sph.umich.edu/abecasis/Metal/download/) using the pre-compiled binaries for you. An updated version of METAL that should work on machines requiring 64-bit is available [here](https://github.com/statgen/METAL/blob/master/README.md).    

Create a .bed file file from GLGC-LDL-preMeta.txt using Linux tools `awk` and `sed`. A [BED file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) is not to be confused with the binary PLINK format .bed, but is a frequently used standard format for genetic data which is required to have chromosome, position start, and position end columns.   

**In the terminal:**
```
awk 'NR > 1 {print $2"\t"$3"\t"$4"\t"$5}' GLGC-LDL-preMeta.txt | sed 's/:/\t/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$4"\t"$5}' > GLCG.hg19.bed
```   

**Web option:**
Upload the `GLGC.hg.bed` file you made [here](http://genome.ucsc.edu/cgi-bin/hgLiftOver) if it's less than 500 mb. Select the genome you're coming from and the genome you're lifting over to.     

**Command line option (ONLY FOR YOUR REFERENCE):**
[Download liftOver](https://hgdownload.soe.ucsc.edu/admin/exe/). 
You can use `wget` like so: `wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver`  
Turn on the executable bit `chmod +x ./filePath/utility_name`. 
Now `./filePath/utility_name` is executable.  

`chmod +x ./liftOver`      

We have installed lifOver here:     
`/mnt/scratch/software/liftOver`   

[Download the map.chain](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) for hg19 to hg38       

`wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz`    

The liftover command requires 4 parameters in this order: 
1) oldFile (in .bed format) 
2) map.chain 
3) newFile (just the name) 
4) unMapped
Execute this command:     
`/mnt/scratch/software/liftOver GLCG.hg19.bed hg19ToHg38.over.chain GLGC.h38.bed GLGC.hg38.unmapped`


Example code to create a file with compatible header is here:     
`join -1 4 -2 2 <(sort -k 4 GLGC.h38.bed) <(sort -k 2 GLGC-LDL-preMeta.txt) | awk -v OFS='\t' '{$5=toupper($5);$6=toupper($6)}1' | awk '{print $0"\t"substr($2, 4)"\t"$4":"$5":"$6}' | awk '{print $0"\t"$16":"$17}'| awk -v OFS='\t' '{print $16, $4, $18, $5, $6, $15, $13, $11, $12, $14}'| sed  '1i\CHR\tPOS38\tSNPID\tAllele1\tAllele2\tAF_Allele2\tN\tBETA\tSE\tp.value' > GLGC-LDL-hg38-preMeta-v2.txt`   

We noted that `join` fails in workbench and this needs to be run directly in the terminal, so we have created and subset the GLGC file so it is quicker to run in the meta-analysis, so use the following file moving forward:
`GLGC-LDL-hg38-preMeta-U.txt`



