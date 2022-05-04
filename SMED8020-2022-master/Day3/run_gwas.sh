#!/bin/bash

## run SAIGE Step 1 for trait ${trait}
function step1() {
  trait=$1
  #input
  phenoFile=${phenoDir}${trait}_pheno.txt
  plinkFile=${plinkPrefix}
  #output
  outpath1=${outpath}/step1/
  outPrefix=${outpath1}${trait}
  timefile=${outPrefix}_step1.time.txt
  logfile=${outPrefix}_step1.log.txt
  mkdir -p ${outpath1}

  /usr/bin/time -o ${timefile} -v \
    Rscript step1_fitNULLGLMM.R \
    --plinkFile=${plinkFile} \
    --phenoFile=${phenoFile} \
    --phenoCol=${trait} \
    --covarColList=${covars} \
    --sampleIDColinphenoFile=IID \
    --traitType=binary \
    --outputPrefix=${outPrefix} \
    --nThreads=4 \
    --LOCO=FALSE \
  > ${logfile}
}

## run SAIGE Step 2 for traít ${trait}, on chunk ${chunk}
function step2_1chunk() {
  trait=$1
  chunk=$2
  #input
  inpPrefix=${outpath}/step1/${trait}  
  IID=${doseDir}Iid.txt
  dosageFileIn=${doseDir}PART_${chunk}.dose.gz
  #output
  outpath2=${outpath}/step2/
  outPrefix=${outpath2}${trait}
  dosageResults=${outPrefix}_SAIGE_PART_${chunk}.imputed.results.txt
  logfile=${outPrefix}_SAIGE_PART_${chunk}.log.txt
  timefile=${outPrefix}_SAIGE_PART_${chunk}.time.txt
  mkdir -p ${outpath2}

  /usr/bin/time -o ${timefile} -v \
    Rscript --vanilla step2_SPAtests.R \
    --dosageFile ${dosageFileIn} \
    --dosageFileNrowSkip 0 \
    --dosageFileNcolSkip 5 \
    --dosageFilecolnamesSkip "SNPID,CHR,POS,Allele0,Allele1" \
    --sampleFile ${IID} \
    --minMAC 1 \
    --minMAF 0 \
    --GMMATmodelFile ${inpPrefix}.rda \
    --varianceRatioFile ${inpPrefix}.varianceRatio.txt \
    --SAIGEOutputFile ${dosageResults} \
    --IsOutputAFinCaseCtrl T \
    > ${logfile}
}

## run SAIGE Step 2 for traít ${trait}, on 4 chunks in parallel
function step2() {
  trait=$1
  for chunk in {01..04} ; do
    echo $chunk
    step2_1chunk $trait $chunk &
  done
  wait
}

## check that the 4 chunk result files exist
function checkChunksExist() {
  trait=$1
  outPrefix=${outpath}/step2/${trait}
  for chunk in {01..04} ; do
    dosageResults=${outPrefix}_SAIGE_PART_${chunk}.imputed.results.txt
    if [ ! -f ${dosageResults} ] ; then
      echo "Error: Chunk $chunk is missing"
      exit 2
    fi
  done
}

## merge the 4 chunk result files
## output *_SAIGE.imputed.results.txt
function mergeChunks() {
  trait=$1
  #input
  inpPrefix=${outpath}/step2/${trait}
  #ouput
  outpath3=${outpath}/step3/
  mergedfile=${outpath3}${trait}_SAIGE.imputed.results.txt
  mkdir -p ${outpath3}
  
  awk 'FNR==1 && NR!=1{next;}{print}' ${inpPrefix}_SAIGE_PART_{01..04}.imputed.results.txt | sed "s/ /\t/g" > ${mergedfile}
}

## check that the 4 chunk result files exist and merge them
function step3() {
  trait=$1
  checkChunksExist $trait
  mergeChunks $trait
}

## filter on MAF and MAF, correct lambda GC
## Plot manhattan and QQplot
## Output: *SAIGE.imputed.finalResults.txt
## *SAIGE.imputed_Manhattan.png
## *SAIGE.imputed.QQ_Plot.png
function plotGwas() {
  trait=$1
  #input
  inpPrefix=${outpath}/step3/${trait}
  #output
  outpath4=${outpath}/step4/
  outPrefix=${outpath4}${trait}
  mkdir -p ${outpath4}

  Rscript --vanilla plot_GWAS_Results.r \
          --analysis SAIGE \
          --minMAC 3 \
          --minMAF 0 \
          --prefix_pheno ${phenoDir}${trait} \
          --prefix_inp ${inpPrefix}  \
          --prefix_out ${outPrefix}  \
          --format dose \
          &> ${outPrefix}.plot.log.txt
}

## make final results, plot Manhattan and qqplot
function step4() {
  plotGwas $trait
}

#export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/3.6

plinkPrefix=/mnt/scratch/day2/data/genotype/plink/genotyped_PIDS802019_pruned
doseDir=/mnt/scratch/day2/data/genotype/chunks/

phenoDir=/mnt/scratch/day2/data/phenotype/

trait=HighBmiH3
covars="batch,Sex,PC1,PC2,PC3,PC4,BirthYear"

outpath=/mnt/scratch/day2/output/
mkdir -p ${outpath}

#step1 $trait
#step2 $trait
#step3 $trait
#step4 $trait

