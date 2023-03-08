#!/bin/bash

## admixture
## following https://speciationgenomics.github.io/ADMIXTURE/
## 21 feb 2023
## ADMIXTURE is a clustering software similar to STRUCTURE with the aim to infer populations and individual ancestries.

# 1. ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1=0;print $0}' ldPrunedVCF.bim > ldPrunedVCF.bim.tmp
mv ldPrunedVCF.bim.tmp ldPrunedVCF.bim

# 2. run ADMIXTURE in loop
# default is 5-fold cross validation
# loop for 8 populations -- we have 8 we sampled from
    # example with just k=1
    # /programs/admixture/admixture --cv=5 SOSP_admixture_out_plink.bed 1
    # -j[number of cores]
    # -- cv default is 5

[ -f /workdir/kcarbeck/admix/admixtureCommands.txt ] && rm /workdir/kcarbeck/admix/admixtureCommands.txt

cd /workdir/kcarbeck/admix

for i in {1..21}
do
 echo "/programs/admixture/admixture --cv=5 -j2 ldPrunedVCF.bed $i > log${i}.out" >> /workdir/kcarbeck/admix/admixtureCommands.txt
done

parallel -j 10 < /workdir/kcarbeck/admix/admixtureCommands.txt

    # uses a lot of RAM, ran on 64 core/256 GB RAM cluster



## identify the best value of k clusters (lowest cross-validation error)
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > admixture.cv.error

    # (K=1): 0.30681
    # (K=2): 0.44270
    # (K=3): 0.47874
    # (K=4): 0.43771
    # (K=5): 0.48460
    # (K=6): 0.39155
    # (K=7): 0.38747
    # (K=8): 0.40383
    # (K=9): 0.58485
    # (K=10): 0.58806
    # (K=11): 0.60077
    # (K=14): 0.63035
    # (K=20): 0.29449

grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > $FILE.cv.error
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > $FILE.cv.error


## To make plotting easier, we can make a file with the individual names in one column and the species names in the second column. As the species name is in the individual name, it is very easy to extract the species name from the individual name:
awk '{split($1,name,"_"); print $1,name[1]}' filtered_SOSP_352Samples_122322.nosex > filtered_SOSP_352Samples_122322.list

paste filtered_SOSP_352Samples_122322.list filtered_SOSP_352Samples_122322.20.Q
    # All the pops are in different groups.


## NOW PLOT 
# using plotAdmixture.R script:

Rscript plotADMIXTURE.r -p $FILE -i $FILE.list -k 5 -l PunNyerMak,PunPundMak,PunNyerPyt,PunHybrPyt,PunPundPyt
