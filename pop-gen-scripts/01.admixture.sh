#!/bin/bash

## admixture
## following https://speciationgenomics.github.io/ADMIXTURE/
## 14 july 2023
## ADMIXTURE is a clustering software similar to STRUCTURE with the aim to infer populations and individual ancestries.



############ get rid of missing data -- maybe don't need this step ###############
# get rid of missing data
# bcftools view -e 'GT[*]="mis"' filtered6x_SOSP_316Samples_final_071423.vcf.gz > filtered6x_SOSP_316Samples.NoMiss.vcf.gz
# check num of snps
# bcftools stats filtered6x_SOSP_316Samples.NoMiss.vcf.gz > snps.stats &
    # 34,919

#bgzip SOSP_316Samples_final.NoMiss.vcf.gz &
#tabix SOSP_316Samples_final.NoMiss.vcf.gz.gz 



############ LD pruning ##############
bash ldpruning.sh filtered6x_SOSP_316Samples_final_071423.vcf.gz & 
bcftools stats filtered6x_SOSP_316Samples_final_071423.LDpruned.vcf.gz > snps.stats 
    # 6,059,018


# Generate the input file in plink format
vcftools --gzvcf filtered6x_SOSP_316Samples_final_071423.LDpruned.vcf.gz --plink --out LDpruned_plink 
    # After filtering, kept 6539274 out of a possible 6539274 Sites
    # Run Time = 717.00 seconds

/programs/plink-1.9-x86_64-beta3.30/plink --file LDpruned_plink --make-bed --allow-extra-chr --out LDpruned_plink &
#LDpruned_plink-temporary.fam written.
    # 6539274 variants loaded from .bim file.
    # 316 people (0 males, 0 females, 316 ambiguous) loaded from .fam.
    # Ambiguous sex IDs written to LDpruned_plink.nosex .
    # Using 1 thread (no multithreaded calculations invoked.
    # Before main variant filters, 316 founders and 0 nonfounders present.
    # Calculating allele frequencies... done.
    # Total genotyping rate is 0.973735.
    # 6539274 variants and 316 people pass filters and QC.
    # Note: No phenotypes present.
######



# 1. ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1=0;print $0}' LDpruned_plink.bim > LDpruned_plink.bim.tmp
mv LDpruned_plink.bim.tmp LDpruned_plink.bim

# 2. run ADMIXTURE in loop
# default is 5-fold cross validation
    # example with just k=1
    # /programs/admixture/admixture --cv=5 SOSP_admixture_out_plink.bed 1
    # -j[number of cores]
    # -- cv default is 5

[ -f /workdir/kcarbeck/admixture/admixtureCommands.txt ] && rm /workdir/kcarbeck/admixture/admixtureCommands.txt

cd /workdir/kcarbeck/admixture

for i in {1..29}
do
 echo "/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed $i > log${i}.out" >> /workdir/kcarbeck/admixture/admixtureCommands.txt
done

parallel -j 5 < /workdir/kcarbeck/admixture/admixtureCommands.txt

    # uses a lot of RAM, ran on 24 core/128 GB RAM cluster
    # took about a week to run, but stalled for k > 22

#try to run k=23-29 seperately to see if theyll complete (25 july 2023)
/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed 23 > log23.out &
/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed 24 > log24.out &
/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed 25 > log25.out 
/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed 26 > log26.out &
/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed 27 > log27.out &
/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed 28 > log28.out 
/programs/admixture/admixture --cv=5 -j2 LDpruned_plink.bed 29 > log29.out &

## identify the best value of k clusters (lowest cross-validation error) 
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > admixture.cv.error
    # 1 0.57670
    # 2 0.54944
    # 3 0.52711
    # 4 0.51932
    # 5 0.51369
    # 6 0.50779
    # 7 0.50794
    # 8 0.50542
    # 9 0.50409
    # 10 0.50549
    # 11 0.50957
    # 12 0.51479
    # 13 0.51728
    # 14 0.51907
    # 15 0.52249
    # 16 0.52772
    # 17 0.53622
    # 18 0.54955
    # 19 0.54804
    # 20 0.55295
    # 21 0.56145
    # 22 0.58397



## To make plotting easier, we can make a file with the individual names in one column and the species names in the second column. As the species name is in the individual name, it is very easy to extract the species name from the individual name:
awk '{split($1,name,"_"); print $1,name[1]}' LDpruned_plink.nosex > LDpruned_plink.list


paste LDpruned_plink.list LDpruned_plink.22.Q

## OUTPUT STORED: /home/lc736_0001/song_sparrow/final_vcf/admixture/admixture

##! NOW PLOT using plotAdmixture.R script
# copy output (.Q files, .list file) from terminal to personal computer
