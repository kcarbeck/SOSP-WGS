#!/bin/bash

## linkage pruning for future analysis
## following https://speciationgenomics.github.io/pca/
## 14 march 2023


cd pca

VCF=/workdir/kcarbeck/pca/filtered_SOSP_319Samples_final.vcf.gz

## 1. Generate the input file in plink format 
# commands:
    # --vcf - specified the location of our VCF file.
    # --double-id - told plink to duplicate the id of our samples (this is because plink typically expects a family and individual id - i.e. for pedigree data - this is not necessary for us.
    # --allow-extra-chr - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink expects chromosomes 1-22 and the human X chromosome.
    # --new-id-max-allele-len - whole-genome sequencing frequently contain variants which have not been assigned standard IDs. If you don't want to throw out all of that data, you'll usually want to assign them chromosome-and-position-based IDs.
    # --indep-pairwise - finally we are actually on the command that performs our linkage pruning! The first argument, 50 denotes we have set a window of 50 Kb. The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage. Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.1.
    # --out Produce the prefix for the output data.
/programs/plink2_linux_x86_64_20221024/plink2 --vcf $VCF --allow-extra-chr --double-id --rm-dup --new-id-max-allele-len 51 missing --set-missing-var-ids @:#\$r:\$a --indep-pairwise 50 10 0.1 --out SOSP_319Samples_LDpruning
    # Using up to 24 threads (change this with --threads).
    # --vcf: 13576969 variants scanned.
    # --vcf: SOSP_319Samples_LDpruning-temporary.pgen +
    # SOSP_319Samples_LDpruning-temporary.pvar.zst +
    # SOSP_319Samples_LDpruning-temporary.psam written.
    # 319 samples (0 females, 0 males, 319 ambiguous; 319 founders) loaded from
    # SOSP_319Samples_LDpruning-temporary.psam.
    # 13576969 variants loaded from SOSP_319Samples_LDpruning-temporary.pvar.zst.
    # Note: No phenotype data present.
    # Note: Skipping --rm-dup since no duplicate IDs are present.
    # Calculating allele frequencies... done.
    # --indep-pairwise (13 compute threads): 10534429/13576969 variants removed.
    # Variant lists written to SOSP_319Samples_LDpruning.prune.in and
    # SOSP_319Samples_LDpruning.prune.out .
    # 
    # End time: Tue Mar 14 18:51:20 2023

## 2. prune and create pca
# commands:
    # --extract - this just lets plink know we want to extract only these positions from our VCF - in other words, the analysis will only be conducted on these.
    # --make-bed - this is necessary to write out some additional files for another type of population structure analysis - a model based approach with admixture.
    # --pca - fairly self explanatory, this tells plink to calculate a principal components analysis.

/programs/plink2_linux_x86_64_20221024/plink2 --vcf $VCF --double-id --allow-extra-chr --extract SOSP_319Samples_LDpruning.prune.in --make-bed --pca --out ldPruned_out



vcftools --vcf Orioles_filtered_final_093020_55individuals.recode.vcf --max-missing 1 --not-chr chrz --recode --stdout | gzip > Orioles.noN.autosomes.vcf.gz
bash ldPruning.sh Orioles.noN.autosomes.vcf.gz

plink --bfile DATA --indep-pairwise 50 10 0.8 --out OUTPUT --noweb
plink --bfile DATA --exclude OUTPUT.prune.out --noweb --make-bed --out DATA_FILTERED