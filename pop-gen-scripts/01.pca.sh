#!/bin/bash

## generate files for pca and admixture
## following https://speciationgenomics.github.io/pca/
## 21 feb 2023

# One of the major assumptions of PCA is that the data we use is indpendent - i.e. there are no spurious correlations among the measured variables. This is obviously not the case for most genomic data as allele frequencies are correlated due to physical linkage and linkage disequilibrium. So as a first step, we need to prune our dataset of variants that are in linkage.


## 1. Generate the input file in plink format 
# commands:
    # --vcf - specified the location of our VCF file.
    # --double-id - told plink to duplicate the id of our samples (this is because plink typically expects a family and individual id - i.e. for pedigree data - this is not necessary for us.
    # --allow-extra-chr - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink expects chromosomes 1-22 and the human X chromosome.
    # --indep-pairwise - finally we are actually on the command that performs our linkage pruning! The first argument, 50 denotes we have set a window of 50 Kb. The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage. Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.1.
    # --out Produce the prefix for the output data.

VCF=/workdir/kcarbeck/filtered_SOSP_352Samples_122322.vcf.gz

/programs/plink_linux_x86_64_20221210/plink --vcf $VCF --double-id --allow-extra-chr --indep-pairwise 50 10 0.1 --out ldPrunedVCF

    # As well as being versatile, plink is very fast. It will quickly produce a linkage analysis for all our data and write plenty of information to the screen. When complete, it will write out two files: prune.in and prune.out. The first of these is a list of sites which fell below our linkage threshold - i.e. those we should retain. The other file is the opposite of this. In the next step, we will produce a PCA from these linkage-pruned sites.


## 2. prune and create pca
# commands:
    # --extract - this just lets plink know we want to extract only these positions from our VCF - in other words, the analysis will only be conducted on these.
    # --make-bed - this is necessary to write out some additional files for another type of population structure analysis - a model based approach with admixture.
    # --pca - fairly self explanatory, this tells plink to calculate a principal components analysis.

/programs/plink_linux_x86_64_20221210/plink --vcf $VCF --double-id --allow-extra-chr --extract ldPrunedVCF.prune.in --make-bed --pca --out ldPrunedVCF
    
    # Logging to ldPrunedVCF.log.
    # Options in effect:
    #   --allow-extra-chr
    #   --double-id
    #   --extract ldPrunedVCF.prune.in
    #   --make-bed
    #   --out ldPrunedVCF
    #   --pca
    #   --vcf /workdir/kcarbeck/filtered_SOSP_352Samples_122322.vcf.gz
    # 
    # 255921 MB RAM detected; reserving 127960 MB for main workspace.
    # --vcf: ldPrunedVCF-temporary.bed + ldPrunedVCF-temporary.bim +
    # ldPrunedVCF-temporary.fam written.
    # 13663084 variants loaded from .bim file.
    # 352 people (0 males, 0 females, 352 ambiguous) loaded from .fam.
    # Ambiguous sex IDs written to ldPrunedVCF.nosex .
    # --extract: 13663084 variants remaining.
    # Warning: At least 3221691 duplicate IDs in --extract file.
    # Using up to 63 threads (change this with --threads).
    # Before main variant filters, 352 founders and 0 nonfounders present.
    # Calculating allele frequencies... done.
    # Total genotyping rate is 0.951422.
    # 13663084 variants and 352 people pass filters and QC.

    # output: 
    # .eigenval and .eigenvec for pca
    # .bed = binary file necessary for admixture analysis. It is essentially the genotypes of the pruned dataset recoded as 1s and 0s.
    # .bim = a map file (i.e. information file) of the variants contained in the bed file.
    # .fam = a map file for the individuals contained in the bed file.


### 3. PCA in R







