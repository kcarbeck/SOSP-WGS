#!/bin/bash

## LD Decay
## following https://speciationgenomics.github.io/ld_decay/
## 21 feb 2023


# move to your home directory
cd ~
# make a plink directory
mkdir ld_decay
# move into it
cd ld_decay

# 1. calculate LD using plink
# commands:
    # --vcf - specified the location of our VCF file.
    # --double-id - told plink to duplicate the id of our samples (this is because plink typically expects a family and # individual id - i.e. for pedigree data - this is not necessary for us.
    # --allow-extra-chr - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink # expects chromosomes 1-22 and the human X chromosome.
    # --maf - this filters based on minor allele frequency - 0.01 in this case.
    # --geno - this filters out any variants where more than X proportion of genotypes are missing data. Here we throw out # anything with >10% missing dataset.
    # --thin - thin randomly thins out the data - i.e. it randomly retains p proportion of the data. Here we set that to 0.1 or 10%. This means that each site has a 10% probability of remaining in the dataset. This is done to ensure the analysis runs quickly and our output isn’t too big as it can quickly get out of hand!
    # --r2 - finally we’re on to the options for LD! This tells plink to produce squared correlation coefficients. We also provide the argument gz in order to ensure the output is compressed. This is very important as it is easy to produce EXTREMELY large files.
    # --ld-window - this allows us to set the size of the lower end of the LD window. In this case, we set it to 100 bp - i.e. any sites with < 100 sites between them are ignored.
    # --ld-window-kb - this is the upper end of the LD window. Here we set it to 1000, meaning that we ignore any two sites more than 1 Mb apart in the genome.
    # --ld-window-r2 - the final LD command - this sets a filter on the final output but we want all values of LD to be written out, so we set it to 0.
    # --make-bed - this just makes a plink bedfile for future analyses
    # --out Produce the prefix for the output data.

VCF=/workdir/kcarbeck/filtered_SOSP_352Samples_122322.vcf.gz

/programs/plink_linux_x86_64_20221210/plink --vcf $VCF --double-id --allow-extra-chr --maf 0.01 --geno 0.1 --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out ldDecay
    #OUTPUT:
    #255921 MB RAM detected; reserving 127960 MB for main workspace.
    #--vcf: ldDecay-temporary.bed + ldDecay-temporary.bim + ldDecay-temporary.fam
    #written.
    #13663084 variants loaded from .bim file.
    #352 people (0 males, 0 females, 352 ambiguous) loaded from .fam.
    #Ambiguous sex IDs written to ldDecay.nosex .
    #--thin: 12297023 variants removed (1366061 remaining).
    #Using up to 63 threads (change this with --threads).
    #Before main variant filters, 352 founders and 0 nonfounders present.
    #Calculating allele frequencies... done.
    #Total genotyping rate is 0.951407.
    #81169 variants removed due to missing genotype data (--geno).
    #0 variants removed due to minor allele threshold(s)
    #(--maf/--max-maf/--mac/--max-mac).
    #1284892 variants and 352 people pass filters and QC.
    #Note: No phenotypes present.
    #--make-bed to ldDecay.bed + ldDecay.bim + ldDecay.fam ... done.
    #--r2 to ldDecay.ld.gz ... done.


# 2. Calculating average LD across set distances
# using ld_decay_calc.py script from Joanna Meier

# run mark python script
module load python/2.7.15
python 03.ldDecayCalc.py -i ldDecay.ld.gz -o ldDecayOut
     # This will take a few moments to run and will produce two files, the pairwise SNPs arranged into distance classes and the average LD across 100 Kb bins. We’ll take the latter and then plot it in R to see the LD decay.

     #OUTPUT FILES:
        # ldDecayOut.ld_decay
        # ldDecayOut.ld_decay_bins



### R ###
# 3. LD decay plotting script
