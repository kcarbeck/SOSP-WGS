#!/bin/bash

## 22 August 2023
## katherine carbeck
## calculate LD decay using PopLDdecay
## https://github.com/hewm2008/PopLDdecay/blob/main/Manual.pdf
## DOI:10.1093/bioinformatics/bty875

# install PopLDdecay
git clone https://github.com/hewm2008/PopLDdecay.git 
cd PopLDdecay; chmod 755 configure; ./configure;
make;
mv PopLDdecay  bin/;    #     [rm *.o]


####!  step 1: core program where VCF is the input file and LD stat file is output
# path to the dir containing the population vcf files
subset_vcf_dir="/workdir/kcarbeck"

# clear the command file
[ -f /workdir/kcarbeck/ldDecay/ld_decay_commands.txt ] && rm /workdir/kcarbeck/ldDecay/ld_decay_commands.txt

# calc ld
for vcf_file in "$subset_vcf_dir"/*_subset.vcf.gz; do
    population=$(basename "${vcf_file%_subset.vcf.gz}")
    # gnerate the command for each population
    echo "/workdir/kcarbeck/ldDecay/PopLDdecay/bin/PopLDdecay -InVCF \"$vcf_file\" -MAF 0.05 -MaxDist 300 -OutType 1 -OutStat \"/workdir/kcarbeck/ldDecay/out300kb/${population}_ldDecay.stat.gz\"" >> /workdir/kcarbeck/ldDecay/ld_decay_commands.txt
done

# run in parallel
parallel -j 29 < /workdir/kcarbeck/ldDecay/ld_decay_commands.txt
 

####!  step 2: plot using plot_OnePop.pl 
# the plots from this step don't look great so see R script for better plots
ld_dir="/workdir/kcarbeck/ldDecay/out300kb"

# clear the command file
[ -f /workdir/kcarbeck/ldDecay/ld_decay_plot_commands.txt ] && rm /workdir/kcarbeck/ldDecay/ld_decay_plot_commands.txt

for stat_file in "$ld_dir"/*_ldDecay.stat.gz; do
    population=$(basename "${stat_file%_ldDecay.stat.gz}")
    # gnerate the command for each population
    echo "perl /workdir/kcarbeck/ldDecay/PopLDdecay/bin/Plot_OnePop.pl -inFile \"/workdir/kcarbeck/ldDecay/out300kb/${population}_ldDecay.stat.gz\" -keepR -output \"/workdir/kcarbeck/ldDecay/out300kb/output/${population}\"" >> /workdir/kcarbeck/ldDecay/ld_decay_plot_commands.txt
done

parallel -j 29 < /workdir/kcarbeck/ldDecay/ld_decay_plot_commands.txt

# imgcat /workdir/kcarbeck/ldDecay/out300kb/output/nominate_ON.png

####! Step 3: Better plots in R - see PopLDdecayPlot.R script
