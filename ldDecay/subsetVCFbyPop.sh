#!/bin/bash

## 21 August 2023
## katherine carbeck
## subset main VCF to population VCFs
## create 29 separate VCF files for the 29 populations using bcftools


# zip and index if not already done
bgzip filtered6x_SOSP_316Samples_final_071423.vcf
tabix -p vcf filtered6x_SOSP_316Samples_final_071423.vcf.gz


####  loop to create population txt files based on pop.txt (also saved as 316samples_pop.txt)  ####
# path to pop.txt file
pop_file="/workdir/kcarbeck/pop.txt"
# create directory to store the population-specific files
mkdir -p pop_files

# loop for each line in the pop.txt file
while read -r sample population; do
    # append the sample to the population-specific file
    echo "$sample" >> "pop_files/${population}.txt"
done < "$pop_file"


###  subsetting loop  ####
[ -f /workdir/kcarbeck/vcf_commands.txt ] && rm /workdir/kcarbeck/vcf_commands.txt

# path to input VCF file
input_vcf="/workdir/kcarbeck/filtered6x_SOSP_316Samples_final_071423.vcf.gz"
# path to dir containing population-specific .txt files
pop_txt_dir="/workdir/kcarbeck/pop_files"

# loop through each population-specific .txt file and generate commands to commands file
for txt_file in "$pop_txt_dir"/*.txt; do
    population=$(basename "${txt_file%.txt}")
    output_vcf="${population}_subset.vcf.gz"
    
    # generate the command for each population
    echo "bcftools view -Oz --samples-file \"$txt_file\" \"$input_vcf\" > \"$output_vcf\"" >> /workdir/kcarbeck/vcf_commands.txt
done

# run in parallel
parallel -j 15 < /workdir/kcarbeck/vcf_commands.txt
    #probably could have run all 29 together; took maybe about an hour?



####  tabix loop  ####
subset_vcf_dir="/workdir/kcarbeck"

[ -f /workdir/kcarbeck/tabix_commands.txt ] && rm /workdir/kcarbeck/tabix_commands.txt
# loop through each subset VCF file and generate tabix commands
for vcf_file in "$subset_vcf_dir"/*_subset.vcf.gz; do
    echo "tabix -p vcf \"$vcf_file\"" >> /workdir/kcarbeck/tabix_commands.txt
done

# run tabix indexing in parallel
parallel -j 29 < /workdir/kcarbeck/tabix_commands.txt 