# random code 

#####  check num of snps in vcf files ####

### 1. all individuals, LD pruned
bcftools stats filtered_SOSP_352Samples_122322.LDpruned.vcf.gz > file.stats
    # number of SNPs = 7,032,298

# show the whole file from the top (use show just the first 10 lines)
bcftools view filtered_SOSP_352Samples_122322.LDpruned.vcf.gz | head

# show just the header with -h.  Here look at just the last 10 lines of the header
bcftools view -h filtered_SOSP_352Samples_122322.LDpruned.vcf.gz | tail


### 2. ON indiv removed, LD pruned
bcftools stats filtered_SOSP_351Samples.LDpruned.vcf.gz > file.stats.LD351
    # number of snps = 12,606,921

# show the whole file from the top (use show just the first 10 lines)
bcftools view filtered_SOSP_351Samples.LDpruned.vcf.gz | head

# show just the header with -h.  Here look at just the last 10 lines of the header
bcftools view -h filtered_SOSP_351Samples.LDpruned.vcf.gz | tail


### 3. full vcf
bcftools stats filtered_SOSP_352Samples_122322.vcf.gz > file.stats.S352
    # number of snps = 12,579,961

# show the whole file from the top (use show just the first 10 lines)
bcftools stats filtered_SOSP_352Samples_122322.vcf.gz | head

# show just the header with -h.  Here look at just the last 10 lines of the header
bcftools view -h stats filtered_SOSP_352Samples_122322.vcf.gz | tail






###

bcftools index data.vcf 
bcftools index -n data.vcf



# from jen
vcftools --vcf Orioles_filtered_final_093020_55individuals.recode.vcf --max-missing 1 --not-chr chrz --recode --stdout | gzip > Orioles.noN.autosomes.vcf.gz
bash ldPruning.sh Orioles.noN.autosomes.vcf.gz

