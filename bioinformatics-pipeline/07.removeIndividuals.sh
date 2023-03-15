

# remove low quality individuals
 vcftools --remove removeSF.txt --gzvcf filtered_SOSP_358Samples_final.vcf.gz --recode --out filtered_SOSP_319Samples
    # Excluding individuals in 'exclude' list
    # After filtering, kept 319 out of 358 Individuals
    # Outputting VCF file...
    # After filtering, kept 13576969 out of a possible 13576969 Sites
    # Run Time = 8863.00 seconds

bgzip filtered_SOSP_319Samples.recode.vcf
tabix filtered_SOSP_319Samples.recode.vcf.gz

#check sample names of final vcf
bcftools query -l filtered_SOSP_319Samples.recode.vcf.gz > final319Names.txt &
 #all good

# check distribution of ref vs. alt alleles
bcftools view -v snps filtered_SOSP_319Samples.recode.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn > Ref319Alt.txt &
    # 2127106 C       T
    # 2124254 G       A
    # 1885091 T       C
    # 1883773 A       G
    #  579824 G       T
    #  579783 C       A
    #  567556 A       T
    #  567189 T       A
    #  554891 T       G
    #  553942 A       C
    #  541340 C       G
    #  540098 G       C

#check number of SNPs
bcftools stats filtered_SOSP_319Samples.recode.vcf.gz > final319file.stats &
  # 12,504,847

