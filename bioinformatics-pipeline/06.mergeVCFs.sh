# merge vcfs to create final vcf of all 358 individuals
# katherine carbeck
# 18 march 2023
# had a few trials of merging vcfs that resulted in a lot of missing data for the added samples. realized i needed to use the --missing-to-ref argument to deal with the missingness problem

# merging requires that the files be indexed
bcftools index filtered_SOSP_6Samples_031323.vcf.gz
bcftools index filtered_SOSP_352Samples_122322.vcf.gz

# bcftools merge 
bcftools merge --merge none --missing-to-ref -Oz filtered_SOSP_6Samples_031323.vcf.gz filtered_SOSP_352Samples_122322.vcf.gz > filtered_SOSP_358Samples_031823.vcf.gz

bcftools index filtered_SOSP_358Samples_031823.vcf.gz

# Filter again
vcftools --gzvcf filtered_SOSP_358Samples_031823.vcf.gz --max-missing 0.8 --maf 0.05 --min-meanDP 2 --max-meanDP 50 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --out filtered_SOSP_358Samples_031823
  # After filtering, kept 358 out of 358 Individuals
  # Outputting VCF file...
  # After filtering, kept 13561802 out of a possible 14774144 Sites
  # Run Time = 10227.00 seconds

bgzip filtered_SOSP_358Samples_031823.recode.vcf
bcftools index filtered_SOSP_358Samples_031823.recode.vcf.gz

#calculate missingness
vcftools --gzvcf filtered_SOSP_358Samples_031823.recode.vcf.gz --missing-indv --out missingness
  # After filtering, kept 358 out of 358 Individuals
  # Outputting Individual Missingness
  # After filtering, kept 13561802 out of a possible 13561802 Sites
  # Run Time = 1027.00 seconds


# distribution of ref vs. alt alleles
bcftools view -v snps filtered_SOSP_358Samples_031823.recode.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn > RefAlt.txt &
  # 2120917 C       T
  # 2117811 G       A
  # 1882601 T       C
  # 1881275 A       G
  # 578995 G       T
  # 578626 C       A
  # 566780 A       T
  # 566435 T       A
  # 555003 T       G
  # 553990 A       C
  # 541283 C       G
  # 540017 G       C


bcftools stats filtered_SOSP_358Samples_031823.recode.vcf.gz > march18file.stats &
  # 12,483,733

