# Filtering VCF
# more stringent filtering with 316 sample file to reduce number of SNPs for GF analysis and only keep high quality/coverage SNPs
# katherine carbeck
# 06 July 2023

# increase filtering for 6x coverage
vcftools --gzvcf SOSP_316Samples_final_070523.vcf.gz --max-missing 0.8 --maf 0.05 --min-meanDP 6 --max-meanDP 50 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --out filtered6x_SOSP_316Samples_final_070623
    # After filtering, kept 316 out of 316 Individuals
    # Outputting VCF file...
    # After filtering, kept 12,264,573 out of a possible 13,561,802 Sites
    # Run Time = 8024.00 seconds


bgzip filtered6x_SOSP_316Samples_final_070623.recode.vcf
bcftools index filtered6x_SOSP_316Samples_final_070623.recode.vcf.gz


# distribution of ref vs. alt alleles
bcftools view -v snps filtered6x_SOSP_316Samples_final_070623.recode.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn > RefAlt.txt &
  GNU nano 5.6.1                                             RefAlt.txt
    # 1942672 C       T
    # 1939896 G       A
    # 1722480 T       C
    # 1720414 A       G
    #  517103 G       T
    #  516832 C       A
    #  509948 A       T
    #  509894 T       A
    #  493643 A       C
    #  491691 T       G
    #  487122 C       G
    #  486317 G       C


bcftools stats filtered6x_SOSP_316Samples_final_070623.recode.vcf.gz > Jul06file.stats &
  # 11,338,012 SNPs


