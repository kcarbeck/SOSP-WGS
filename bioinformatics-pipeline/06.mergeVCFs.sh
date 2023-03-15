# merge vcfs to create final vcf of all 358 individuals

# merge vcf of 352 individuals with new batch of 6 individuals
tabix -f filtered_SOSP_6Samples_031323.vcf.gz &
tabix -f filtered_SOSP_352Samples_122322.vcf.gz 
bcftools merge --merge all filtered_SOSP_6Samples_031323.vcf.gz filtered_SOSP_352Samples_122322.vcf.gz -o merged_filtered_SOSP_358Samples_031323.vcf


# Filter AGAIN
vcftools --vcf merged_filtered_SOSP_358Samples_031323.vcf --max-missing 0.8 --maf 0.05 --min-meanDP 2 --max-meanDP 50 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --out filtered_SOSP_358Samples_031323

bgzip filtered_SOSP_358Samples_031323.recode.vcf
tabix filtered_SOSP_358Samples_031323.recode.vcf.gz

# distribution of ref vs. alt alleles
bcftools view -v snps filtered_SOSP_358Samples_031323.recode.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn > RefAlt.txt
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


########## change sample names in VCF ###########
# first create list of individuals/sample names in 358 individual vcf
bcftools query -l filtered_SOSP_358Samples_031323.recode.vcf.gz > vcfNames.txt

### then edit txt file for new sample names (vcfNames.txt -> vcfNamesNew.txt)

# easiest to use bcftools, where vcfNamesNew.txt has the new sample names:
# s = new sample names, one name per line, in the same order as they appear in the VCF file. Alternatively, only samples which need to be renamed can be listed as "old_name new_name\n" pairs separated by whitespaces, each on a separate line. If a sample name contains spaces, the spaces can be escaped using the backslash character, for example "Not\ a\ good\ sample\ name"
bcftools reheader -s vcfNamesNew.txt --threads 8 filtered_SOSP_358Samples_031323.recode.vcf.gz > filtered_SOSP_358Samples_final.vcf.gz

#check sample names of final vcf
bcftools query -l filtered_SOSP_358Samples_final.vcf.gz > finalNames.txt
   #all good

#check number of SNPs
bcftools stats filtered_SOSP_358Samples_final.vcf.gz > final358file.stats &
  #12,504,847


tabix filtered_SOSP_358Samples_final.vcf.gz &







####### merging with bcftools #########
# merging requires that the files be indexed
bcftools index filtered_SOSP_6Samples_031323.vcf.gz
bcftools index filtered_SOSP_352Samples_122322.vcf.gz

# merge those into a file 
bcftools merge -Oz filtered_SOSP_6Samples_031323.vcf.gz filtered_SOSP_352Samples_122322.vcf.gz > bcftools_merge.vcf.gz

bcftools index bcftools_merge.vcf.gz

# Filter AGAIN
vcftools --gzvcf bcftools_merge.vcf.gz --max-missing 0.8 --maf 0.05 --min-meanDP 2 --max-meanDP 50 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --out filtered_SOSP_358Samples_031523 

tabix bcftools_merge.vcf.gz




bgzip filtered_SOSP_358Samples_031323.recode.vcf