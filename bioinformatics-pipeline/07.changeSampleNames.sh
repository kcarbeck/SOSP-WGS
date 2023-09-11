# change sample names in VCF
# katherine carbeck
# 18 march 2023

# first create list of individuals/sample names in 358 individual vcf
bcftools query -l filtered_SOSP_358Samples_031823.recode.vcf.gz > vcfNames.txt

### then edit txt file for new sample names (vcfNames.txt -> vcfNamesNew.txt)

# easiest to use bcftools, where vcfNamesNew.txt has the new sample names:
# s = new sample names, one name per line, in the same order as they appear in the VCF file. Alternatively, only samples which need to be renamed can be listed as "old_name new_name\n" pairs separated by whitespaces, each on a separate line. If a sample name contains spaces, the spaces can be escaped using the backslash character, for example "Not\ a\ good\ sample\ name"
bcftools reheader -s vcfNamesNew.txt --threads 8 filtered_SOSP_358Samples_031823.recode.vcf.gz > SOSP_358Samples_final.vcf.gz &

#check sample names of final vcf
bcftools query -l SOSP_358Samples_final.vcf.gz > finalNames.txt &
   #all good

#check number of SNPs
bcftools stats SOSP_358Samples_final.vcf.gz > final358file.stats &
  #12,504,847

bcftools index SOSP_358Samples_final.vcf.gz &


###############################################################

# change sample names in VCF
# katherine carbeck
# 05 July 2023 - Phred sent updated list of sample IDs 


# first create list of individuals/sample names in vcf
bcftools query -l SOSP_316Samples_final.vcf.gz > vcfNames.txt

### then edit txt file for new sample names (vcfNames.txt -> vcfNamesNew.txt)

# easiest to use bcftools, where vcfNamesNew.txt has the new sample names:
# s = new sample names, one name per line, in the same order as they appear in the VCF file. Alternatively, only samples which need to be renamed can be listed as "old_name new_name\n" pairs separated by whitespaces, each on a separate line. If a sample name contains spaces, the spaces can be escaped using the backslash character, for example "Not\ a\ good\ sample\ name"
bcftools reheader -s vcfNamesNew.txt --threads 8 SOSP_316Samples_final.vcf.gz > SOSP_316Samples_final_070523.vcf.gz  &

#check sample names of final vcf
bcftools query -l SOSP_316Samples_final_070523.vcf.gz > finalNames.txt &
   #all good

#check number of SNPs
bcftools stats SOSP_316Samples_final_070523.vcf.gz > final316.stats &
  # 12,483,733

bcftools index SOSP_316Samples_final_070523.vcf.gz &



###########        14 July 2023         #############


# updating sample names again because montana and cleonensis 
bcftools query -l SOSP_316Samples_final_070523.vcf.gz > vcfNames.txt
bcftools query -l filtered6x_SOSP_316Samples_final_070623.recode.vcf.gz > vcfNames2.txt

bcftools reheader -s vcfNamesNew.txt --threads 8 SOSP_316Samples_final_070523.vcf.gz > SOSP_316Samples_final_071423.vcf.gz &
bcftools reheader -s vcfNamesNew2.txt --threads 8 filtered6x_SOSP_316Samples_final_070623.recode.vcf.gz > filtered6x_SOSP_316Samples_final_071423.vcf.gz &


#check sample names of final vcf
bcftools query -l SOSP_316Samples_final_071423.vcf.gz > finalNames.txt &
bcftools query -l filtered6x_SOSP_316Samples_final_071423.vcf.gz > finalNames2.txt &
   #all good

#check number of SNPs
bcftools stats SOSP_316Samples_final_071423.vcf.gz > final316.stats &
  # 12,483,733
bcftools stats filtered6x_SOSP_316Samples_final_071423.vcf.gz > filtered6x_final316.stats &
  # 11,338,012

bcftools index SOSP_316Samples_final_071423.vcf.gz &
bcftools index filtered6x_SOSP_316Samples_final_071423.vcf.gz &


# copy files to storage folder
cp SOSP_316Samples_final_071423.vcf.gz* /lustre2/home/lc736_0001/song_sparrow/final_vcf &
cp filtered6x_SOSP_316Samples_final_071423.vcf.gz* /lustre2/home/lc736_0001/song_sparrow/final_vcf &