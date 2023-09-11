# remove low quality individuals (SF samples and AK_kenaiensis_S23, OR_montana_S49, ON_melodia_S77)
 vcftools --remove removeSamps.txt --gzvcf SOSP_358Samples_final.vcf.gz --recode --out SOSP_316Samples &
    # After filtering, kept 316 out of 358 Individuals
    # Outputting VCF file...

# compress and index vcf
# can use this to index as well but seems like it takes longer: tabix filtered_SOSP_312Samples.recode.vcf.gz
bgzip SOSP_316Samples.recode.vcf 
bcftools index SOSP_316Samples.recode.vcf.gz





