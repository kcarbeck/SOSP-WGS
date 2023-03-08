# prep vcf for gradient forest
# author: katherine carbeck
# 25 jan 2023

/usr/bin/ssh kcarbeck@cbsulogin.biohpc.cornell.edu

# first need to phase and impute missing data using beagle:
# https://faculty.washington.edu/browning/beagle/beagle.html
# tuitorial: https://github.com/adrianodemarino/Imputation_beagle_tutorial

# Other things about which to be careful:
    # Beagle might not like compressed VCFs - you will have to decompress .gz files. But if Beagle accepts compressed or uncompressed VCFs, then your compressed VCFs should be compressed with bgzip and then tab-indexed with tabix -p vcf Variant.vcf.gz


cd /workdir/kcarbeck #workdir

#java -Xmx150g -jar /programs/beagle41/beagle41.jar gt=Scaffold_1_filtered.recode.vcf out=Scaffold_1_imputed

# Run beagle for each file in vcf list (split by region)
for i in $(cat vcf_file_list.txt)
do
STEM=$(echo ${i} | cut -f 1 -d ".")
java -Xmx240g -jar /programs/beagle41/beagle41.jar gt=${i} out=${STEM}_imputed
done

#####
# index first
for imputedVCF in *imputed.vcf.gz; do
    echo $imputedVCF
    tabix $imputedVCF &
done

# then, bcftools concatenate
bcftools concat regions_01_imputed.vcf.gz regions_02_imputed.vcf.gz regions_03_imputed.vcf.gz regions_04_imputed.vcf.gz regions_05_imputed.vcf.gz regions_06_imputed.vcf.gz regions_07_imputed.vcf.gz regions_08_imputed.vcf.gz regions_09_imputed.vcf.gz regions_10_imputed.vcf.gz regions_11_imputed.vcf.gz regions_12_imputed.vcf.gz regions_13_imputed.vcf.gz regions_14_imputed.vcf.gz regions_15_imputed.vcf.gz regions_16_imputed.vcf.gz regions_17_imputed.vcf.gz regions_18_imputed.vcf.gz regions_19_imputed.vcf.gz regions_20_imputed.vcf.gz regions_21_imputed.vcf.gz regions_22_imputed.vcf.gz regions_23_imputed.vcf.gz regions_24_imputed.vcf.gz regions_25_imputed.vcf.gz regions_26_imputed.vcf.gz regions_27_imputed.vcf.gz regions_28_imputed.vcf.gz regions_29_imputed.vcf.gz regions_30_imputed.vcf.gz -Oz -o imputed_SOSP_352Samples_012923.vcf.gz &

# convert to 0,1,2 format for R
vcftools --gzvcf imputed_SOSP_352Samples_012923.vcf.gz --012 --out snps_for_R &




