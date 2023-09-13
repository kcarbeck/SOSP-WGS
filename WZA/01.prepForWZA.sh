## katherine carbeck
## 12 sept 2023
## prep input for the wza

####### *** ▼ TEST VCF ▼ *** #######
## subset vcf for 200 lines using awk sort to test code first
bcftools view --header-only filtered6x_SOSP_316Samples_final_071423.vcf.gz > test.vcf
bcftools view --no-header filtered6x_SOSP_316Samples_final_071423.vcf.gz | awk '{printf("%f\t%s\n",rand(),$0);}' | sort -t $'\t'  -T . -k1,1g | head -n 200 | cut -f 2- >> test.vcf

## subset vcf using shuf for 20,000 lines
bcftools view --header-only filtered6x_SOSP_316Samples_final_071423.vcf.gz > test.vcf
bcftools view --no-header filtered6x_SOSP_316Samples_final_071423.vcf.gz | shuf -n 20000 >> test.vcf

# exceute python script (computeCorrelationsForWZA.py) to compute correlations for the WZA
# files needed:
    # vcf = filtered6x_SOSP_316Samples_final_071423.vcf.gz
    # population-file = 316samples_pops.txt in the format of 2 columns "Sample_ID" "Population_ID"
    # environments = scaledPopEnvAnnual.csv in the format of columns "Population" followed by n environmental vars
    # minAlleleCount = min number of alleles used to compute population allele frequnecy. Default = 4
    # output = output file for WZA
python computeCorrelationsForWZA.py --vcf test.vcf --population-file 316samples_pops.txt --environments scaledPopEnvAnnual.csv --minAlleleCount 4 --verbose --output /workdir/kcarbeck/out/output_corr
    # took a few seconds for the test file with 20,000 SNPs

####### *** ▲ TEST VCF ▲ *** #######




# install cyvcf2 bc it's not already on server
pip install cyvcf2

# exceute python script (computeCorrelationsForWZA.py) to compute correlations for the WZA
# files needed:
    # vcf = filtered6x_SOSP_316Samples_final_071423.vcf.gz
    # population-file = 316samples_pops.txt in the format of 2 columns "Sample_ID" "Population_ID"
    # environments = scaledPopEnvAnnual.csv in the format of columns "Population" followed by n environmental vars
    # minAlleleCount = min number of alleles used to compute population allele frequnecy. Default = 4
    # output = output file for WZA
python computeCorrelationsForWZA.py --vcf filtered6x_SOSP_316Samples_final_071423.vcf.gz --population-file 316samples_pops.txt --environments scaledPopEnvAnnual.csv --minAlleleCount 4 --verbose --output /workdir/kcarbeck/out/output_corr
    # started running 12 sept 2023 around 4:00pm
    # was finished in the morning 13 sept 

#output files:
    # output_corr.correlations.csv = 2.9G
    # output_corr.freqs.csv = 4.0G

# files are quite large so gzip compress:
gzip *.csv
    # gzipped files are now 1.2G and 490M
