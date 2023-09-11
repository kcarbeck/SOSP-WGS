## katherine carbeck
## 1 sept 2023
## prep input for the wza

# exceute python script (computeCorrelationsForWZA.py) to compute correlations for the WZA
# files needed:
    # vcf = filtered6x_SOSP_316Samples_final_071423.vcf.gz
    # population-file (316samples_pops.txt) in the format of 2 columns "Sample_ID" "Population_ID"
    # environments = 
    # output =
python computeCorrelationsForWZA.py --vcf filtered6x_SOSP_316Samples_final_071423.vcf.gz --population-file 316samples_pops.txt --environments environment_data.csv --output output_results
