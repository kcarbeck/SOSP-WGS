## katherine carbeck
## 1 sept 2023
## prep input for the wza

####### *** ▼ TEST VCF ▼ *** #######
## subset vcf to test code first
bcftools view --header-only filtered6x_SOSP_316Samples_final_071423.vcf.gz > test.vcf
bcftools view --no-header filtered6x_SOSP_316Samples_final_071423.vcf.gz | awk '{printf("%f\t%s\n",rand(),$0);}' | sort -t $'\t'  -T . -k1,1g | head -n 200 | cut -f 2- >> test.vcf

# exceute python script (computeCorrelationsForWZA.py) to compute correlations for the WZA
# files needed:
    # vcf = filtered6x_SOSP_316Samples_final_071423.vcf.gz
    # population-file = 316samples_pops.txt in the format of 2 columns "Sample_ID" "Population_ID"
    # environments = scaledPopEnvAnnual.csv in the format of columns "Population" followed by n environmental vars
    # output = output file for WZA
python computeCorrelationsForWZA.py --vcf test.vcf --population-file 316samples_pops.txt --environments scaledPopEnvAnnual.csv --output /workdir/kcarbeck/out/output_corr
    #! Error:
   #  Traceback (most recent call last):
   #    File "/usr/local/lib64/python3.9/site-packages/pandas/core/indexes/base.py", line 3803, in get_loc
   #      return self._engine.get_loc(casted_key)
   #    File "pandas/_libs/index.pyx", line 138, in pandas._libs.index.IndexEngine.get_loc
   #    File "pandas/_libs/index.pyx", line 165, in pandas._libs.index.IndexEngine.get_loc
   #    File "pandas/_libs/hashtable_class_helper.pxi", line 5745, in pandas._libs.hashtable.PyObjectHashTable.get_item
   #    File "pandas/_libs/hashtable_class_helper.pxi", line 5753, in pandas._libs.hashtable.PyObjectHashTable.get_item
   #  KeyError: 'Sample_ID'

   #  The above exception was the direct cause of the following exception:
 
   #  Traceback (most recent call last):
   #    File "/local/workdir/kcarbeck/computeCorrelationsForWZA.py", line 144, in <module>
   #      main()
   #    File "/local/workdir/kcarbeck/computeCorrelationsForWZA.py", line 96, in main
   #      freqs, pop_names = parse_VCF(args.vcf, args.population_file)
   #    File "/local/workdir/kcarbeck/computeCorrelationsForWZA.py", line 26, in parse_VCF
   #      sample_to_population = dict(zip(population_df['Sample_ID'], population_df['Population_ID']))
   #    File "/usr/local/lib64/python3.9/site-packages/pandas/core/frame.py", line 3804, in __getitem__
   #      indexer = self.columns.get_loc(key)
   #    File "/usr/local/lib64/python3.9/site-packages/pandas/core/indexes/base.py", line 3805, in get_loc
   #      raise KeyError(key) from err
   #  KeyError: 'Sample_ID'


####### *** ▲ TEST VCF ▲ *** #######




# install cyvcf2 bc it's not already on server
pip install cyvcf2

# exceute python script (computeCorrelationsForWZA.py) to compute correlations for the WZA
# files needed:
    # vcf = filtered6x_SOSP_316Samples_final_071423.vcf.gz
    # population-file = 316samples_pops.txt in the format of 2 columns "Sample_ID" "Population_ID"
    # environments = scaledPopEnvAnnual.csv in the format of columns "Population" followed by n environmental vars
    # output = output file for WZA
python computeCorrelationsForWZA.py --vcf filtered6x_SOSP_316Samples_final_071423.vcf.gz --population-file 316samples_pops.txt --environments scaledPopEnvAnnual.csv --output /workdir/kcarbeck/out/output_corr














##
# bcftools view --header-only /workdir/kcarbeck/filtered6x_SOSP_316Samples_final_071423.vcf.gz > /workdir/kcarbeck/out/test.vcf
# bcftools view --no-header /workdir/kcarbeck/filtered6x_SOSP_316Samples_final_071423.vcf.gz | shuf -n 200 >> /workdir/kcarbeck/out/test.vcf