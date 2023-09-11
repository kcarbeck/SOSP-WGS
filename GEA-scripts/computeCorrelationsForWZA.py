## 1 sept 2023
## based on code from Tom Booker


import pandas as pd
from scipy import stats
import numpy as np
import sys, glob, argparse, math
from cyvcf2 import VCF
from collections import Counter

# functions to round up to nearest 100 or 10,000
def roundup(x):
	return int(math.ceil(x / 100.0)) * 100

def roundup_10000(x):
	return int(math.ceil(x / 10000.0)) * 10000

# function to read VCF 
def parse_VCF(input_data, population_file):
    vcf_file = VCF(input_data)

    # Read the population file into a pandas DataFrame
    population_df = pd.read_csv(population_file, delimiter='\t')

    # Create a dictionary mapping Sample IDs to Population IDs
    sample_to_population = dict(zip(population_df['Sample_ID'], population_df['Population_ID']))

    output_lines = []
    populations = []  # Initialize a list to store populations
    for variant in vcf_file:
#		if len(output_lines) == 1000: break
        populations.append(sample_to_population.get(variant.samples[0]['sample']))  # Get the population for the sample

	    # Enter the filter bed
        pop_freqs = {}
		for g, pop_p in zip(variant.genotypes, populations):

			alt = sum(g[0:2])
			ref = 2-alt

			try:
				pop_freqs[pop_p][0] += ref
				pop_freqs[pop_p][1] += alt
			except KeyError:
				pop_freqs[pop_p] = [0,0]
				pop_freqs[pop_p][0] += ref
				pop_freqs[pop_p][1] += alt

		pop_freqs_list = []
		for pf in pop_freqs.keys():
			freq = pop_freqs[pf][1]/(pop_freqs[pf][0]+pop_freqs[pf][1])
			#	if variant.POS== 113266:print(freq)
			pop_freqs_list.append(freq)
			
	#	print(pop_freqs_list)
        output_lines.append([variant.CHROM, variant.POS, int(roundup_10000(variant.POS)), pop_p] + pop_freqs_list)

	output_header = ["CHROM","POS","WINDOW"]+["p_"+str(p) for p in set(populations)]

	return pd.DataFrame( output_lines ,
							columns = output_header), vcf_file.samples


def main():
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--vcf", "-v",
	                    required = True,
	                    dest = "vcf",
	                    type = str,
	                    help = "The VCF file containing the genotypes")

	parser.add_argument("--population-file", "-p",
                        required=True,
                        dest="population_file",
                        type=str,
                        help="The text file containing sample to population mapping")

	parser.add_argument("--environments","-e",
	                    required = True,
	                    dest = "environments",
	                    help = "The file of environments that you are comparing to the freq data")

	parser.add_argument("--output",
	                    required = True,
	                    dest = "output",
	                    type = str,
	                    help = "The name of the output file (the environment will be prepended to the file name so be sure to write to this dir!)")

	args = parser.parse_args()

	#print(roundup(19600))

#	def roundup(x):
#		return int(math.ceil(x / 100.0)) * 100

	freqs, pop_names = parse_VCF(args.vcf, args.population_file)
	print(pop_names)

	freqs.to_csv(args.output + ".freqs.csv", index = False)

	envs = pd.read_csv(args.environments)

## Assumes that each optimum is unique, that's not sensivle
	env_vector = envs.optimum.unique()
	count = 0
	print("calculating correlations for each SNP")

	output_file = open(args.output+".correlations.csv", "w")
	output_file.write("chrom,pos,window,maf,tau,pVal\n")

## Iterate over each SNP0
	for k, row in freqs.iterrows():

## Make a container for the output
		snp_output = []
## Add info about this SNP to the output...
		snp_output.append(row.CHROM)
		snp_output.append(str(row.POS))
		snp_output.append(str(row.WINDOW))

# Grab the allele ferquncy data and store as a vector
		freq_vector = row[3:]
#		env_vector[[int(i.split("_")[1])-1 for i in freq_vector.index]]

## Calculate the minor allele frequency
		maf = min(freq_vector.mean(),1-freq_vector.mean())
		if maf ==0:
			continue
## Add it to the output list
		snp_output.append(str(maf))

## The correlation is computed here...
		tau, p_value = stats.kendalltau(freq_vector, env_vector)
## Add the correlation data to the output list
		snp_output.append(str(tau))
		snp_output.append(str(p_value))

		output_file.write( ",".join(snp_output) + "\n" )
		count += 1
	output_file.close()


if __name__ == "__main__":
	main()