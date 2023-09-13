## 12 sept 2023
## based on code from Tom Booker

import pandas as pd
from scipy import stats
import numpy as np
import sys, glob, argparse, math
from cyvcf2 import VCF
from collections import Counter

# function to read VCF
def parse_VCF(input_data, population_file, output_file, minimum_allele_count):
	vcf_file = VCF(input_data)

# Read the population file into a pandas DataFrame
	population_df = pd.read_csv(population_file, sep='\s+')

# Create a dictionary mapping Sample IDs to Population IDs
	sample_to_population = dict(zip(population_df['Sample_ID'], population_df['Population_ID']))
	sample_populations= [sample_to_population[samp] for samp in vcf_file.samples]

# Make an ordered list of the non-redundant populations from the sample name list
	sample_populations_unique = []
	for sp in sample_populations:
		if sp not in sample_populations_unique:
			sample_populations_unique.append(sp)

	output_temp = open(output_file + ".freqs.csv", "w")
	header = ",".join(["CHROM","POS"]+ sample_populations_unique)
	output_temp.write(header+"\n")

## Now iterate over all variants in the VCF
	for variant in vcf_file:

		# Initialise the dictionary of allele count vectors
		pop_freqs = {}
		for spu in sample_populations_unique:
			pop_freqs[spu] = [0,0]


		for g, pop_p in zip(variant.genotypes, sample_populations):
			# If a sample's genotype contains "-1" do not count it
			if -1 in g[0:2]:continue

## Count up the number of alt and ref alleles:
			alt = sum(g[0:2])
			ref = 2-alt
		# E.g. homozygous refeerence would be 0
		# E.g. heterozygous would be 1
		# E.g. homozygous alt would be 2

			pop_freqs[pop_p][0] += ref
			pop_freqs[pop_p][1] += alt

# Obtain an ordered list of the populaton allele frequencies
		pop_freqs_list = []
		# For each population...
		for pf in sample_populations_unique:
			# Check if the total number of alleles for the population is greater than the minimum...
			if sum(pop_freqs[pf]) <= minimum_allele_count:
				pop_freqs_list.append(np.nan)
				continue
			# Compute the allele frequency
			freq = pop_freqs[pf][1]/(pop_freqs[pf][0]+pop_freqs[pf][1])

			# Add the result to the frequency list
			pop_freqs_list.append(freq)

		# Store the output in the allele frequency file
		output_line= ",".join(map(str,[variant.CHROM, variant.POS] + pop_freqs_list))
		output_temp.write(output_line+"\n")

# Close the output file after finising the loop
	output_temp.close()

	return sample_populations_unique


def main():
	## Set up command line arguments
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
						help = "The name prefix of the output file")

	parser.add_argument("--minAlleleCount",
					required = False,
					dest = "minAlleleCount",
					type = int,
					help = "The minimum number of alleles used to compute population allele frequnecy. Default = 4",
					default = 4)

	parser.add_argument("--verbose",
						required = False,
						action = "store_true",
						dest = "verbose",
						help = "Give this flag if you want to print info to screen")

	args = parser.parse_args()

# Read in the VCF and compute allele frequencies for each population...
	if args.verbose: print("Computing allele freqs from VCF file and population name file")
	ordered_pop_names = parse_VCF(args.vcf,
							args.population_file,
							args.output,
							args.minAlleleCount)

	if args.verbose: print("Reading environmental data")
	envs = pd.read_csv(args.environments)

	# sort envs Population column into the same order as pop_names
	sorter_index = dict(zip(ordered_pop_names, range(len(ordered_pop_names))))
	envs["rank"] = envs["Population"].map(sorter_index)
	envs_sorted = envs.sort_values(["rank"])

	envs_sorted["rank"]=range(29)

## Define the list of environments that will be analysed
	env_list = ["TD","MSP","SHM","DD_0","DD18"]

	if args.verbose: print("Calculating correlations for each SNP")

	output_file = open(args.output+".correlations.csv", "w")
	header_line = "chrom,pos,numPops,maf,"
	header_line += "".join([env+"_tau,"+env+"_pVal," for env in env_list])
	header_line += "\n"

	output_file.write(header_line)

	freqs = open(args.output + ".freqs.csv")

	count = 0
## Iterate over each SNP
	for raw_row in freqs:
# Process the plain text of the freq file
		row = raw_row.strip("\n").split(",")
# If you're on the header, move one
		if row[0]=="CHROM":continue
		count +=1
		if args.verbose:
			if count%10000==0:print("\tProcessed ",count," variants")

## Make a container for the output
		snp_output = []

## Add info about this SNP to the output...
		snp_output.append(row[0])
		snp_output.append(row[1])

# Grab the allele frequency data and store as a vector
		freq_vector = np.array( row[2:], dtype ="float"  )

# Build a masking array - to remove NaN values
		mask_array = np.invert(np.isnan(freq_vector))

# Add the total number of pops to the output
		snp_output.append(str(mask_array.sum()))

## Calculate the minor allele frequency
		maf = min(freq_vector[mask_array].mean(),1-freq_vector[mask_array].mean())

## If the SNP is invariant after all the filtering and whatnot, ignore it
		if maf ==0:
			continue

## Add the MAF to the output list
		snp_output.append(str(maf))

# For each of the environmental variables...
		for env in env_list:
			env_vector = envs_sorted[env]

## The correlation is computed here...
## Masking any nan entries in the freq vector
			tau, p_value = stats.kendalltau(freq_vector[mask_array], env_vector[mask_array])
## Add the correlation data to the output list
			snp_output.append(str(tau))
			snp_output.append(str(p_value))
# Write the output to the correlation file
		output_file.write( ",".join(snp_output)+"\n" )
	output_file.close()


if __name__ == "__main__":
	main()
