# split data by each population for leave one out cross validation
# katherine carbeck
# 17 jan 2024

# maf file = WZA_Outlier_SNPs_MAFs.csv
# env file = scaledPopEnvAnnualUpdatedwLATLONG.csv

#! bash
export PYTHONPATH=/home/kcarbeck/.local/lib64/python3.9/site-packages:$PYTHONPATH
python
#! 

import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneOut
import os

# load dfs
snp = pd.read_csv("WZA_Outlier_SNPs_MAFs.csv")
env = pd.read_csv("scaledPopEnvAnnualUpdatedwLATLONG.csv")

# set 'population' column as the index
snp.set_index('population', inplace=True)
env.set_index('Population', inplace=True)

# Assuming 'population' column is common to both dataframes
populations = snp.index.unique()

# Initialize LeaveOneOut
loo = LeaveOneOut()

# Directory for output files
output_dir = "loocv"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Iterate through each population
for train_index, test_index in loo.split(populations):
    train_populations = populations[train_index]
    test_population = populations[test_index][0]
    print(f"Leave one out for population: {test_population}")
    print("Training populations:")
    print(train_populations)
    # Subset your dataframes for training and testing based on populations
    train_fold_snp = snp.loc[snp.index.isin(train_populations)]
    test_fold_snp = snp.loc[snp.index == test_population]
    train_fold_env = env.loc[env.index.isin(train_populations)]
    test_fold_env = env.loc[env.index == test_population]
    # Save training and testing sets to files in the loocv directory
    train_fold_snp.to_csv(f"{output_dir}/snp_train_pop_{test_population}.csv")
    test_fold_snp.to_csv(f"{output_dir}/snp_test_pop_{test_population}.csv")
    train_fold_env.to_csv(f"{output_dir}/env_train_pop_{test_population}.csv")
    test_fold_env.to_csv(f"{output_dir}/env_test_pop_{test_population}.csv")
    # For simplicity, print a separator line
    print('-' * 40)

# Save the full dataframes with the index as the population column in the loocv directory
snp.to_csv(f"{output_dir}/full_snp.csv")
env.to_csv(f"{output_dir}/full_env.csv")





###! data wrangling below (1.3)


###*#########################*###
###*  data wrangling (full) *###
###*#########################*###
#! R 

library(tidyr)
library(dplyr)
setwd("/workdir/kcarbeck/data/loocv")

df <- read.csv("full_snp.csv", header = TRUE)
subset_df <- df[, c('population', 'SNP', 'MAF')]

# Create a new data frame with 'population' as rows, 'SNP' as columns, and 'maf' as values
snpfile <- pivot_wider(subset_df, id_cols = population, names_from = SNP, values_from = MAF, values_fn = mean)
head(snpfile)

snpfile_df <- data.frame(snpfile, row.names = NULL)


# set the row names to be the values in the 'population' column
rownames(snpfile_df) <- snpfile_df$population
# remove 'population' column
snpfile_df <- snpfile_df[, -1]
str(snpfile_df)
# save df
write.table(snpfile_df, file = "snpfile_full.txt", sep = "\t", row.names = TRUE)


###*#################################*###
 ###*  data wrangling (LOOCV FILES) *###
###*#################################*###

# get a list of files matching the pattern of file structure
file_list <- list.files(pattern = "^snp_(train|test)_pop_.*\\.csv$")

for (file_name in file_list) {
    # Read the file
    df <- read.csv(file_name, header = TRUE)
    subset_df <- df[, c('population', 'SNP', 'MAF')]

    # Pivot the data
    snpfile <- pivot_wider(subset_df, id_cols = population, names_from = SNP, values_from = MAF, values_fn = mean)

    # Create a new data frame without row names
    snpfile_df <- data.frame(snpfile, row.names = NULL)

    # Set the row names to be the values in the 'population' column and remove the 'population' column
    rownames(snpfile_df) <- snpfile_df$population
    snpfile_df <- snpfile_df[, -1]

    # Extract the train/test and population name from the filename
    # This regex captures the train/test part and the population name separately
    parts <- strsplit(file_name, "_pop_")[[1]]
    train_test <- parts[1]
    population_name <- sub(".csv", "", parts[2])

    # Save the processed dataframe
    save_file_name <- paste0("snpfile_", train_test, "_", population_name, ".txt")
    write.table(snpfile_df, file = save_file_name, sep = "\t", row.names = TRUE)
}

