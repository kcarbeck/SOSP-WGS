# wrangle data for input into GF
# katherine carbeck
# 6 december 2023

#! Bash
# copy files from home dir to workdir
cp -R /home/lc736_0001/song_sparrow/final_vcf/GF/datasets /workdir/kcarbeck/data &
mv datasets/* .

# copy raster files from home to workdir
cp -R /home/lc736_0001/song_sparrow/final_vcf/GF/climNAdat/Normal_1991_2020Y /workdir/kcarbeck/data

# copy R imports to workdir
cp /home/lc736_0001/song_sparrow/final_vcf/GF/imports.R /workdir/kcarbeck/r_imports

#!


###*#########################*###
###*  data wrangling (full) *###
###*#########################*###
#! R 

library(tidyr)
setwd("/workdir/kcarbeck/data")

df <- read.csv("full_df1.txt", header = TRUE)
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
write.table(snpfile_df, file = "snpdir/snpfile_full.txt", sep = "\t", row.names = TRUE)


###*############################*###
###*  data wrangling (k-folds) *###
###*############################*###

# get a list of files matching the pattern of file structure
file_list <- list.files(pattern = "df1_(training|testing)-k\\d+\\.txt")

# loop through each file in file_list
for (file in file_list) {
  # read in data
  df <- read.csv(file, header = TRUE)
  
  # extract info from the file name
  file_info <- strsplit(file, "_|-|\\.")[[1]]
  population_type <- file_info[2]
  k_value <- file_info[3]
  
  # subset the data
  subset_df <- df[, c('population', 'SNP', 'MAF')]
  
  # create 'population' as rows, 'SNP' as columns, and 'maf' as values
  snpfile <- pivot_wider(subset_df, id_cols = population, names_from = SNP, values_from = MAF, values_fn = mean)
  
  # create file name for output
  output_file <- paste("snpdir/" , "snpfile", "_", population_type, "_", k_value, ".txt", sep = "")
  
  # save df
  write.table(snpfile, file = output_file, sep = "\t", row.names = TRUE)
  
  # print a message for completion of each file
  cat("Processed", file, "\n")
}



###*#######################*###
   ###*  create test df *###
###*#######################*###
# Set the seed for reproducibility
set.seed(123)

# Randomly select 500 SNP column names
selected_columns <- sample(setdiff(colnames(snpfile_df), "population"), 500)

# Include the 'population' column in the selected columns
selected_columns <- c("population", selected_columns)

# Subset the data frame with the selected columns
random_columns_subset <- snpfile_df[selected_columns]
str(random_columns_subset)

# Set the row names to be the values in the 'population' column
rownames(random_columns_subset) <- random_columns_subset$population
# Remove the 'population' column
random_columns_subset <- random_columns_subset[, -1]
# Display the structure of the modified data frame
str(random_columns_subset)

write.table(random_columns_subset, file = "testdir/test_snps_delete.txt", sep = "\t", row.names = TRUE)




###*###############################*###
  ###*  rename and move env data *###
###*###############################*###
#! bash
# Move and rename testing files
for file in df2_testing-k*.txt; do
  k_value=$(echo "$file" | sed -n 's/.*k\([0-9]\+\)\.txt/\1/p')
  mv "$file" "envdir/envfile_testing_k$k_value.txt"
done

# Move and rename training files
for file in df2_training-k*.txt; do
  k_value=$(echo "$file" | sed -n 's/.*k\([0-9]\+\)\.txt/\1/p')
  mv "$file" "envdir/envfile_training_k$k_value.txt"
done


mv full_df2.txt envdir/envfile_full.txt
#!