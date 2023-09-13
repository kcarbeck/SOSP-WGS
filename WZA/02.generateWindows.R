## katherine carbeck
## 13 sept 2023
## script to add a column to the correlation file containing information on analysis windows 


# using output file from 01.prepForWZA.sh + computeCorrelationsForWZA.py = output_corr.correlations.csv

# local wd
setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/WZA")
# server wd
setwd("/workdir/kcarbeck/out")

library(tidyverse)
library(readr)

# format required:
# CHROM,POS,window,
# 1,14,1:1-1000, ...
# 1,567,1:1-1000, ...


corr <- read.csv("output_corr.correlations.csv")
#if gzipped:
corr <- read.csv(gzfile("output_corr.correlations.csv.gz"))

windows <- corr %>%
  mutate(pos = as.numeric(pos),  # convert "pos" to numeric
         bin_start = 1 + floor((pos - 1) / 1000) * 1000,  # calc bin start
         bin_end = bin_start + 999,  # calc bin end
         window = paste(chrom, ":", bin_start, "-", bin_end, sep = ""))  # create the "window" column

# save new df to a new csv
write.csv(windows, "output_corr_windows.correlations.csv", row.names = FALSE)

