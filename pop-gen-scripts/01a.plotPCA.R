### katherine carbeck
## 14 march 2023
## first PCA on final 319 individual VCF file (no LD pruning yet)

# /Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen-scripts/01a.PCA.md
# .gds file brought over from terminal 

setwd("/workdir/kcarbeck/pca")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("SNPRelate"))
install.packages("scales")
library(RColorBrewer)
library(SNPRelate)
library(gdsfmt)
library(scales)
library(tidyverse)

