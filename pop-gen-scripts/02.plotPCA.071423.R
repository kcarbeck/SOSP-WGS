### katherine carbeck
## 17 july 2023
## PCA on final 316 individual VCF file (LD pruned)

# /Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen-scripts/01a.PCA.md
# .gds file brought over from terminal 

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen-data/PCA")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("SNPRelate"))
install.packages("scales")
library(SNPRelate)
library(gdsfmt)
library(scales)
library(tidyverse)



snpgdsSummary("316Samples6xLD_071623_out.gds")
  # The file name: /Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen-data/PCA/316Samples6xLD_071623_out.gds 
  # The total number of samples: 316 
  # The total number of SNPs: 6059018 
  # SNP genotypes are stored in SNP-major mode (Sample X SNP).

genofile <- snpgdsOpen("316Samples6xLD_071623_out.gds")
miss <- snpgdsSampMissRate(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss
write.csv(miss, file = "output/316Samples_missingness.csv")
pca <- snpgdsPCA(gdsobj = genofile, autosome.only=FALSE)
  # Principal Component Analysis (PCA) on genotypes:
  #   Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
  # # of samples: 316
  # # of SNPs: 6,059,018
  # using 1 thread
  # # of principal components: 32
  # PCA:    the sum of all selected genotypes (0,1,2) = 2622476591
  # CPU capabilities:
  #   Mon Jul 17 13:21:24 2023    (internal increment: 1552)
  # [==================================================] 100%, completed, 4.5m
  # Mon Jul 17 13:25:57 2023    Begin (eigenvalues and eigenvectors)
  # Mon Jul 17 13:25:57 2023    Done.


pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# [1] 4.87 4.02 2.44 2.14 1.93 1.24
tab <- data.frame(sample = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],    # the second eigenvector
                  PC3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)
head(tab)
tab
write.csv(tab, file = "output/316samples0714_PCA.csv")
tab <- read.csv("output/316samples0714_PCA.csv",stringsAsFactors = F)
str(tab)
tab$Population <- as.factor(tab$Population)
tab$Subspecies <- as.factor(tab$Subspecies)
str(tab)


#########       plot       ##########

library(randomcoloR)
n <- 29
palette <- distinctColorPalette(n)

library(plotly)

# by population 
ggplot(data=tab,
       mapping=aes(x=PC1,
                   y=PC2, 
                   color=Population,
                   label=sample)) + 
  xlab("PC1 4.87%")+
  ylab("PC2 4.02%")+
  geom_point() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values=palette) 
# geom_text_repel(direction="both", size=1, max.overlaps = 50)
# geom_text(hjust=.5, vjust=.5)
ggplotly()


# color plot by subspecies 
# hex colors to match subspecies map
ggplot(data=tab,
       mapping=aes(x=PC1,
                   y=PC2, 
                   color=Subspecies)) + 
  xlab("PC1 4.87%")+
  ylab("PC2 4.02%")+
  geom_point() +
  scale_color_manual(values=c("maxima"="#FAD114",
                              "sanaka"="#F57B7A",
                              "insignis"="#E7D8EF",
                              "kenaiensis"="#FEFF7F",
                              "caurina"="#F8AA03",
                              "rufina"="#39AFB6",
                              "merrilli"="#C19ED6",
                              "morphna"="#89CD66",
                              "cleonensis"="#FACB10",
                              "montana"="#D6C29D",
                              "samuelis"="#BFE9E8",
                              "pusillula"="#61C4CA",
                              "maxillaris"="#1B4865",
                              "gouldii"="#5FA8D3",
                              "heermanni"="#CAE9FE",
                              "fallax"="#D69E9D",
                              "graminea"="#C55A11",
                              "rivularis"="#FCDD9B",
                              "adusta"="#CCCCCC",
                              "mexicana"="#FBCFBE",
                              "nominate"="#D8D79E"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

