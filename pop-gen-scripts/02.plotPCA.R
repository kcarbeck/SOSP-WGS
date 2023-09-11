### katherine carbeck
## 27 march 2023
## PCA on final 316 individual VCF file (*not* LD pruned)

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



snpgdsSummary("316Samples_032123_out.gds")
  # The file name: /Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen-data/PCA/316Samples_032123_out.gds 
  # The total number of samples: 316 
  # The total number of SNPs: 12483733 
  # SNP genotypes are stored in SNP-major mode (Sample X SNP).

genofile <- snpgdsOpen("316Samples_032123_out.gds")
miss <- snpgdsSampMissRate(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
miss
write.csv(miss, file = "output/316Samples_missingness.csv")
pca <- snpgdsPCA(gdsobj = genofile, autosome.only=FALSE)
  # Principal Component Analysis (PCA) on genotypes:
  #   Excluding 2,850,400 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
  # # of samples: 316
  # # of SNPs: 9,633,333
  # using 1 thread
  # # of principal components: 32
  # PCA:    the sum of all selected genotypes (0,1,2) = 4215170425
  # CPU capabilities:
  #   Mon Mar 27 15:32:32 2023    (internal increment: 1552)
  # [==================================================] 100%, completed, 7.2m 
  # Mon Mar 27 15:39:44 2023    Begin (eigenvalues and eigenvectors)
  # Mon Mar 27 15:39:44 2023    Done.


pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# [1] 6.18 4.57 2.67 2.22 1.95 1.19
tab <- data.frame(sample = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],    # the second eigenvector
                  PC3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)
head(tab)
tab
write.csv(tab, file = "output/316samples_PCA.csv")
tab <- read.csv("output/316samples_PCA.csv",stringsAsFactors = T)
str(tab)



# color plot by subspecies 
# hex colors to match subspecies map
ggplot(data=tab,
       mapping=aes(x=PC1,
                   y=PC2, 
                   color=subspecies)) + 
  xlab("PC1 6.18%")+
  ylab("PC2 4.57%")+
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
                              "nominate"="#D8D79E",
                              "unknown"="#404040"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




######################################
########## LD pruned PCA ############
snpgdsSummary("316Samples_LDpruned_out.gds")
  # The file name: /Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen-data/PCA/316Samples_LDpruned_out.gds 
  # The total number of samples: 316 
  # The total number of SNPs: 6841816 
  # SNP genotypes are stored in SNP-major mode (Sample X SNP).

genofile <- snpgdsOpen("316Samples_LDpruned_out.gds")
pca <- snpgdsPCA(gdsobj = genofile, autosome.only=FALSE)
  #Principal Component Analysis (PCA) on genotypes:
  #  Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
  ## of samples: 316
  ## of SNPs: 6,841,816
  #using 1 thread
  ## of principal components: 32
  #PCA:    the sum of all selected genotypes (0,1,2) = 2958790577
  #CPU capabilities:
  #  Mon Mar 27 15:57:55 2023    (internal increment: 1552)
  #[==================================================] 100%, completed, 5.2m
  #Mon Mar 27 16:03:06 2023    Begin (eigenvalues and eigenvectors)
  # Mon Mar 27 16:03:06 2023    Done.


pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# [1] 4.77 3.92 2.42 2.08 1.92 1.21
tab <- data.frame(sample = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],    # the second eigenvector
                  PC3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)
head(tab)
tab
write.csv(tab, file = "output/316samples_PCA_LDpruned.csv")
tab <- read.csv("output/316samples_PCA_LDpruned.csv",stringsAsFactors = T)
str(tab)



# color plot by subspecies 
# hex colors to match subspecies map
ggplot(data=tab,
       mapping=aes(x=PC1,
                   y=PC2, 
                   color=subspecies)) + 
  xlab("PC1 4.77%")+
  ylab("PC2 3.92%")+
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
                              "nominate"="#D8D79E",
                              "unknown"="#404040"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


