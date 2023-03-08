# 23 Feb 2023
# author: katherine carbeck
# linkage pruned PCA with 352 SOSP WGS samples based on script from https://speciationgenomics.github.io/pca/

library(tidyverse)
setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen/PCA")
pca <- read_table("./ldPrunedVCF.eigenvec", col_names = FALSE)
eigenval <- scan("./ldPrunedVCF.eigenval")
head(pca)

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
write.csv(pca, "eigenvec.csv", row.names = FALSE)
pca<-read.csv("eigenvec.csv")


# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
pve$pve
# 23.904044
# 17.766158

# color plot by subspecies 
# hex colors to match subspecies map
ggplot(data=pca,
       mapping=aes(x=PC1,
                   y=PC2, 
                   color=subspecies)) + 
  xlab("PC1 23.90%")+
  ylab("PC2 17.77%")+
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




