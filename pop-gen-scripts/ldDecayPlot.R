# 23 Feb 2023
# author: katherine carbeck
# LD decay plot based on script from https://speciationgenomics.github.io/ld_decay/

library(tidyverse)
setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen/ldDecay")


# set path
my_bins <- "./ldDecayOut.ld_decay_bins"

# read in data
ld_bins <- read_tsv(my_bins)

# plot LD decay
ggplot(ld_bins, aes(distance, avg_R2)) +
  geom_line() + 
  xlab("Distance (bp)") + ylab(expression(italic(r)^2))


# ld block density
ggplot(ld_bins,aes(distance))+
  geom_density(size=0.5,colour="grey40")+
  labs(x="LD block length (bp)",y="Density")+
  theme_bw()

# ggsave("snp-thin-ld-blocks.png",p,height=8,width=8,units="cm",dpi=250)
