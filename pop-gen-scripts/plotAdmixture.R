#!/usr/bin/Rscript

# plot structure from ADMIXTURE 
# 24 feb 2023
# Author: katherine carbeck
library(tidyverse)

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen/ADMIXTURE")

#samplelist <- read_delim("filtered_SOSP_352Samples_122322.list", delim=" ",
#                         col_names = c("sample", "population"))
samplelist <- read.csv("listFile.csv")

# CV from first run:
  # (K=1): 0.30681
  # (K=2): 0.44270
  # (K=3): 0.47874
  # (K=4): 0.43771
  # (K=5): 0.48460
  # (K=6): 0.39155
  # (K=7): 0.38747
  # (K=8): 0.40383
  # (K=9): 0.58485
  # (K=10): 0.58806
  # (K=11): 0.60077
  # (K=14): 0.63035
  # (K=20): 0.29449

# We're using the command _read\_delim_ which requires column names. Since the data file doesn't have any column names, we have to specify them and we're using a combination of paste and seq to produce "Q1", "Q2". We could have hard coded it c("Q1","Q2") but this way it works for an arbitrary number of columns just by changing the second value in seq(). 

# Now we could work on each value of K individually, but its easier to load all of them at once. One problem is that they each have different numbers of columns. The solution is converting from a wide to long format. In a long format, each row has a single data point instead of multiple. The tidyverse tools are set up to prefer long data ((and there are other reasons)[https://sejdemyr.github.io/r-tutorials/basics/wide-and-long/#a-case-for-long-data]) so lets do that. 

# Its possible to load multiple files in a single call, but for transparency lets use a loop. We first make an empty dataframe that we're going to fill, then we loop through our output files, convert them to long format and add them to the master set.


all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

#vec<- c(1,2,3,4,5,6,7,8,9,10,11,)
for (k in 1:8){
  data <- read_delim(paste0("filtered_SOSP_352Samples_122322.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$sample
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data


'# A tibble: 23,232 × 4
   sample             k Q     value
   <chr>          <int> <chr> <dbl>
 1 AK_caurina_S11     1 Q1        1
 2 AK_caurina_S12     1 Q1        1
 3 AK_caurina_S13     1 Q1        1
 4 AK_caurina_S14     1 Q1        1
 5 AK_caurina_S15     1 Q1        1
 6 AK_caurina_S16     1 Q1        1
 7 AK_caurina_S17     1 Q1        1
 8 AK_caurina_S18     1 Q1        1
 9 AK_caurina_S19     1 Q1        1
10 AK_caurina_S20     1 Q1        1
'

all_data <- as.data.frame(all_data)

dat <- left_join(x=all_data, y=samplelist, by = "sample", all.x = TRUE)
write.csv(dat, "admixture_data.csv")
head(dat)


# Plot
#positions <- (limits=c("maxima", "sanaka", "rufina", "merrilli"))

#sub <- subset(dat, k< 8)

tiff(file="admixture_plot.tiff",width = 10000, height = 5000,res=300)  
dat %>%
  arrange(factor(population, 
                 #levels = c("maxima", "sanaka", "rufina", "merrilli")
                 )) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  #scale_fill_manual(values = c("#C3C1C1","#A4A2A2","#777777", "#3E3E3E"), name="K",
                  #  labels=seq(1:8)) +
  facet_wrap(~k,ncol=1)
dev.off()


###############################################################################################################

### look at the k=20 now
samplelist <- read.csv("listFile.csv")
all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 20){
  data <- read_delim(paste0("filtered_SOSP_352Samples_122322.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$sample
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data

# add subspp and pop columns
all_data <- as.data.frame(all_data)

dat <- left_join(x=all_data, y=samplelist, by = "sample", all.x = TRUE)
head(dat)
#

tiff(file="admixture_k20_plot.tiff",width = 10000, height = 2000,res=300)
dat %>%
  arrange(factor(subspecies, levels = c("maxima", "sanaka", "insignis", "kenaiensis", "caurina", "rufina", "merrilli", "morphna", "cleonensis", "montana","samuelis","pusillula","maxillaris","gouldii","heermanni","fallax","graminea","rivularis","adusta","mexicana","nominate","unknown"))) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Population") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#FAD114", "#F57B7A","#E7D8EF","#FEFF7F", "#F8AA03", "#39AFB6","#C19ED6", "#89CD66","#FACB10", "#D6C29D",
                              "#BFE9E8", "#61C4CA", "#1B4865", "#5FA8D3", "#CAE9FE", "#D69E9D", "#C55A11", "#FCDD9B", "#CCCCCC", "#FBCFBE",
                              "#D8D79E", "#404040"), 
                    name="K",
                    labels=seq(1:20))  +
  facet_wrap(~k,ncol=1)
dev.off()










