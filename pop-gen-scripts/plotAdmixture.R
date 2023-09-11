#!/usr/bin/Rscript

# plot structure from ADMIXTURE 
# 27 july 2023
# Author: katherine carbeck
library(tidyverse)

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/pop-gen-data/ADMIXTURE")


# use sample/pop list prepared for baypass 
samplelist <- read.table(file='popfile.tsv', sep='\t', quote='',
                 col.names = c('sample', 'population') )


## old:
# samplelist <- read_delim("LDpruned_plink.list", delim=" ",
#                         col_names = c("sample", "population"))

## CV errors:
  # 1 0.57670
  # 2 0.54944
  # 3 0.52711
  # 4 0.51932
  # 5 0.51369
  # 6 0.50779
  # 7 0.50794
  # 8 0.50542
  # 9 0.50409 **
  # 10 0.50549
  # 11 0.50957
  # 12 0.51479
  # 13 0.51728
  # 14 0.51907
  # 15 0.52249
  # 16 0.52772
  # 17 0.53622
  # 18 0.54955
  # 19 0.54804
  # 20 0.55295
  # 21 0.56145
  # 22 0.58397

# We're using the command _read\_delim_ which requires column names. Since the data file doesn't have any column names, we have to specify them and we're using a combination of paste and seq to produce "Q1", "Q2". We could have hard coded it c("Q1","Q2") but this way it works for an arbitrary number of columns just by changing the second value in seq(). 

# Now we could work on each value of K individually, but its easier to load all of them at once. One problem is that they each have different numbers of columns. The solution is converting from a wide to long format. In a long format, each row has a single data point instead of multiple. The tidyverse tools are set up to prefer long data ((and there are other reasons)[https://sejdemyr.github.io/r-tutorials/basics/wide-and-long/#a-case-for-long-data]) so lets do that. 

# Its possible to load multiple files in a single call, but for transparency lets use a loop. We first make an empty dataframe that we're going to fill, then we loop through our output files, convert them to long format and add them to the master set.


all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:22){
  data <- read_delim(paste0("LDpruned_plink.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$sample
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data

# A tibble: 79,948 × 4
#    sample                k Q     value
#    <chr>             <int> <chr> <dbl>
#  1 "ON_melodia_81  "     1 Q1        1
#  2 "ON_melodia_83  "     1 Q1        1
#  3 "ON_melodia_88  "     1 Q1        1
#  4 "AK_caurina_S11 "     1 Q1        1
#  5 "AK_caurina_S12 "     1 Q1        1
#  6 "AK_caurina_S13 "     1 Q1        1
#  7 "AK_caurina_S14 "     1 Q1        1
#  8 "AK_caurina_S15 "     1 Q1        1
#  9 "AK_caurina_S16 "     1 Q1        1
# 10 "AK_caurina_S17 "     1 Q1        1


dat <- samplelist %>% 
  left_join(x=all_data, 
            y=samplelist, 
            by = c("sample"="sample"))

write.csv(dat, "admixture_data.csv")
head(dat)


## GGPLOTLY - plot by population w/ all ks 
library(randomcoloR)
n <- 29
palette <- distinctColorPalette(n)

library(plotly)
dat %>%
  arrange(factor(population)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Population") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = palette, 
                    name="K",
                    labels=seq(1:22))  +
  facet_wrap(~k,ncol=1)
ggplotly()


##### 


## REG GGPLOT - plot by population w/ all ks 
jpeg(file="output/316Samples6x_k1-22.jpeg",width = 3000, height = 4000)
dat %>%
  arrange(factor(population)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Population") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = palette, 
                    name="K",
                    labels=seq(1:22))  +
  facet_wrap(~k,ncol=1)
dev.off()



## REG GGPLOT - plot by population w/ ONLY *k=9* most supported

#filter
k09 <- dat %>%
  dplyr::filter(k == 9)
#plot
k09 %>%
  arrange(factor(population)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Population") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = palette, 
                    name="K",
                    labels=seq(1:22))  +
  facet_wrap(~k,ncol=1)


## REG GGPLOT - plot by population w/ ONLY *k=22* just for fun

#filter
k22 <- dat %>%
  dplyr::filter(k == 22)
#plot
k22 %>%
  arrange(factor(population)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Population") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = palette, 
                    name="K",
                    labels=seq(1:22))  +
  facet_wrap(~k,ncol=1)












###############################################################################################################


# old



# Plot
#positions <- (limits=c("maxima", "sanaka", "rufina", "merrilli"))

pdf(file="admixture_plot.pdf",width = 10000, height = 5000,res=300)  
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







###########################



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

pdf(file="admixture_k20_plot.pdf",width = 10000, height = 2000,res=300)
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










