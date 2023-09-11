# summarize data for committee meeting
# 21 july 2023
# katherine carbeck

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/GEA-data")

library(tidyverse)

dat <- read.csv("env-seasonal2.csv")
str(dat)
dat$population<- as.factor(dat$population)

# summarize number of samples by population
dat.summarized <- 
  dat %>%
  count(population) %>%
  arrange(desc(n))

write.csv(dat.summarized, "summarized_number_of_samples_populuation.csv")


# calculate the average of each column based on groups
by.pop <- dat %>%
  group_by(Population) %>%
  summarize(across(everything(), mean))