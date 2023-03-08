# summarize climateNA data by population
# 28 feb 2023
# katherine carbeck

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/GEA")
library(tidyverse)

dat <- read.csv("input/env.csv")
str(dat)
dat$POP<- as.factor(dat$POP)

# summarize by population
by.pop <- 
  dat %>% group_by(POP) %>%
  summarise(n = n(),
            POP = mean(POP, na.rm = TRUE)
  )

# calculate the average of each column based on groups
by.pop <- dat %>%
  group_by(POP) %>%
  summarize(across(everything(), mean))

write.csv(by.pop, "input/popenv.csv")

# read in data for pca
by.pop <- read.csv("input/popenv.csv")

# Subset the data to include only the environmental variables
env_vars <- by.pop[,2:ncol(by.pop)]

# Scale the environmental variables
env_vars_scaled <- scale(env_vars)

# Perform PCA
pca_result <- prcomp(env_vars_scaled)

# Print the summary of the PCA
summary(pca_result)

# Plot the scree plot
plot(pca_result, type="l")

# Plot the biplot
biplot(pca_result)

# Project the data onto the principal components
pca_data <- predict(pca_result, newdata=env_vars_scaled)

# Add the PCA results back into the original data frame
by.pop[,paste0("PC",1:ncol(pca_data))] <- pca_data

# View the updated data frame
head(by.pop)

write.csv(by.pop, "input/popenvPCA.csv")

by.pop <- read.csv("input/popenvPCA.csv")
head(by.pop)
write.table(by.pop, file='input/popenvPCA.tsv', quote=FALSE, sep='\t', row.names = FALSE)
