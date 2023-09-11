# summarize climateNA data by population
# 25 july 2023
# katherine carbeck

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/GEA-data")
library(tidyverse)
library(dplyr)

# *step 1* july 20 - create popfile.tsv in the right format for baypass
by.pop <- read.csv("pop.csv",header = FALSE)
head(by.pop)
write.table(by.pop, file='popfile.tsv', quote=FALSE, sep='\t', col.names=FALSE, row.names = FALSE)


# july 25 - summarize climate vars (climate NA output) by population for baypass input
clm<-read.csv("ClimateNA_output2.csv", stringsAsFactors = TRUE)
str(clm)
sub<-subset(clm, select=-c(ID1,long,lat,el))
str(sub)

by.pop <- sub %>%
  group_by(ID2) %>%
  summarize(across(everything(), mean))

str(sub)
write.csv(by.pop, "popenv.csv")



###   standardize vars   ###
clm<-read.csv("popenv.csv")
str(clm)

sub<-subset(clm, select=-(Population))
str(sub)

scaled <- sub %>% 
  mutate_all(~(scale(.) %>% 
                 as.vector))
str(scaled)

# cbind pop and scaled vars
pop <- subset(clm, select=(Population))
scaledpop<- cbind(pop, scaled)

write.csv(scaledpop, "scaledPopEnv.csv")


## change format for baypass input
write.table(scaledpop, file='scaledPopEnv.tsv', quote=FALSE, sep='\t', row.names = FALSE)











### old ####

dat <- read.csv("env-seasonal.csv")
str(dat)
dat$Population<- as.factor(dat$Population)

# summarize by population
# by.pop <- 
#   dat %>% group_by(Population) %>%
#   summarise(n = n(),
#            Population = mean(Population, na.rm = TRUE)
#  )

# calculate the average of each column based on groups
by.pop <- dat %>%
  group_by(Population) %>%
  summarize(across(everything(), mean))

write.csv(by.pop, "popenv.csv")

df4 <- by.pop %>% 
  mutate_all(~(scale(.) %>% as.vector))
df4

# read in data for pca
# by.pop <- read.csv("popenv.csv")

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

write.csv(by.pop, "popenvPCA.csv")

by.pop <- read.csv("popenvPCA.csv")
head(by.pop)
write.table(by.pop, file='popenvPCA.tsv', quote=FALSE, sep='\t', row.names = FALSE)






