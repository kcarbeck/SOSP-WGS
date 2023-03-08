# gradientForest with pilot data from body size manuscript (Carbeck et al, in review)
# 30 Jan 2023
# katherine carbeck

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS")
library(dplyr)
library(sf)
library(tidyr)
library(stringr)

#gradientForest was a massive pain to install. Errors: 
#ERROR: compilation failed for package ‘extendedForest’
#Warning in install.packages : installation of package ‘extendedForest’ had non-zero exit status
#ERROR: dependency ‘extendedForest’ is not available for package ‘gradientForest’
#Warning in install.packages : installation of package ‘gradientForest’ had non-zero exit status
#R was not able to compile packages...you have to install homebrew and gcc (https://formulae.brew.sh/formula/gcc#default)
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")
#install.packages("extendedForest", repos="http://R-Forge.R-project.org")
library(gradientForest)


#### DATA ####
#load environmental data of 80 predictors (measured at each site or interpolated) 
clm<-read.csv(file="data/climateNA_output.csv")  #removed solar radiation vars because a lot of missing data
clm<-clm %>%
  rename(c(sample_id=ID1, 
           subspecies=ID2)) %>%
  dplyr::select(!c(subspecies, Latitude, Longitude, Elevation))
head(clm);dim(clm) 
# 354  81

#        sample_id   Tmax_wt Tmax_sp Tmax_sm Tmax_at Tmin_wt Tmin_sp Tmin_sm Tmin_at Tave_wt Tave_sp Tave_sm Tave_at PPT_wt PPT_sp PPT_sm PPT_at DD_0_wt
# 1 AK_caurina_S11     0.3     6.8    15.5     7.8    -7.4    -1.2     7.6     0.5    -3.6     2.8    11.5     4.1    676    557    585   1001     411
# 2 AK_caurina_S12     0.3     6.8    15.5     7.8    -7.4    -1.2     7.6     0.5    -3.6     2.8    11.5     4.1    676    557    585   1001     411
# 3 AK_caurina_S13     0.5     7.0    15.5     8.1    -6.8    -0.7     8.0     1.1    -3.1     3.2    11.8     4.6    663    541    562    975     377
# 4 AK_caurina_S14     0.5     7.0    15.5     8.1    -6.8    -0.7     8.0     1.1    -3.1     3.2    11.8     4.6    663    541    562    975     377
# 5 AK_caurina_S15   -12.8     5.5    19.2     1.1   -23.5    -8.5     5.5    -9.6   -18.2    -1.5    12.3    -4.2     58     49    149     87    1643
# 6 AK_caurina_S16     0.3     6.8    15.5     7.8    -7.4    -1.2     7.6     0.5    -3.6     2.8    11.5     4.1    676    557    585   1001     411

clm2 <- clm[,-1]
rownames(clm2) <- clm[,1]
clm <-clm2



### --- ###            run below on terminal            ### --- ###
setwd("/workdir/kcarbeck")

# load subsset data of 352 individuals from 457 SNPs 
# each species is treated as dependent var in usual use of GF, but in my dataset each SNPs or MAF will be used
library(vcfR)
my_vcf <- read.vcfR("imputed_SOSP_352Samples_012923.vcf.gz")
# Scanning file to determine attributes.
# File attributes:
#   meta lines: 3007
# header_line: 3008
# variant count: 13,663,084
# column count: 361
# Meta line 3000 read in.




getFIX(my_vcf)
my_vcf@gt
?vcfR2tidy





snp<-read.csv(file="data/pilot_snpsforR.csv")
snp<-snp %>%
  select(!c(population, mass))
dim(snp)
# 79 458 (458 including sample_id)

snp2 <- snp[,-1]
rownames(snp2) <- snp[,1]
snp <-snp2
head(snp)


#### GRADIENT FOREST ANALYSIS ####

# good resource: https://github.com/pgugger/LandscapeGenomics/blob/master/2018_China/Exercise3.md

#can account for correlated predictors by implementing conditional permutation [Ellis et~al., 2010], following the strategy outlined by Strobl et~al. [2008] - the predictor to be assessed is permuted only within blocks of the dataset defined by splits in the given tree on any other predictors correlated above a certain threshold (e.g. r = 0:5) and up to a maximum number of splits set by the maxLevel option (if required).
nSites <- dim(snp)[1]
nSpecs <- dim(snp)[2] 

# In GF, a maximum number of splits can be defined following the developers suggestion
lev <- floor(log2(nSites * 0.368/2))
lev

dat<-cbind(clm,snp) 
head(dat)


#paramaters:
#typically set trees=500
#whether splits should be compact into bins = prevents memory problems for large datasets
#num bins
#correlation threshold for conditional permutation
gf <- gradientForest(cbind(clm, snp),
                     predictor.vars = colnames(clm), response.vars = colnames(snp),
                     ntree = 500, transform = NULL, compact = T,
                     nbin = 201, maxLevel = lev, corr.threshold = 0.5)
# *****try making snps a factor to run classification rather than regression ******** # 
# The input is the combined climate, spatial and SNP data as input (cbind(env.gf, snp)), and the subsequent parts of the command define which variables are predictors and response variables, as well as a number of other parameters that I have left as suggested. When it finishes, there will be warnings about having less than five values for response variables, which is because we have only three: 0, 1, or 2. You can ignore them.

#50 or more warnings:
#In randomForest.default(x = X, y = spec_vec, maxLevel = maxLevel,  ... :
#The response has five or fewer unique values.  Are you sure you want to do regression?

gf #summary 
# A forest of 500 regression trees for each of 441 species
# Call:

# gradientForest(data = cbind(clm, snp), predictor.vars = colnames(clm), 
#                 response.vars = colnames(snp), ntree = 500, transform = NULL, 
#                 maxLevel = 0, corr.threshold = 0, compact = T, nbin = 201)

# Important variables:
# [1] MAP  CMI  AHM  RH   MWMT


names(gf)
#[1] "X"               "Y"               "result"          "overall.imp"     "overall.imp2"    "ntree"           "imp.rsq"         "species.pos.rsq"
#[9] "ranForest.type"  "res"             "res.u"           "dens"            "call" 


#### GRADIENT FOREST PLOTS ####

# 1. predictor overall importance plot (shows mean accuracy importance and mean importance weighted by SNP R2). In this example, both are conditional importance.
plot(gf, plot.type = "O")
most_important <- names(importance(gf))[1:25]
par(mgp = c(2, 0.75, 0))


# 1a. which vars are most important?
# We can also plot the "turnover functions" showing how allelic composition changes along the spatial or environmental gradients. The shapes are nonlinear and large jumps show steep genetic changes along certain portions of the environmental gradient. The height that the function acheives on the right side of the plot is the total importance and should match the barplot. First, organize the variables by importance and then plot:
by.importance <- names(importance(gf))

plot(gf, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))

# For individual loci: 
# species cumulative plot (shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardized by density of observations). Shows cumulative change in abundance of individuals species, where changes occur on the gradient, and the species changing most on each gradient. 
plot(gf, plot.type = "C", imp.vars = most_important,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
     cex.axis = 0.6, line.ylab = 0.9, 
     par.args = list(mgp = c(1.5,0.5, 0), mar = c(2.5, 1, 0.1, 0.5), 
                     omi = c(0, 0.3, 0, 0)))
# Each line within each panel represents allelic change at a single SNP. Notice that in each panel some SNPs show very steep changes along the environmental gradient. One might consider these SNPs as              candidates for involvement in local adaptation along the gradient. This approach to "outlier" detection      is still being tested, but Fitzpatrick & Keller (2015) show a promising example.




# 2. splits density plot (shows binned split importance and location on each gradient (spikes), kernel density of splits (black lines), of observations (red lines) and of splits standardized by observations density (blue lines)). Each distribution integrates to predictor importance. These plots show where important changes in the abundance of multiple species are occurring along the gradient; they indicate composition change rate.
plot(gf, plot.type = "S", imp.vars = most_important,
     leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                             0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))



# 4. predictor cumulative plot (shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations, averaged over all species). These show cumulative change in overall composition of the community, where changes occur on the gradient.
#common.scale=T ensures that plots for all predictors have the same y-scale as the most important predictor

plot(gf, plot.type = "C", imp.vars = most_important,
     show.species = F, common.scale = T, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                             0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))


# 5. shows the R2 measure of the fit of the random forest model for each species, ordered in various ways
plot(gf, plot.type = "P", show.names = F, horizontal = F,
     cex.axis = 1, cex.labels = 0.7, line = 2.5)

plot(gf, plot.type = "P", show.names = T, horizontal = F,
     cex.axis = 1, cex.labels = 0.7, line = 2.5)

plot(gf, plot.type = "P", show.names = F, horizontal = T,
     cex.axis = 1, cex.labels = 0.6, line = 2.5)

plot(gf, plot.type = "P", show.names = T, horizontal = T,
     cex.axis = 1, cex.labels = 0.6, line = 2.5)



##############################################################
#### subset to only run with important vars:
# MAP  CMI  AHM  RH   MWMT
head(clm)
clm2<-clm %>%
  select(MAP, CMI, AHM, RH, MWMT)
clm2

gf <- gradientForest(cbind(clm2, snp),
                     predictor.vars = colnames(clm2), response.vars = colnames(snp),
                     ntree = 500, transform = NULL, compact = T,
                     nbin = 201, maxLevel = lev, corr.threshold = 0.5)
gf

