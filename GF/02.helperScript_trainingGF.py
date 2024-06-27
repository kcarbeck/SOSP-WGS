# use '02.gradient_training.R' script to fit GF model (based on Brandon Lind scripts)
# katherine carbeck
# 6 december 2023

##!

##! Bash
git clone https://github.com/brandonlind/pythonimports.git
##!


###*###################################################*###
  ###* test GF training script w small set of snps *###
###*###################################################*###
##! python

##files from 01.splitTrainingAndTesting.py & 01a.dataWrangling.R
# maf file = WZA_Outlier_SNPs_MAFs.csv = df1 (output from 01.splitTrainingAndTesting.py)
    # data wrangled version of this file = snpfile_full.txt
# env file = scaledPopEnvAnnualUpdatedwLATLONG.csv = df2 (output from 01.splitTrainingAndTesting.py)
    # moved to envdir and renamed in dataWrangling 
# range_file= combined_raster_data.txt (output from 00.clipClimateToRange.R)


python
DIR = '/workdir/kcarbeck'
from pythonimports import *
import subprocess

#install.packages("hash")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

script_file = '/workdir/kcarbeck/02.gradient_training.R'
snpfile = '/workdir/kcarbeck/data/testdir/test_snps_delete.txt'
envfile = '/workdir/kcarbeck/data/envdir/envfile_full.txt'
range_file = '/workdir/kcarbeck/data/2001_2011/combined_raster_data_2001_2011.txt'
basename = 'test'
save_dir = '/workdir/kcarbeck/data/testdir'
imports_path = '/workdir/kcarbeck/r_imports'

command = [
    'Rscript',
    script_file,
    snpfile,
    envfile,
    range_file,
    basename,
    save_dir,
    imports_path
]

subprocess.run(command, check=True)



###*############################################*###
        ###*  run GF training script  *###
###*############################################*###
##! python
python
DIR = '/workdir/kcarbeck'
from pythonimports import *
import subprocess

#see imported modules:
import sys
for module_name in sys.modules:
     print(module_name)

script_file = '/workdir/kcarbeck/02.gradient_training.R'
snpfile = '/workdir/kcarbeck/data/snpdir/snpfile_full.txt'
envfile = '/workdir/kcarbeck/data/envdir/envfile_full.txt'
range_file = '/workdir/kcarbeck/data/2001_2011/combined_raster_data_2001_2011.txt'
basename = 'full'
save_dir = '/workdir/kcarbeck/data/training'
imports_path = '/workdir/kcarbeck/r_imports'

command = [
    'Rscript',
    script_file,
    snpfile,
    envfile,
    range_file,
    basename,
    save_dir,
    imports_path
]

subprocess.run(command, check=True)

# Warning message:
# In randomForest.default(x = X, y = spec_vec, maxLevel = # maxLevel,  :
#   The response has five or fewer unique values.  Are you sure you # want to do regression?
#                   user                 system                # elapsed
# 3556.51099999999996726  754.31799999999998363 4320.# 91800000000057480
# 
# 
# Interpolating gradient forests model ...
#                   user                 system                # elapsed
# 44.4629999999997380655  1.8730000000000472937 46.# 4369999999998981366
# 
# 
# Saving files ...
# [1] "/workdir/kcarbeck/data/training/# full_gradient_forest_training.RDS"
# [1] "/workdir/kcarbeck/data/training/full_gradient_forest_predOut.# RDS"
# 
# DONE!
# [1] "December 16 2023 19:08:35"
# CompletedProcess(args=['Rscript', '/workdir/kcarbeck/02.# gradient_training.R', '/workdir/kcarbeck/data/snpdir/snpfile_full.# txt', '/workdir/kcarbeck/data/envdir/envfile_full.txt', '/workdir/# kcarbeck/data/2001_2011/combined_raster_data_2001_2011.txt', # 'full', '/workdir/kcarbeck/data/training', '/workdir/kcarbeck/# r_imports'], returncode=0)




###*############################################*###
     ###*  WEST run GF training script  *###
###*############################################*###
##! python
python
DIR = '/workdir/kcarbeck'
from pythonimports import *
import subprocess

#see imported modules:
import sys
for module_name in sys.modules:
     print(module_name)

script_file = '/workdir/kcarbeck/02.gradient_training.R'
snpfile = '/workdir/kcarbeck/data/snpdir/snpfile_full.txt'
envfile = '/workdir/kcarbeck/data/envdir/envfile_full.txt'
range_file = '/workdir/kcarbeck/data/2001_2011/west_combined_raster_data_2001_2011.txt'
basename = 'west_full'
save_dir = '/workdir/kcarbeck/data/training'
imports_path = '/workdir/kcarbeck/r_imports'

command = [
    'Rscript',
    script_file,
    snpfile,
    envfile,
    range_file,
    basename,
    save_dir,
    imports_path
]

subprocess.run(command, check=True)


'''
Warning message:
In randomForest.default(x = X, y = spec_vec, maxLevel = maxLevel,  :
  The response has five or fewer unique values.  Are you sure you want to do regression?
                  user                 system                elapsed
3736.12699999999995271  939.53100000000006276 4688.72000000000025466


Interpolating gradient forests model ...
                  user                 system                elapsed
45.6030000000000654836  7.0249999999999772626 52.8539999999993597157


Saving files ...
[1] "/workdir/kcarbeck/data/training/west_full_gradient_forest_training.RDS"
[1] "/workdir/kcarbeck/data/training/west_full_gradient_forest_predOut.RDS"

DONE!
[1] "December 20 2023 03:17:13"
CompletedProcess(args=['Rscript', '/workdir/kcarbeck/02.gradient_training.R', '/workdir/kcarbeck/data/snpdir/snpfile_full.txt', '/workdir/kcarbeck/data/envdir/envfile_full.txt', '/workdir/kcarbeck/data/2001_2011/west_combined_raster_data_2001_2011.txt', 'west_full', '/workdir/kcarbeck/data/training', '/workdir/kcarbeck/r_imports'], returncode=0)
'''



###*############################################*###








#!R
setwd("/workdir/kcarbeck/data/training")
library(gradientForest)
library(dplyr)
library(sf)
library(tidyr)
library(stringr)

trainOut <- readRDS("full_gradient_forest_training.RDS")


# 1. predictor overall importance plot (shows mean accuracy importance and mean importance weighted by SNP R2). In this example, both are conditional importance.
png("importance.png", width = 800, height = 600, units = "px", res = 100)
plot(trainOut, plot.type = "O")
most_important <- names(importance(trainOut))[1:8]
par(mgp = c(2, 0.75, 0))
dev.off()
#! imgcat /workdir/kcarbeck/data/training/importance.png


# 1a. which vars are most important?
# We can also plot the "turnover functions" showing how allelic composition changes along the spatial or environmental gradients. The shapes are nonlinear and large jumps show steep genetic changes along certain portions of the environmental gradient. The height that the function acheives on the right side of the plot is the total importance and should match the barplot. First, organize the variables by importance and then plot:
by.importance <- names(importance(trainOut))

png("turnover.png", width = 800, height = 600, units = "px", res = 100)
plot(trainOut, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
dev.off()
#! imgcat /workdir/kcarbeck/data/training/turnover.png


# For individual loci: 
# species cumulative plot (shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardized by density of observations). Shows cumulative change in abundance of individuals species, where changes occur on the gradient, and the species changing most on each gradient. 
# png("turnover_indiv.png", width = 800, height = 600, units = "px", res = 100)
# plot(trainOut, plot.type = "C", imp.vars = most_important,
#      show.overall = F, legend = T, leg.posn = "topleft",
#      leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
#      cex.axis = 0.6, line.ylab = 0.9, 
#      par.args = list(mgp = c(1.5,0.5, 0), mar = c(2.5, 1, 0.1, 0.5), 
#                      omi = c(0, 0.3, 0, 0)))
# dev.off()
##DID NOT PLOT TAKES FOREVER
#! imgcat /workdir/kcarbeck/data/training/turnover_indiv.png
# Each line within each panel represents allelic change at a single SNP. Notice that in each panel some SNPs show very steep changes along the environmental gradient. One might consider these SNPs as              candidates for involvement in local adaptation along the gradient. This approach to "outlier" detection      is still being tested, but Fitzpatrick & Keller (2015) show a promising example.



# 2. splits density plot (shows binned split importance and location on each gradient (spikes), kernel density of splits (black lines), of observations (red lines) and of splits standardized by observations density (blue lines)). Each distribution integrates to predictor importance. These plots show where important changes in the abundance of multiple species are occurring along the gradient; they indicate composition change rate.
png("splits_density.png", width = 800, height = 600, units = "px", res = 100)
plot(trainOut, plot.type = "S", imp.vars = most_important,
     leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                             0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
dev.off()
#! imgcat /workdir/kcarbeck/data/training/splits_density.png


# 4. predictor cumulative plot (shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations, averaged over all species). These show cumulative change in overall composition of the community, where changes occur on the gradient.
#common.scale=T ensures that plots for all predictors have the same y-scale as the most important predictor
png("pred_cumulative.png", width = 800, height = 600, units = "px", res = 100)
plot(trainOut, plot.type = "C", imp.vars = most_important,
     show.species = F, common.scale = T, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                             0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))
dev.off()
#! imgcat /workdir/kcarbeck/data/training/pred_cumulative.png



# 5. shows the R2 measure of the fit of the random forest model for each species, ordered in various ways
png("r2_fit.png", width = 800, height = 600, units = "px", res = 100)
par(mar = c(5,8,4,2) + 0.1)
plot(trainOut, plot.type = "P", show.names = T, horizontal = T,
     cex.axis = 1, cex.labels = 0.6, line = 2.5)
dev.off()
#! imgcat /workdir/kcarbeck/data/training/r2_fit.png







###*############################################*###
   ###*  run GF training script for K-folds *###
###*############################################*###

python
DIR = '/workdir/kcarbeck'
from concurrent.futures import ProcessPoolExecutor
from pythonimports import *
import subprocess

#see imported modules:
import sys
for module_name in sys.modules:
     print(module_name)
     
# Function to run the R script
def run_r_script(k):
    script_file = '/workdir/kcarbeck/02.gradient_training_kfold.R'
    snpfile = f'/workdir/kcarbeck/data/snpdir/snpfile_training_k{k}.txt'
    envfile = f'/workdir/kcarbeck/data/envdir/envfile_training_k{k}.txt'
    range_file = '/workdir/kcarbeck/data/2001_2011/west_combined_raster_data_2001_2011.txt'
    basename = f'west_kfold_{k}'
    save_dir = '/workdir/kcarbeck/data/training'
    imports_path = '/workdir/kcarbeck/r_imports'
    command = [
        'Rscript',
        script_file,
        snpfile,
        envfile,
        range_file,
        basename,
        save_dir,
        imports_path
    ]
    subprocess.run(command, check=True)

# Number of folds
num_folds = 5

# Using ProcessPoolExecutor to run subprocesses in parallel
with ProcessPoolExecutor() as executor:
    # Submit a job for each fold
    futures = [executor.submit(run_r_script, k) for k in range(1, num_folds + 1)]
    # Wait for all jobs to complete
    for future in futures:
        future.result()


