# use '04.gradient_fitting_script.R' script to fit GF model (based on Brandon Lind scripts)
# katherine carbeck
# 14 december 2023

###*############################################*###
        ###*  test GF fitting script  *###
###*############################################*###

##files from 01.splitTrainingAndTesting.py > 01a.
# maf file = WZA_Outlier_SNPs_MAFs.csv = df1 (output from 01.splitTrainingAndTesting.py)
    # data wrangled version of this file = snpfile.txt
# env file = scaledPopEnvAnnualUpdatedwLATLONG.csv = df2 (output from 01.splitTrainingAndTesting.py)
# range_file= combined_raster_data.txt (output from 00.clipClimateToRange.R)


python
from pythonimports import *
import xarray, rioxarray
import subprocess
DIR = '/workdir/kcarbeck'

# test fitting script on a very small snps file
fitting_script = '/workdir/kcarbeck/04.gradient_fitting_script.R'
gfOut_trainingfile = '/workdir/kcarbeck/data/testdir/test_gradient_forest_training.RDS'
range_file = '/workdir/kcarbeck/data/2012_2022/combined_raster_data_2012_2022.txt'
predfile = '/workdir/kcarbeck/data/testdir/test_gradient_forest_predOut.RDS'
maskfile = '/workdir/kcarbeck/data/maskdir/rangeWGS84_mask.RDS'
basename = 'test_fit'
save_dir = '/workdir/kcarbeck/data/testdir'


command = [
    'Rscript',
    fitting_script,
    gfOut_trainingfile,
    range_file,
    predfile,
    maskfile,
    basename,
    save_dir
]

subprocess.run(command, check=True)


# warning after get_offset()
#   Saved netCDF file to:
#           /workdir/kcarbeck/data/testdir/test_fit_offset.ncWarning message:
#   In x@data@values[i] <- value :
#     number of items to replace is not a multiple of replacement length






###*############################################*###
        ###*  GF fitting script  *###
###*############################################*###

##files from 01.splitTrainingAndTesting.py > 01a.
# maf file = WZA_Outlier_SNPs_MAFs.csv = df1 (output from 01.splitTrainingAndTesting.py)
    # data wrangled version of this file = snpfile.txt
# env file = scaledPopEnvAnnualUpdatedwLATLONG.csv = df2 (output from 01.splitTrainingAndTesting.py)
# range_file= combined_raster_data.txt (output from 00.clipClimateToRange.R)


python
from pythonimports import *
import xarray, rioxarray
import subprocess
DIR = '/workdir/kcarbeck'


#see imported modules:
import sys
for module_name in sys.modules:
     print(module_name)
     

# test fitting script on a very small snps file
fitting_script = '/workdir/kcarbeck/04.gradient_fitting_script.R'
gfOut_trainingfile = '/workdir/kcarbeck/data/training/west_full_gradient_forest_training.RDS'
range_file = '/workdir/kcarbeck/data/2012_2022/west_combined_raster_data_2012_2022.txt'
predfile = '/workdir/kcarbeck/data/training/west_full_gradient_forest_predOut_noNA.RDS'
maskfile = '/workdir/kcarbeck/data/maskdir/westWGS84_mask.RDS'
basename = 'west_full'
save_dir = '/workdir/kcarbeck/data/fitting'


command = [
    'Rscript',
    fitting_script,
    gfOut_trainingfile,
    range_file,
    predfile,
    maskfile,
    basename,
    save_dir
]

subprocess.run(command, check=True)




###*############################################*###
   ###*  run GF fitting script for K-folds *###
###*############################################*###

python
DIR = '/workdir/kcarbeck'
from concurrent.futures import ProcessPoolExecutor
import xarray, rioxarray
from pythonimports import *
import subprocess

#see imported modules:
import sys
for module_name in sys.modules:
     print(module_name)
     
# Function to run the R script
def run_r_script(k):
    fitting_script = '/workdir/kcarbeck/04.gradient_fitting_script.R'
    gfOut_trainingfile = f'/workdir/kcarbeck/data/training/west_kfold_{k}_gradient_forest_training.RDS'
    range_file = '/workdir/kcarbeck/data/2012_2022/west_combined_raster_data_2012_2022.txt'
    predfile = f'/workdir/kcarbeck/data/training/west_kfold_{k}_gradient_forest_predOut_noNA.RDS'
    maskfile = '/workdir/kcarbeck/data/maskdir/westWGS84_mask.RDS'
    basename = f'west_kfold_{k}'
    save_dir = '/workdir/kcarbeck/data/fitting'
    command = [
        'Rscript',
        fitting_script,
        gfOut_trainingfile,
        range_file,
        predfile,
        maskfile,
        basename,
        save_dir
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
