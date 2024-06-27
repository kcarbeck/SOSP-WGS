# use '04.gradient_fitting_script.R' script to fit GF model (based on Brandon Lind scripts)
# katherine carbeck
# 18 jan 2024

###*############################################*###
   ###*  run GF fitting script for LOOCV *###
###*############################################*###

python
DIR = '/workdir/kcarbeck'
from concurrent.futures import ProcessPoolExecutor
import xarray, rioxarray
from pythonimports import *
import subprocess
import os
import re

#see imported modules:
import sys
for module_name in sys.modules:
     print(module_name)


def get_population_names(directory, pattern):
    # get a list of all files in the directory
    all_files = os.listdir(directory)
    # filter files based on the pattern
    filtered_files = [f for f in all_files if re.match(pattern, f)]
    # extract the population names
    population_names = [re.sub(pattern, r'\1', f) for f in filtered_files]
    # remove duplicates and return
    return list(set(population_names))


# directory where the files are stored
directory = '/workdir/kcarbeck/data/loocv/'

# regex pattern to match filenames and get population names (assumes names start with 'snpfile_snp_train_' followed by the pop name and '.txt')
pattern = r'snpfile_snp_train_(.+)\.txt'

# get the list of population names
population_names = get_population_names(directory, pattern)
print(population_names) 

# function to run GF fitting R script
def run_r_script(population_name):
    fitting_script = '/workdir/kcarbeck/04.gradient_fitting_script.R'
    gfOut_trainingfile = f'/workdir/kcarbeck/data/loocv/training/west_loocv_{population_name}_gradient_forest_training.RDS'
    range_file = '/workdir/kcarbeck/data/2012_2022/west_combined_raster_data_2012_2022.txt'
    predfile = f'/workdir/kcarbeck/data/loocv/training/west_loocv_{population_name}_gradient_forest_predOut_noNA.RDS'
    maskfile = '/workdir/kcarbeck/data/maskdir/westWGS84_mask.RDS'
    basename = f'west_loocv_{population_name}'
    save_dir = '/workdir/kcarbeck/data/loocv/fitting'
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


# use ProcessPoolExecutor to run subprocesses in parallel
#max_workers = os.cpu_count()
max_workers = 5
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    # submit a job for each population name
    futures = [executor.submit(run_r_script, name) for name in population_names]
    # wait for all jobs to complete
    for future in futures:
        future.result()

