# use '02.gradient_training.R' script to fit GF model (based on Brandon Lind scripts)
# katherine carbeck
# 17 jan 2024


###*############################################*###
   ###*  run GF training script for LOOCV *###
###*############################################*###

python
DIR = '/workdir/kcarbeck'
from concurrent.futures import ProcessPoolExecutor
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


# function to run the GF training R script
def run_r_script(population_name):
    script_file = '/workdir/kcarbeck/02.gradient_training.R'
    snpfile = f'/workdir/kcarbeck/data/loocv/snpfile_snp_train_{population_name}.txt'
    envfile = f'/workdir/kcarbeck/data/loocv/env_train_pop_{population_name}.csv'
    range_file = '/workdir/kcarbeck/data/2001_2011/west_combined_raster_data_2001_2011.txt'
    basename = f'west_loocv_{population_name}'
    save_dir = '/workdir/kcarbeck/data/loocv/training'
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


# use ProcessPoolExecutor to run subprocesses in parallel
#max_workers = os.cpu_count()
max_workers = 5
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    # submit a job for each population name
    futures = [executor.submit(run_r_script, name) for name in population_names]
    # wait for all jobs to complete
    for future in futures:
        future.result()




###*############################################*###
    ###*  check to see everything finished *###
###*############################################*###

import os

# dir containing output files
output_dir = '/workdir/kcarbeck/data/loocv/training'

# list of pop names
population_names = [
    'montana_OR', 'caurina_AK', 'sanaka_AK', 'pusillula_CA', 'merrilli_AK', 
    'montana_N_CA', 'gouldii_CA', 'merrilli_WA', 'montana_NV', 'cleonensis_CA', 
    'maxima_AK', 'rivularis_MX', 'rufina_BC', 'insignis_AK', 'kenaiensis_AK', 
    'morphna_BC', 'fallax_CA', 'nominate_VA', 'adusta_MX', 'heermanni_N_CA', 
    'nominate_ON', 'samuelis_CA', 'heermanni_S_CA', 'mexicana_MX', 
    'fallax_UT', 'montana_S_CA', 'graminea_CA', 'maxillaris_CA', 'fallax_AZ'
]

# expected file patterns
file_patterns = [
    'west_loocv_{}_gradient_forest_predOut_noNA.RDS',
    'west_loocv_{}_gradient_forest_predOut.RDS',
    'west_loocv_{}_gradient_forest_training.RDS'
]

# generate a list of all expected files
expected_files = [pattern.format(pop) for pop in population_names for pattern in file_patterns]

# get a list of actual files in the directory
actual_files = os.listdir(output_dir)

# identify missing files
missing_files = [file for file in expected_files if file not in actual_files]

# print missing files
print("Missing Files:")
for file in missing_files:
    print(file)


# west_loocv_merrilli_WA_gradient_forest_predOut_noNA.RDS
# west_loocv_merrilli_WA_gradient_forest_predOut.RDS
# west_loocv_merrilli_WA_gradient_forest_training.RDS
# west_loocv_insignis_AK_gradient_forest_predOut_noNA.RDS
# west_loocv_insignis_AK_gradient_forest_predOut.RDS
# west_loocv_insignis_AK_gradient_forest_training.RDS
# west_loocv_morphna_BC_gradient_forest_predOut_noNA.RDS
# west_loocv_morphna_BC_gradient_forest_predOut.RDS



###*############################################*###
        ###*  redo failed populations *###
###*############################################*###
# merrilli_WA
# insignis_AK
# morphna_BC

from concurrent.futures import ProcessPoolExecutor
from pythonimports import *
import subprocess
import os
import re

# specify failed pop names
failed_populations = ['merrilli_WA', 'insignis_AK', 'morphna_BC']

# function to run GF training R script
def run_r_script(population_name):
    script_file = '/workdir/kcarbeck/02.gradient_training.R'
    snpfile = f'/workdir/kcarbeck/data/loocv/snpfile_snp_train_{population_name}.txt'
    envfile = f'/workdir/kcarbeck/data/loocv/env_train_pop_{population_name}.csv'
    range_file = '/workdir/kcarbeck/data/2001_2011/west_combined_raster_data_2001_2011.txt'
    basename = f'west_loocv_{population_name}'
    save_dir = '/workdir/kcarbeck/data/loocv/training'
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


# use ProcessPoolExecutor to run subprocesses in parallel
max_workers = 22
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    # submit a job only for the failed population names
    futures = [executor.submit(run_r_script, name) for name in failed_populations]
    # wait for all jobs to complete
    for future in futures:
        future.result()
