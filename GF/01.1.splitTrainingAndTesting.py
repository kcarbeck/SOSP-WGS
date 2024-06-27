# split into training and test sets to use for cross validation
# katherine carbeck
# 5 dec 2023

# maf file = WZA_Outlier_SNPs_MAFs.csv
# env file = scaledPopEnvAnnualUpdatedwLATLONG.csv


# create training data sets - one for each kfold, and one 'full' model using all pops

#! bash
export PYTHONPATH=/home/kcarbeck/.local/lib64/python3.9/site-packages:$PYTHONPATH
python
#! 


import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

# Load your dataframes
df1 = pd.read_csv("WZA_Outlier_SNPs_MAFs.csv")
df2 = pd.read_csv("scaledPopEnvAnnualUpdatedwLATLONG.csv")

# Set 'population' column as the index
df1.set_index('population', inplace=True)
df2.set_index('Population', inplace=True)

# Assuming 'population' column is common to both dataframes
populations = df1.index.unique()

# Set the number of desired folds (K)
k_fold = 5
kf = KFold(n_splits=k_fold, shuffle=True, random_state=42)

# Iterate through each fold
for fold, (train_index, test_index) in enumerate(kf.split(populations)):
    train_populations = populations[train_index]
    test_populations = populations[test_index]
    print(f"Fold {fold + 1}:")
    print("Training populations:")
    print(train_populations)
    print("Testing populations:")
    print(test_populations)
    # Subset your dataframes for training and testing based on populations
    train_fold_df1 = df1.loc[df1.index.isin(train_populations)]
    test_fold_df1 = df1.loc[df1.index.isin(test_populations)]
    train_fold_df2 = df2.loc[df2.index.isin(train_populations)]
    test_fold_df2 = df2.loc[df2.index.isin(test_populations)]
    # Save training and testing sets to files for df1
    train_fold_df1.to_csv(f"df1_training-k{fold + 1}.txt")
    test_fold_df1.to_csv(f"df1_testing-k{fold + 1}.txt")
    # Save training and testing sets to files for df2
    train_fold_df2.to_csv(f"df2_training-k{fold + 1}.txt")
    test_fold_df2.to_csv(f"df2_testing-k{fold + 1}.txt")
    # For simplicity, print a separator line
    print('-' * 40)

# Save the full dataframes with the index as the population column
df1.to_csv("full_df1.txt")
df2.to_csv("full_df2.txt")


quit()


#! bash
# copy files to home dir
cp df1* /home/lc736_0001/song_sparrow/final_vcf/GF/datasets &
cp df2* /home/lc736_0001/song_sparrow/final_vcf/GF/datasets &
cp full* /home/lc736_0001/song_sparrow/final_vcf/GF/datasets &

# maf file = WZA_Outlier_SNPs_MAFs.csv = df1
# env file = scaledPopEnvAnnualUpdatedwLATLONG.csv = df2

#############!
Fold 1:
Training populations:
Index(['adusta_MX', 'caurina_AK', 'cleonensis_CA', 'fallax_AZ', 'fallax_CA',
       'fallax_UT', 'gouldii_CA', 'graminea_CA', 'insignis_AK',
       'kenaiensis_AK', 'maxima_AK', 'merrilli_AK', 'merrilli_WA',
       'montana_NV', 'montana_N_CA', 'montana_OR', 'montana_S_CA',
       'morphna_BC', 'nominate_VA', 'pusillula_CA', 'rivularis_MX',
       'rufina_BC', 'sanaka_AK'],
      dtype='object', name='population')
Testing populations:
Index(['heermanni_N_CA', 'heermanni_S_CA', 'maxillaris_CA', 'mexicana_MX',
       'nominate_ON', 'samuelis_CA'],
      dtype='object', name='population')
----------------------------------------
Fold 2:
Training populations:
Index(['caurina_AK', 'cleonensis_CA', 'fallax_AZ', 'fallax_CA', 'fallax_UT',
       'gouldii_CA', 'graminea_CA', 'heermanni_N_CA', 'heermanni_S_CA',
       'insignis_AK', 'maxillaris_CA', 'merrilli_AK', 'mexicana_MX',
       'montana_NV', 'montana_N_CA', 'montana_OR', 'montana_S_CA',
       'nominate_ON', 'nominate_VA', 'pusillula_CA', 'rivularis_MX',
       'samuelis_CA', 'sanaka_AK'],
      dtype='object', name='population')
Testing populations:
Index(['adusta_MX', 'kenaiensis_AK', 'maxima_AK', 'merrilli_WA', 'morphna_BC',
       'rufina_BC'],
      dtype='object', name='population')
----------------------------------------
Fold 3:
Training populations:
Index(['adusta_MX', 'fallax_AZ', 'gouldii_CA', 'graminea_CA', 'heermanni_N_CA',
       'heermanni_S_CA', 'insignis_AK', 'kenaiensis_AK', 'maxillaris_CA',
       'maxima_AK', 'merrilli_AK', 'merrilli_WA', 'mexicana_MX',
       'montana_N_CA', 'montana_OR', 'montana_S_CA', 'morphna_BC',
       'nominate_ON', 'nominate_VA', 'rivularis_MX', 'rufina_BC',
       'samuelis_CA', 'sanaka_AK'],
      dtype='object', name='population')
Testing populations:
Index(['caurina_AK', 'cleonensis_CA', 'fallax_CA', 'fallax_UT', 'montana_NV',
       'pusillula_CA'],
      dtype='object', name='population')
----------------------------------------
Fold 4:
Training populations:
Index(['adusta_MX', 'caurina_AK', 'cleonensis_CA', 'fallax_CA', 'fallax_UT',
       'gouldii_CA', 'graminea_CA', 'heermanni_N_CA', 'heermanni_S_CA',
       'insignis_AK', 'kenaiensis_AK', 'maxillaris_CA', 'maxima_AK',
       'merrilli_AK', 'merrilli_WA', 'mexicana_MX', 'montana_NV', 'montana_OR',
       'morphna_BC', 'nominate_ON', 'pusillula_CA', 'rufina_BC',
       'samuelis_CA'],
      dtype='object', name='population')
Testing populations:
Index(['fallax_AZ', 'montana_N_CA', 'montana_S_CA', 'nominate_VA',
       'rivularis_MX', 'sanaka_AK'],
      dtype='object', name='population')
----------------------------------------
Fold 5:
Training populations:
Index(['adusta_MX', 'caurina_AK', 'cleonensis_CA', 'fallax_AZ', 'fallax_CA',
       'fallax_UT', 'heermanni_N_CA', 'heermanni_S_CA', 'kenaiensis_AK',
       'maxillaris_CA', 'maxima_AK', 'merrilli_WA', 'mexicana_MX',
       'montana_NV', 'montana_N_CA', 'montana_S_CA', 'morphna_BC',
       'nominate_ON', 'nominate_VA', 'pusillula_CA', 'rivularis_MX',
       'rufina_BC', 'samuelis_CA', 'sanaka_AK'],
      dtype='object', name='population')
Testing populations:
Index(['gouldii_CA', 'graminea_CA', 'insignis_AK', 'merrilli_AK',
       'montana_OR'],
      dtype='object', name='population')

----------------------------------------