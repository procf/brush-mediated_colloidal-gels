# export all Z distribution data 

import numpy as np
import pandas as pd
import math
import re
import os
import glob



##############
""" INPUTS """
##############

K_scale = 1
hb = 20

data_filename = 'Z_data/zdist_allsim_datasets_binned.csv'

save_dir = "zdist_K1xUC"
if os.path.exists(save_dir) == False:
  os.mkdir(save_dir)

brushes = ['H', 'L']
salts = ['2.5', '3.5', '5', '10']
depletants = ['0.5', '1', '2', '3', '4', '6', '8']


#################
""" LOAD DATA """
#################

# bins,counts,probs,brush,salt,depletant,tag
df = pd.read_csv(data_path+'/'+data_filename)

# load and process simulation data from LAST FRAME ONLY (networkx data)

for brush in brushes:
  brush_df = df.loc[df['brush'] == int(brush)]
  for salt in salts:
    salt_df = brush_df.loc[brush_df['salt'] == float(salt)]
    for depletant in depletants:
      depletant_df = salt_df.loc[salt_df['depletant'] == float(depletant)]

      # create the dataframe : Contact Numbers,Small,Small (Normalized),Total Number of Particles
      out_df = pd.DataFrame({
          "Contact Numbers": depletant_df['bins'],
          "Small": depletant_df['counts'],
          "Small (Normalized)": depletant_df['probs'],
          #"Total Number of Particles": [ncolloids] * len(depletant_df['bins']),
      })

      out_filename = 'distribution_10S_'+brush+'_'+salt+'mM_'+depletant+'mg_simK'+str(K_scale)+'xUC.csv'
      out_df.to_csv(save_dir+'/'+out_filename, index=False)
      print("...",out_filename,"created")

print(f"All data files saved to {save_dir}")
