import numpy as np
import pandas as pd
import math
import re
import os
import glob
from pathlib import Path

brushes = ['H', 'L']
salts = ['2.5', '3.5', '5', '10']
depletants = ['0.5', '1', '2', '3', '4', '6', '8']
hb_values = [ '20' ]
K_scale = 1

# choose which frames to plot (default last frame)
framechoice = [-1]


# load and process simulation data from LAST FRAME ONLY (networkx data)
df_list = []
for brush in brushes:

  if brush == 'H':
    base = Path(sim_path_base+'K_'+str(K_scale)+'xUC/brush'+str(brush))
  if brush == 'L':
    base = Path(sim_path_base+'K_'+str(K_scale)+'xUC/brush'+str(brush))

  for salt in salts:
    for depletant in depletants:
      for hb_val in hb_values:
          void_filepath = 'void_data/'+brush+'_'+salt+'mM'+'/voidsize-avg-microns.csv'
          lastframe_voidsize_avg = pd.read_csv(void_filepath)

          allframes_sim_df = None
          # format the parts exactly as in your folder names
          if brush == 'L':
            prefix = f"data_{depletant}mgmL-{salt}mM-K"
          elif brush == 'H':
            prefix = f"data_{depletant}mgmL-{salt}mM-K"
          suffix = f"_{brush}"

          # scan the base directory
          print(base)
          for folder in base.iterdir():
            if folder.is_dir() and folder.name.startswith(prefix) and folder.name.endswith(suffix):
              sim_filepath = folder / "networkx-allframes.csv"
              #print("Checking:", sim_filepath)
              if sim_filepath.exists():
                allframes_sim_df = pd.read_csv(sim_filepath)
                #data_df["depletant"] = f"{depletant}mg/mL"
                #print("Loaded:", sim_filepath)
                break  # stop after first match

          if allframes_sim_df is None:
            print("No matching CSV found!")
            break

          # only use last frame
          frames = pd.unique(allframes_sim_df['frame'])
          lastframe = frames[-1]
          sim_df = allframes_sim_df.loc[allframes_sim_df['frame'] == lastframe]

          # only use current voidsize value
          voidsize_df = lastframe_voidsize_avg.loc[lastframe_voidsize_avg['depletant'] == float(depletant)]

          print(brush, depletant, salt)

          # create the dataframe : brush,salt,depletant,Zavg,lcc,percolation,void_size
          out_df = pd.DataFrame({
            "brush": [brush],
            "salt": [salt],
            "depletant": [depletant],
            "Zavg": [sim_df['avg_degree'].values[0]],
            "lcc": [sim_df['lcc_size'].values[0] / sim_df['ncolloids'].values[0]],
            "percolation": [sim_df['percolation'].values[0]],
            "void_size": [voidsize_df['G_avg'].values[0]],
          })

          df_list.append(out_df)

all_data_df = pd.concat(df_list, ignore_index=True)
out_filename = 'sim_Zavg_lcc_data_microns_K_'+str(K_scale)+'xUC.csv'
all_data_df.to_csv(out_filename, index=False)

print("All data files created.")


