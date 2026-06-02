# convert Z counts to probabilities

########################
""" MODULE LIBRARY """
########################
import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path

##########################
""" INPUT PARAMETERS """
##########################

hb = '20' #str(sys.argv[2])
R_C = 0.55 #[micron]
K_scale = 1

if os.path.exists('Z_data/') == False:
  os.mkdir('Z_data/')

brushes = ['H', 'L']
salts = ['2.5', '3.5', '5', '10']
depletants = ['0.5', '1', '2', '3', '4', '6', '8']


#################################
""" LOOP THROUGH ALL DATASETS """
#################################

tag_list = []
all_counts_list = []
for brush in brushes:

  if brush == 'H':
    base = Path('K_'+str(K_scale)+'xUC/brush'+str(brush))
  if brush == 'L':
    base = Path('K_'+str(K_scale)+'xUC/brush'+str(brush))

  for salt in salts:
    salt_folder = brush+'_'+salt+'mM'

    # create "salt" subfolder if it doesn't exit
    if os.path.exists('Z_data/'+salt_folder) == False:
      os.mkdir('Z_data/'+salt_folder)

    for depletant in depletants:
        tag = brush+'_'+depletant+'mgmL_'+salt+'mM'

        data_df = None
        # format the parts exactly as in your folder names
        if brush == 'L':
          prefix = f"data_{depletant}mgmL-{salt}mM-K"
        elif brush == 'H':
          prefix = f"data_{depletant}mgmL-{salt}mM-K"
        suffix = f"_{brush}"

        print("Looking for folders starting with:", prefix)
        print("And ending with:", suffix)

        # scan the base directory
        for folder in base.iterdir():
          if folder.is_dir() and folder.name.startswith(prefix) and folder.name.endswith(suffix):
            csv_path = folder / "Z-counts.csv"
            #print("Checking:", csv_path)
            if csv_path.exists():
              data_df = pd.read_csv(csv_path)
              #data_df["depletant"] = f"{depletant}mg/mL"
              #print("Loaded:", csv_path)
              break  # stop after first match

        if data_df is None:
          print("No matching CSV found!")
          continue

        frames = pd.unique(data_df['frame'])
        lastframe = max(frames)
        lastframe_df = data_df.loc[data_df['frame'] == lastframe].copy()
        
        lastframe_df['brush'] = brush
        lastframe_df['salt'] = salt
        lastframe_df['depletant'] = depletant
        lastframe_df['tag'] = tag

        tag_list.append(tag)
        all_counts_list.append(lastframe_df)

all_counts_df = pd.concat(all_counts_list, ignore_index=True)

# make a list of bins for the x-axis
max_Z_list = []
for tag in tag_list:
  tag_df = all_counts_df.loc[all_counts_df['tag'] == tag]
  tag_maxZ = max(tag_df['Z'].to_numpy())
  max_Z_list.append(tag_maxZ)
true_max_Z = int(max(max_Z_list))
bins = np.linspace(0,true_max_Z,true_max_Z+1).astype(int)

# bin all the data!
binned_data_list = []
for brush in brushes:

  for salt in salts:
    salt_folder = brush+'_'+salt+'mM'

    # create "salt" subfolder if it doesn't exit
    if os.path.exists('Z_data/'+salt_folder) == False:
      os.mkdir('Z_data/'+salt_folder)

    for depletant in depletants:
        tag = brush+'_'+depletant+'mgmL_'+salt+'mM'

        dataset_df = all_counts_df.loc[all_counts_df['tag'] == tag]
        ncolloids = len(dataset_df['colloidID'])
        #print('# ',ncolloids)
        df = dataset_df.copy()

        histbins = np.zeros(int(true_max_Z)+1)
        probs = np.zeros(int(true_max_Z)+1)

        # bin the Zs
        for i in range(true_max_Z+1):
          Zs = df['Z'].to_numpy()
          for z in range(len(df['Z'])):
            if Zs[z] == i:
              histbins[i] += 1

        # convert to probability
        probs[:] = histbins[:]/ncolloids 

        hist_df = pd.DataFrame(bins, columns=['bins'])
        hist_df['counts'] = histbins
        hist_df['probs'] = probs
        hist_df['brush'] = brush
        hist_df['salt'] = salt
        hist_df['depletant'] = depletant
        hist_df['tag'] = tag
        binned_data_list.append(hist_df)

binned_data_df = pd.concat(binned_data_list, ignore_index=True)
binned_data_df.to_csv('Z_data/zdist_allsim_datasets_binned.csv', index=False)

print("\nAll analyses complete for")
# record params:
print("R_C conversion:",R_C,"micron")
print("brushes:",brushes)
print("salts:",salts,"mM")
print("depletants:",depletants,"mg/mL")
print("hb:",hb,"nm")
