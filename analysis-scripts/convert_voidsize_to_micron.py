# convert voidsize from sim units to microns

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

hb = 20 #[micron]
R_C = 0.55 #[micron]
K_scale = 1

if os.path.exists('void_data/') == False:
  os.mkdir('void_data/')

brushes = ['H', 'L']
salts = ['2.5', '3.5', '5', '10']
depletants = ['0.5', '1', '2', '3', '4', '6', '8']


#################################
""" LOOP THROUGH ALL DATASETS """
#################################

all_data_list = []
for brush in brushes:

  if brush == 'H':
    base = Path('K_'+str(K_scale)+'xUC/brush'+str(brush))
  if brush == 'L':
    base = Path('K_'+str(K_scale)+'xUC/brush'+str(brush))

  for salt in salts:
    salt_folder = brush+'_'+salt+'mM'

    # create "salt" subfolder if it doesn't exit
    if os.path.exists('void_data/'+salt_folder) == False:
      os.mkdir('void_data/'+salt_folder)

    nvoids_list = []
    volumes_list = []
    T_avg_list = []
    G_avg_list = []
    depletant_list = []
    for depletant in depletants:
        depletant_list.append(depletant)
        tag = depletant+'mg/mL'

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
            csv_path = folder / "voidsize.csv"
            #print("Checking:", csv_path)
            if csv_path.exists():
              data_df = pd.read_csv(csv_path)
              #data_df["depletant"] = f"{depletant}mg/mL"
              #print("Loaded:", csv_path)
              break  # stop after first match

        if data_df is None:
          print("No matching CSV found!")
          T_avg_list.append('N/A')  
          G_avg_list.append('N/A')  
          volumes_list.append('N/A')  
          nvoids_dataset = 'N/A'
          nvoids_list.append(nvoids_dataset)
          continue

        data_df['depletant'] = tag

        print(f"{brush}%, {salt}mM, {depletant}mg/mL")

        # calculate system volume and add to the dataframe and convert to microns
        L_X = int(np.ceil(max(data_df['probe_posx']))-np.floor(min(data_df['probe_posx']))) * R_C
        L_Y = int(np.ceil(max(data_df['probe_posy']))-np.floor(min(data_df['probe_posy']))) * R_C
        L_Z = int(np.ceil(max(data_df['probe_posz']))-np.floor(min(data_df['probe_posz']))) * R_C
        volumes_list.append(L_X*L_Y*L_Z)
        print(f" - System volume: {[L_X, L_Y, L_Z]}")

        # create a single dataset of all sim data and convert to microns
        dataset_df = data_df.copy()
        dataset_df['probe_posx'] = data_df['probe_posx'] * R_C
        dataset_df['probe_posy'] = data_df['probe_posy'] * R_C
        dataset_df['probe_posz'] = data_df['probe_posz'] * R_C
        dataset_df['voidcenter_x'] = data_df['voidcenter_x'] * R_C
        dataset_df['voidcenter_y'] = data_df['voidcenter_y'] * R_C
        dataset_df['voidcenter_z'] = data_df['voidcenter_z'] * R_C
        dataset_df['voiddiameter_T'] = data_df['voiddiameter_T'] * R_C
        dataset_df['voiddiameter_G'] = data_df['voiddiameter_G'] * R_C
        dataset_df['salt'] = salt
        dataset_df['monomer'] = brush
        dataset_df['dims_x'] = L_X
        dataset_df['dims_y'] = L_Y
        dataset_df['dims_z'] = L_Z
        all_data_list.append(dataset_df)

        for voidsize_type in ['voiddiameter_T','voiddiameter_G']:
          print(f' - {voidsize_type}')

          all_void_diameters = dataset_df[voidsize_type].to_numpy()

          # [optional] count NaNs and Infs
          total_count = len(all_void_diameters)
          nan_count = np.isnan(all_void_diameters).sum()
          inf_count = np.isinf(all_void_diameters).sum()
          nan_inf_percentage = (nan_count + inf_count) / total_count * 100
          print(f"   NaNs: {nan_count}, Infs: {inf_count}, Total: {total_count}")
          print(f"   Percentage: {nan_inf_percentage:.2f}%")

          #################################
          """ AVG PORE SIZE CALCULATION """
          #################################

          # get average
          void_diameters = all_void_diameters[np.isfinite(all_void_diameters)] # clean any NaNs and infs
          void_avg = np.mean(void_diameters)
          if voidsize_type == 'voiddiameter_T':
            T_avg_list.append(void_avg)
          elif voidsize_type == 'voiddiameter_G':
            G_avg_list.append(void_avg)
          else:
            print(f"ERROR: sizetype must be 'Torquato' or 'Gubbin', not: {sizetype}")
            exit(1)

          # calc nvoids for future normalization
          nvoids_dataset = len(void_diameters)

        # only add nvoids once, it's the same for Torquato and Gubbin
        nvoids_list.append(nvoids_dataset)

    avgs_df = pd.DataFrame(depletant_list, columns=['depletant'])
    avgs_df['T_avg'] = T_avg_list
    avgs_df['G_avg'] = G_avg_list
    avgs_df['nvoids'] = nvoids_list
    avgs_df['volume'] = volumes_list

    avgs_df.to_csv('void_data'+salt_folder+'/'+'voidsize-avg-microns.csv', index=False)
all_data_df = pd.concat(all_data_list, ignore_index=True)
all_data_df.to_csv('void_data/voidsize_allsim_datasets_microns.csv', index=False)

print("\nAll analyses complete for")
# record params:
print("R_C conversion:",R_C,"micron")
print("brushes:",brushes)
print("salts:",salts,"mM")
print("depletants:",depletants,"mg/mL")
print("hb:",hb,"nm")
