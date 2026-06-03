## Calculate voidsize distribution from non-periodic CSV files 
## NOTE: requires matching Fortran module and solvopt module
## NOTE: this code assumes 1 colloid type
"""
## This code performs the following analyses for data WITHOUT periodic
## boundary conditions (PBC): 
##    * void size distribution using two methods: 
##       - Torquato’s Pore Size Distribution
##       - Gubbins’s Pore Size Distribution 
##       (requires the solvopt algorithm, included as a separate f90 module)
##       NOTE: this version assumes there are multiple replicates 
##             that should be merged
##       NOTE: requires fortran module
##    * trim outlier voids from the edge of a system without PBC 
##    * take the average of the cleaned data
"""
## (Rob Campbell)

########################
""" MODULE LIBRARY """
########################
import numpy as np
import pandas as pd
import os


##########################
""" INPUT PARAMETERS """
##########################

data_path_base = 'csv_data/'

filename_base = '10S'

n_replicates = 3
brushes = ['H', 'L']
salts = ['2.5', '3.5', '5', '10']
depletants = ['0.5', '1', '2', '3', '4', '6', '8']

R_C = 0.55

# filepath to folder where data files will be created
data_outpath = 'void_data'

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)

# filename for cleaned void size data
cleaned_out_filename = 'cleaned_voidsize_all_datasets.csv'


###########################
""" DEFINE SIM CHECKS """
###########################

#######
## calculate void size for all datasets
"""
# calculates a distribution of the void space (approximated as spheres) in between particle clusters
#
# NOTE: without periodic boundary conditions (PBC) the size of voids at the edges of the system will blow up. Trim data below after generating the raw sizes.
#
# There are two methods that are used:
#   - Torquato’s Pore Size Distribution (volume where the center of a particle can fit in between clusters)
#   - Gubbins’s Pore Size Distribution (volume occupied by a whole particle in between clusters) 
#       (Gubbin's PSD requires solvopt, the Solver For Local Nonlinear Optimization Problems, 
#        included as a secondary fortran module) 
# 
# These methods are described in the Section 4.2 and the Appendix of
# Sorichetti, Hugouvieux, and Kob 2020, DOI: 10.1021/acs.macromol.9b02166
#
# This code assumes you are analyzing a porous medium made of uniform particles (like a colloidal gel),
#  and uses particle trajectories in gsd format; It takes particle and probe size as inputs.
#  Then it uses Linked-list method to compute minimum distance of a point in the void space from 
#  nearby porous medium particles. And then uses solvopt non-linear optimization code to compute Gubbin's pore size.
#
# (we usually use Gubbin's PSD and expect a plot of probability vs void diameter to peak at the size of the most common voids)
""" 
def calculate_voidsize(data_path_base, filename_base, n_replicates, brushes, salts, depletants, R_C, data_outpath):

  from fortranmod import void_size_calculation

  for brush in brushes:
    for salt in salts:
      for depletant in depletants:

        ### PREP DATA

        salt_outfolder = brush+'_-'+salt+'mM'

        # create "salt" subfolder if it doesn't exit
        if os.path.exists(data_outpath+'/'+salt_outfolder) == False:
          os.mkdir(data_outpath+'/'+salt_outfolder)


        file_name = filename_base+'_'+brush+'__'+salt+'mM_'+depletant+'mg'

        # compile csv data into a single dataframe (labeled by dataset)
        df_list = []
        for d in range(n_replicates):
          filepath = data_path_base+'/'+salt_outfolder+'/'+file_name+'_'+str(d+1)+'.csv'
          df = pd.read_csv(filepath)
          df['dataset'] = d
          df_list.append(df)

        data_df = pd.concat(df_list, ignore_index=True)
 
        #### VOID SIZE CALCULATION
        # set the number of random points used to explore void size
        nprobe = 10000 # can test quickly at 1
        # set the number of random points used to explore void size
        dcell_init = 5.0

        # check if data file already exists
        if os.path.exists(data_outpath+'/'+salt_outfolder+'/'+file_name+'_voidsize.csv') == True:
          print(f'ERROR: {file_name}_voidsize.csv already exists')
          print('the voidsize calc appends data, it does not overwrite files')
          print('remove/rename the existing file, or comment out this warning to append to it')
          exit()
        else:
          print('Calculating voidsize distribution for '+str(n_replicates)+' '+str(brush)+'_'+str(salt)+'mM-'+str(depletant)+'mg datasets...')

        ## calculate for each dataset individually and then merge (ignore voidsize dataset labels
        ## and treat entire file as one dataset)

        for d in range(n_replicates):
          # select only the data from the current dataset
          curr_data_df = data_df.loc[data_df['dataset'] == d]

          # calculate the number of colloids
          ncolloids = len(curr_data_df)

          # find the radii of every type1 colloid
          radii = curr_data_df['Radius'].to_numpy() #/ R_C

          # get the box size
          L_X = np.ceil(np.max(curr_data_df['x'])) #/ R_C
          L_Y = np.ceil(np.max(curr_data_df['y'])) #/ R_C
          L_Z = np.ceil(np.max(curr_data_df['z'])) #/ R_C
          box_length = np.array([L_X,L_Y,L_Z]) 
          rxi = curr_data_df['x'].to_numpy() - L_X/2 #/ R_C
          ryi = curr_data_df['y'].to_numpy() - L_Y/2 #/ R_C
          rzi = curr_data_df['z'].to_numpy() - L_Z/2 #/ R_C

          print('...calculating dataset '+str(d))

          # calculate the void_size for the current dataset
          # write to the voidsize.csv file created or appended to in fortran module
          void_size_calculation.void_size_calc(data_outpath+'/'+salt_outfolder+'/'+file_name+'_voidsize.csv',d,ncolloids,nprobe,radii,dcell_init,rxi,ryi,rzi,box_length)

      print("Pore size distribution calculated for "+str(nprobe)+" probe points in each of "+str(n_replicates)+" "+str(brush)+'_'+str(salt)+'mM-'+str(depletant)+'mg datasets...')

  print("Poresize calculation complete for all systems\n")


#######
## clean non-PBC void size data
"""
# clean all voidsize data by removing falsely identified voids located outside
# the imaged system volume (or within one colloid radii of the edge of the 
# system volume, similar to how particle tracking data is cleaned)

# assumes the following data structure using brush, salt, and depletant values:
#   data_outpath+'H_-10mM/10S_H__10mM_8mg_voidsize.csv'
"""
def clean_voidsize_data(data_outpath, filename_base, cleaned_out_filename, brushes, salts, depletants, R_C):

  ### LOAD DATA 
  df_list = []
  for brush in brushes:
    salt_df_list = []
    for salt in salts:
      salt_folder = brush+'_-'+salt+'mM'

      depletant_df_list = []
      for depletant in depletants:
        tag = depletant+'mg/mL'

        file_name = filename_base+'_'+brush+'__'+salt+'mM_'+depletant+'mg'
        filepath = data_outpath+'/'+salt_folder+'/'+file_name+'_voidsize.csv'
        df = pd.read_csv(filepath)
        df['depletant'] = tag

        depletant_df_list.append(df)

      salt_df = pd.concat(depletant_df_list, ignore_index=True)
      salt_df['salt'] = salt
      salt_df_list.append(salt_df)

    brush_df = pd.concat(salt_df_list, ignore_index=True)
    brush_df['monomer'] = brush
    df_list.append(brush_df)

  raw_data_df = pd.concat(df_list, ignore_index=True)

  ### CLEAN DATA

  # Function to extend data frames if they are not equal
  def extend_data_frames(base_df, target_len):
    if len(base_df) < target_len:
      extension_df = pd.DataFrame(index=range(len(base_df), target_len))
      for col in base_df.columns:
        unique_values = base_df[col].dropna().unique()
        if len(unique_values) == 1:
          extension_df[col] = unique_values[0]  # Use the single value
        else:
          extension_df[col] = np.nan  # Use NaN for varying columns
      extended_df = pd.concat([base_df, extension_df], ignore_index=True, sort=False)
      return extended_df
    else:
      return base_df

  # remove voids located outside the real system dimensions
  cleaned_data_list = []
  for brush in brushes:
    brush_df = raw_data_df.loc[raw_data_df['monomer'] == brush]

    for salt in salts:
      salt_df = brush_df.loc[brush_df['salt'] == salt]

      for depletant in depletants:
        tag = depletant+'mg/mL'
        dataset_df = salt_df.loc[salt_df['depletant'] == tag].copy()
        print(f"{brush}%, {salt}mM, {depletant}mg/mL")

        # calculate system volume and add to the dataframe
        L_X = int(np.ceil(max(dataset_df['probe_posx']))-np.floor(min(dataset_df['probe_posx']))) # microns
        L_Y = int(np.ceil(max(dataset_df['probe_posy']))-np.floor(min(dataset_df['probe_posy']))) # microns
        L_Z = int(np.ceil(max(dataset_df['probe_posz']))-np.floor(min(dataset_df['probe_posz']))) # microns
        dataset_df['dims_x'] = L_X
        dataset_df['dims_y'] = L_Y
        dataset_df['dims_z'] = L_Z
        print(f" - System volume: {[L_X, L_Y, L_Z]}")

        # clean data so we only consider voids with a center inside the system volume
        # (the same criteria that was used for colloid particle tracking)
        min_x = round(min(dataset_df['probe_posx']))+R_C
        max_x = round(max(dataset_df['probe_posx']))-R_C
        min_y = round(min(dataset_df['probe_posy']))+R_C
        max_y = round(max(dataset_df['probe_posy']))-R_C
        min_z = round(min(dataset_df['probe_posz']))+R_C
        max_z = round(max(dataset_df['probe_posz']))-R_C

        df_filtered = dataset_df.loc[
                       (dataset_df['voidcenter_x'] >= min_x) & (dataset_df['voidcenter_x'] <= max_x) &
                       (dataset_df['voidcenter_y'] >= min_y) & (dataset_df['voidcenter_y'] <= max_y) &
                       (dataset_df['voidcenter_z'] >= min_z) & (dataset_df['voidcenter_z'] <= max_z)
                       ]
        percentage = len(df_filtered) / len(dataset_df) * 100

        outlier_voids = dataset_df.loc[~(
                       (dataset_df['voidcenter_x'] >= min_x) & (dataset_df['voidcenter_x'] <= max_x) &
                       (dataset_df['voidcenter_y'] >= min_y) & (dataset_df['voidcenter_y'] <= max_y) &
                       (dataset_df['voidcenter_z'] >= min_z) & (dataset_df['voidcenter_z'] <= max_z)
                       )]
        percentage_removed = len(outlier_voids) / len(dataset_df) * 100

        print(f" - removed {len(outlier_voids)} void centers outside the box ({round(percentage_removed,2)}%)")
        print(f"   {round(percentage,2)}% remaining: {len(df_filtered)} total")

        cleaned_data_list.append(df_filtered)

        for voidsize_type in ['voiddiameter_T','voiddiameter_G']:
          void_diameters = df_filtered[voidsize_type].to_numpy()
          print(f" {voidsize_type}:")
          # [optional] count NaNs and Infs
          total_count = len(void_diameters)
          nan_count = np.isnan(void_diameters).sum()
          inf_count = np.isinf(void_diameters).sum()
          nan_inf_percentage = (nan_count + inf_count) / total_count * 100
          print(f" - NaNs: {nan_count}, Infs: {inf_count}, Total: {total_count}")
          print(f"   Percentage: {nan_inf_percentage:.2f}%")
        print('----------')

  cleaned_data_df = pd.concat(cleaned_data_list, ignore_index=True)

  ### SAVE DATA
  cleaned_data_df.to_csv(data_outpath+'/'+cleaned_out_filename, index=False)

  print('All void size data trimmed of edge-based outliers (using the same criteria as the particle tracking data cleaning) and saved to '+data_outpath+'/'+cleaned_out_filename)


#######
## average cleaned void size distribution
"""
## Take the average of the cleaned void size distribution 
"""
def avg_voidsize(data_path_base, filename_base, cleaned_out_filename, brushes, salts, depletants, R_C):

  ### LOOP THROUGH ALL DATASETS
  voiddata_df = pd.read_csv(data_path_base+'/'+cleaned_out_filename)

  for brush in brushes:
    brush_data = voiddata_df.loc[voiddata_df['monomer'] == str(brush)]

    for salt in salts:
      salt_folder = brush+'_-'+salt+'mM'
      salt_data = brush_data.loc[brush_data['salt'] == float(salt)]

      nvoids_list = []
      volumes_list = []
      T_avg_list = []
      T_sem_list = []
      G_avg_list = []
      G_sem_list = []
      depletant_list = []
      for depletant in depletants:
        depletant_list.append(depletant)

        # remember to remove the nans that were added to make the dataframes equal in length
        data_df = salt_data.loc[salt_data['depletant'] == str(depletant)+'mg/mL'].dropna()


        ###  AVG PORE SIZE CALCULATION
        ### treat all replicates as one dataset 

        # calc system volume for future normalization
        L_X_list = pd.unique(data_df['dims_x'])
        L_Y_list = pd.unique(data_df['dims_y'])
        L_Z_list = pd.unique(data_df['dims_z'])
        if ( (len(L_X_list) != 1) | (len(L_Y_list) !=1) | (len(L_Z_list) != 1) ):
          if (len(L_X_list) != 1):
             print(f'WARNING: multiple X system dimensions saved for {brush}_{salt}mM_{depletant}mg: {L_X_list}')
          if (len(L_Y_list) != 1):
             print(f'WARNING: multiple Y system dimensions saved for {brush}_{salt}mM_{depletant}mg: {L_Y_list}')
          if (len(L_Z_list) != 1):
             print(f'WARNING: multiple Z system dimensions saved for {brush}_{salt}mM_{depletant}mg: {L_Z_list}')

          print(f'using average volume across these datasets...')
          volumes_dataset_list = []
          for x in L_X_list:
            for y in L_Y_list:
              for z in L_Z_list:
                sys_volume = x * y * z
                volumes_dataset_list.append(sys_volume)
          vol = np.mean(volumes_dataset_list) # use average
        else:
          vol = L_X_list[0] * L_Y_list[0] * L_Z_list[0]
        volumes_list.append(vol)

        for voidsize_type in ['voiddiameter_T','voiddiameter_G']:
          # get the standard error across the experimental replicates
          dataset_avg_list = []
          for dataset in pd.unique(data_df['frame']):
            dataset_df = data_df.loc[data_df['frame'] == dataset]
            dataset_void_diameters = dataset_df[voidsize_type].to_numpy()
            d_void_diameters = dataset_void_diameters[np.isfinite(dataset_void_diameters)]
            dataset_void_avg = np.mean(d_void_diameters)
            dataset_avg_list.append(dataset_void_avg)
          dataset_sem = np.std(dataset_avg_list, ddof=1) / np.sqrt(len(dataset_avg_list))
          if voidsize_type == 'voiddiameter_T':
            T_sem_list.append(dataset_sem)
          elif voidsize_type == 'voiddiameter_G':
            G_sem_list.append(dataset_sem)
          else:
            print(f"ERROR: voidsize_type must be 'Torquato' or 'Gubbin', not: {voidsize_type}")
            exit(1) 

          # get average from the entire distribution
          all_void_diameters = data_df[voidsize_type].to_numpy()
          void_diameters = all_void_diameters[np.isfinite(all_void_diameters)] # clean any infs
          void_avg = np.mean(void_diameters)
          if voidsize_type == 'voiddiameter_T':
            T_avg_list.append(void_avg)
          elif voidsize_type == 'voiddiameter_G':
            G_avg_list.append(void_avg)
          else:
            print(f"ERROR: voidsize_type must be 'Torquato' or 'Gubbin', not: {voidsize_type}")
            exit(1) 

          # calc nvoids for future normalization
          nvoids_dataset = len(void_diameters)

        # only add nvoids once, it's the same for Torquato and Gubbin
        nvoids_list.append(nvoids_dataset)

        #print("Pore size average calculated for "+str(brush)+'_'+str(salt)+'mM-'+str(depletant)+'mg...')

      avgs_df = pd.DataFrame(depletant_list, columns=['depletant'])
      avgs_df['T_avg'] = T_avg_list
      avgs_df['T_sem'] = T_sem_list
      avgs_df['G_avg'] = G_avg_list
      avgs_df['G_sem'] = G_sem_list
      avgs_df['nvoids'] = nvoids_list
      avgs_df['volume'] = volumes_list

      avgs_df.to_csv(data_path_base+'/'+salt_folder+'/'+filename_base+'_'+brush+'__'+salt+'mM_voidsize-avg.csv', index=False)

      print('Averages saved to: '+data_path_base+'/'+salt_folder+'/'+filename_base+'_'+brush+'__'+salt+'mM_voidsize-avg.csv')




####################
""" RUN ANALYSES """
####################

if __name__ == '__main__':
    # calculate void size
    calculate_voidsize(data_path_base, filename_base, n_replicates, brushes, salts, depletants, R_C, data_outpath)
    # trim false outliers from the edges
    clean_voidsize_data(data_outpath, filename_base, cleaned_out_filename, brushes, salts, depletants, R_C)
    # get average per dataset
    avg_voidsize(data_path_base, filename_base, cleaned_out_filename, brushes, salts, depletants, R_C)

    print("\nAll analyses complete for")
    # record params:
    print("brushes:",brushes,"%")
    print("salts:",salts,"mM")
    print("depletants:",depletants,"mg/mL")
    print("R_C:",R_C, "micron")
