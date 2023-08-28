import numpy as np 
import pandas as pd
import os
import sys
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management 

current_folder = os.getcwd()    
data_folder = os.path.join(current_folder,"test_data")  
# data_folder = '/home/jaime/Desktop/hippocampus/external_inputs/baseline/'
data_folder = os.path.join('/home/jaime/Desktop/hippocampus/agosto2023/',sys.argv[1])
# Create a list of all the filenames
# Initial names of all the files we are storing 
# If no file start with these names, nothing happens, no error is raised
# A message tells us that no files were found with this name

# It loads and combines all the files that start with the same name
# If this code is rerun, there will be only one file with matches the name
# and then the lzma file will be deleted and stored again, so nothing happens
filenames = []
filenames.extend(["volt_pyr_ca3_","volt_bas_ca3_","volt_olm_ca3_"])
filenames.extend(["volt_pyr_ca1_","volt_bas_ca1_","volt_olm_ca1_","volt_cck_ca1_"])
filenames.extend(["external_inputs_pyr_ca3_","spikes_bas_ca3_","spikes_olm_ca3_"])
filenames.extend(["spikes_pyr_ca1_","spikes_bas_ca1_","spikes_olm_ca1_","spikes_cck_ca1_"])
filenames.extend(["synapses_pyr_ca3_"])
filenames.extend(["synapses_pyr_ca1_"])
filenames.extend(["lfp_ca3_"])
filenames.extend(["lfp_ca1_"])
filenames.extend(["ica_ca3_"])
filenames.extend(["ica_ca1_"])
# filenames.extend([ ... ]) if necessary

print("Creating the groups of files")
print("----------------------------")
groups = [] 
groups_name = []
k = 1
for name in filenames:
    groups.append([])
    for file in os.listdir(data_folder):
        if name in file:
            groups[-1].append(file)
    if not groups[-1]:
        print(f"{k}. No files found with {name}.")
        k+=1

print(" ")
print("Combining the files and removing them")
print("-------------------------------------")
k=1
for i,group in enumerate(groups):
    if group:
        combined_df = pd.DataFrame()
        for file in group:
            print(file)
            df = file_management.load_lzma(os.path.join(data_folder,file))
            combined_df = pd.concat([combined_df,df],ignore_index=True)
            # os.remove(os.path.join(data_folder,file))
        print(f"{k}.All data in files with {filenames[i][:-1]} were combined and stored.")
        file_management.save_lzma(combined_df, f"{filenames[i][:-1]}.lzma", parent_dir=data_folder)
        print("Proceding to detelete the files...")
        for file in group:
            os.remove(os.path.join(data_folder,file))
        print(f"{k}.All data in files with {filenames[i][:-1]} were deleted.")
        k+=1
        print(" ")
