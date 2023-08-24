import numpy as np
import pandas as pd
import os
import sys
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management

def compute_phases_distribution(data_frame, theta_ryhtm=125.0,binsize=15):
    iseed = np.unique(data_frame["iseed"])
    bseed = np.unique(data_frame["ibseed"])
    binphase = np.arange(0,360,15)

    list_phases = []
    data = dict.fromkeys(["counts","iseed","ibseed"])
    for key in data.keys():
        data[key] = []

    for i in iseed:
        for j in bseed: 
            w = (data_frame["iseed"]==i) & (data_frame["ibseed"]==j)
            spikes = data_frame[w]["tvec"].values
            spikes = spikes[spikes>2000]
            phases = np.mod( spikes-50, theta_ryhtm)*360/theta_ryhtm
            counts, bins = np.histogram( phases, bins=binphase, density="True" )
            list_phases.append(counts)
            # print(i,j, len(counts))

            data["counts"].append(counts)
            data["iseed"].append([i]*len(counts))
            data["ibseed"].append([j]*len(counts))
    
    data["counts"] = np.concatenate(data["counts"])
    data["iseed"] = np.concatenate(data["iseed"])
    data["ibseed"] = np.concatenate(data["ibseed"])

    data = pd.DataFrame(data)

    data2 = dict.fromkeys(["binphase","counts_mean","counts_sem"])
    data2["binphase"] = bins[1:]-np.diff(bins)[0]/2
    data2["counts_mean"] = np.mean(list_phases, axis=0)
    data2["counts_sem"]  = np.std(list_phases, axis=0)/np.sqrt(len(list_phases))

    data2 = pd.DataFrame(data2)
    output = {}
    output["realizations"] = data
    output["average"] = data2
    
    return output

def compute_phases_from_spikes_ca3(folder):
    theta_ryhtm = 125.0
    binsize = 15
    parent_dir = os.path.join(folder, "measurements")
    spikes_files = ["external_inputs_pyr_ca3.lzma","spikes_bas_ca3.lzma","spikes_olm_ca3.lzma"]

    aux = ["pyr","bas","olm"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem") # std/np.sqrt(n)
    data_frame = dict.fromkeys(columns) 

    outputs = []
    # pyramidal cells
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[0]) )
    outputs.append( compute_phases_distribution(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # basket cells
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[1]) )
    outputs.append( compute_phases_distribution(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # olm cells
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[2]) )
    outputs.append( compute_phases_distribution(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    for i, key in enumerate(aux):
        data_frame[key+"_mean"] = outputs[i]["average"]["counts_mean"]
        data_frame[key+"_sem"] = outputs[i]["average"]["counts_sem"]

    file_management.save_lzma(data_frame, "phases_spikes_distribution_ca3.lzma", parent_dir=parent_dir)


def compute_phases_from_spikes_ca1(folder):
    theta_ryhtm = 125.0
    binsize = 15
    parent_dir = os.path.join(folder, "measurements")

    spikes_files = ["spikes_pyr_ca1.lzma","spikes_bas_ca1.lzma","spikes_olm_ca1.lzma","spikes_cck_ca1.lzma"]

    aux = ["pyr","bas","olm","cck"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem") # std/np.sqrt(n)
    data_frame = dict.fromkeys(columns)

    outputs = []
    # pyramidal cells
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[0]) )
    outputs.append( compute_phases_distribution(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # basket cells
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[1]) )
    outputs.append( compute_phases_distribution(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # olm cells
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[2]) )
    outputs.append( compute_phases_distribution(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # cck cells
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[3]) )
    outputs.append( compute_phases_distribution(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    for i, key in enumerate(aux):
        data_frame[key+"_mean"] = outputs[i]["average"]["counts_mean"]
        data_frame[key+"_sem"]  = outputs[i]["average"]["counts_sem"]

    file_management.save_lzma(data_frame,"phases_spikes_distribution_ca1.lzma",parent_dir=parent_dir)