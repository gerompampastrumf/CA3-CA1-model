import numpy as np
import pandas as pd
import os
import sys
from scipy import stats
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management

def compute_phases_distribution_external_inputs(data_frame, theta_ryhtm=125.0,binsize=15):
    iseed = np.unique(data_frame["iseed"].values)
    binphase = np.arange(0,361,15)
    nbins = len(binphase)-1
    list_counts = []
    list_phases = []

    data = dict.fromkeys(["counts","iseed"])

    for key in data.keys():
        data[key] = []

    if data_frame.empty:
        # empty data_frame, THIS SHOULD ME MODIFY MANUALY
        # this happens when EC3 is 0 and then, cck does not spike
        for i in range(iseed):
            data["counts"].append(np.array( np.zeros(nbins) ))
            data["iseed"].append([i]*nbins)

        data["counts"] = np.concatenate(data["counts"])
        data["iseed"]  = np.concatenate(data["iseed"])
        data = pd.DataFrame(data)

        data2 = dict.fromkeys(["binphase","counts_mean","counts_sem"])
        data2["binphase"] = binphase[1:]-np.diff(binphase)[0]/2
        data2["counts_mean"] = np.zeros(nbins)
        data2["counts_sem"]  = np.zeros(nbins)
        data2["kde"] = np.zeros(nbins)

    else: 

        for i in iseed:
            w = (data_frame["iseed"]==i)
            spikes = data_frame[w]["tvec"].values
            spikes = spikes[(spikes>2000) & (spikes<30000)]
            phases = np.mod( spikes-50, theta_ryhtm)*360/theta_ryhtm
            counts, bins = np.histogram( phases, bins=binphase, density="True" )
            list_counts.append(counts)
            list_phases.append(phases)
            # print(i,j, len(counts))
            data["counts"].append(counts)
            data["iseed"].append([i]*len(counts))
        
        if len(np.concatenate(list_phases)) == 0:
            data["counts"] = []
            for i in iseed:
                data["counts"].append(np.array( np.zeros(nbins) ))

            data["counts"] = np.concatenate(data["counts"])
            data["iseed"]  = np.concatenate(data["iseed"])
            data = pd.DataFrame(data)

            data2 = dict.fromkeys(["binphase","counts_mean","counts_sem"])
            data2["binphase"] = binphase[1:]-np.diff(binphase)[0]/2
            data2["counts_mean"] = np.zeros(nbins)
            data2["counts_sem"]  = np.zeros(nbins)
            data2["kde"] = np.zeros(nbins)
        
        else:
            data["counts"] = np.concatenate(data["counts"])
            data["iseed"]  = np.concatenate(data["iseed"])
            data = pd.DataFrame(data)
            
            data2 = dict.fromkeys(["binphase","counts_mean","counts_sem"])
            data2["binphase"] = bins[1:]-np.diff(bins)[0]/2
            data2["counts_mean"] = np.mean(list_counts, axis=0)
            data2["counts_sem"]  = np.std(list_counts, axis=0)/np.sqrt(len(list_counts))
            kde = stats.gaussian_kde(np.concatenate(list_phases))
            yfit = kde(data2["binphase"])
            data2["kde"] = yfit
        
    data2 = pd.DataFrame(data2)
    output = {}
    output["realizations"] = data
    output["average"] = data2
    
    return output

def compute_phases_distribution(data_frame, theta_ryhtm=125.0,binsize=15):
    iseed = np.unique(data_frame["iseed"])
    bseed = np.unique(data_frame["ibseed"])
    binphase = np.arange(0,361,15)
    nbins = len(binphase)-1
    list_counts = []
    list_phases = []

    data = dict.fromkeys(["counts","iseed","ibseed"])

    for key in data.keys():
        data[key] = []

    if data_frame.empty:
        # empty data_frame, THIS SHOULD ME MODIFY MANUALY
        # this happens when EC3 is 0 and then, cck does not spike
        for i in range(10):
            for j in range(10):
                data["counts"].append(np.array( np.zeros(nbins) ))
                data["iseed"].append([i]*nbins)
                data["ibseed"].append([j]*nbins)

        data["counts"] = np.concatenate(data["counts"])
        data["iseed"]  = np.concatenate(data["iseed"])
        data["ibseed"] = np.concatenate(data["ibseed"])
        data = pd.DataFrame(data)

        data2 = dict.fromkeys(["binphase","counts_mean","counts_sem"])
        data2["binphase"] = binphase[1:]-np.diff(binphase)[0]/2
        data2["counts_mean"] = np.zeros(nbins)
        data2["counts_sem"]  = np.zeros(nbins)
        data2["kde"] = np.zeros(nbins)

    else: 

        for i in iseed:
            for j in bseed: 
                w = (data_frame["iseed"]==i) & (data_frame["ibseed"]==j)
                spikes = data_frame[w]["tvec"].values
                spikes = spikes[(spikes>2000) & (spikes<30000)]
                phases = np.mod( spikes-50, theta_ryhtm)*360/theta_ryhtm
                counts, bins = np.histogram( phases, bins=binphase, density="True" )
                list_counts.append(counts)
                list_phases.append(phases)
                # print(i,j, len(counts))
                data["counts"].append(counts)
                data["iseed"].append([i]*len(counts))
                data["ibseed"].append([j]*len(counts))
        
        if len(np.concatenate(list_phases)) == 0:
            data["counts"] = []
            for i in iseed:
                for j in bseed:
                    data["counts"].append(np.array( np.zeros(nbins) ))

            data["counts"] = np.concatenate(data["counts"])
            data["iseed"]  = np.concatenate(data["iseed"])
            data["ibseed"] = np.concatenate(data["ibseed"]) 
            data = pd.DataFrame(data)

            data2 = dict.fromkeys(["binphase","counts_mean","counts_sem"])
            data2["binphase"] = binphase[1:]-np.diff(binphase)[0]/2
            data2["counts_mean"] = np.zeros(nbins)
            data2["counts_sem"]  = np.zeros(nbins)
            data2["kde"] = np.zeros(nbins)
        
        else:
            data["counts"] = np.concatenate(data["counts"])
            data["iseed"]  = np.concatenate(data["iseed"])
            data["ibseed"] = np.concatenate(data["ibseed"])
            data = pd.DataFrame(data)
            
            data2 = dict.fromkeys(["binphase","counts_mean","counts_sem"])
            data2["binphase"] = bins[1:]-np.diff(bins)[0]/2
            data2["counts_mean"] = np.mean(list_counts, axis=0)
            data2["counts_sem"]  = np.std(list_counts, axis=0)/np.sqrt(len(list_counts))
            kde = stats.gaussian_kde(np.concatenate(list_phases))
            yfit = kde(data2["binphase"])
            data2["kde"] = yfit
        
    data2 = pd.DataFrame(data2)
    output = {}
    output["realizations"] = data
    output["average"] = data2
    
    return output


def compute_phases_from_spikes_external_inputs(folder):
    theta_ryhtm = 125.0
    binsize = 15
    parent_dir = os.path.join(folder, "measurements")
    spikes_files = ["external_inputs_sep_180.lzma","external_inputs_sep_360.lzma","external_inputs_ec2_180.lzma",
                    "external_inputs_ec3_360.lzma","external_inputs_dg_regular.lzma","external_inputs_dg_burst.lzma"]
    
    aux = ["sep_180","sep_360","ec2_180","ec3_360","dg_regular","dg_burst"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem")
    data = dict.fromkeys(columns)

    outputs = []
    # sep 180
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[0]) )
    outputs.append( compute_phases_distribution_external_inputs(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # sep 360
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[1]) )
    outputs.append( compute_phases_distribution_external_inputs(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # ec2 180
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[2]) )
    outputs.append( compute_phases_distribution_external_inputs(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # ec3 360
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[3]) )
    outputs.append( compute_phases_distribution_external_inputs(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # dg regular
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[4]) )
    outputs.append( compute_phases_distribution_external_inputs(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    # dg burst
    data_frame = file_management.load_lzma( os.path.join(folder, spikes_files[5]) )
    outputs.append( compute_phases_distribution_external_inputs(data_frame, theta_ryhtm=theta_ryhtm,binsize=binsize) )

    for i, key in enumerate(aux):
        data[key+"_mean"] = outputs[i]["average"]["counts_mean"]
        data[key+"_sem"]  = outputs[i]["average"]["counts_sem"]
        data[key+"_kde"]  = outputs[i]["average"]["kde"]
    data["binphase"] = outputs[0]["average"]["binphase"]
    print(parent_dir)
    file_management.save_lzma(data,"phases_spikes_distribution_external_inputs.lzma", parent_dir=parent_dir)


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
    data = dict.fromkeys(columns) 

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
        data[key+"_mean"] = outputs[i]["average"]["counts_mean"]
        data[key+"_sem"] = outputs[i]["average"]["counts_sem"]
        data[key+"_kde"] = outputs[i]["average"]["kde"]
    data["binphase"] = outputs[0]["average"]["binphase"]

    file_management.save_lzma(data,"phases_spikes_distribution_ca3.lzma", parent_dir=parent_dir)


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
    data = dict.fromkeys(columns)

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
        data[key+"_mean"] = outputs[i]["average"]["counts_mean"]
        data[key+"_sem"]  = outputs[i]["average"]["counts_sem"]
        data[key+"_kde"]  = outputs[i]["average"]["kde"]
    data["binphase"] = outputs[0]["average"]["binphase"]
        
    file_management.save_lzma(data,"phases_spikes_distribution_ca1.lzma",parent_dir=parent_dir)