###################################################################################################
# Computation of PSD of all the data possible
# Files: voltage: olm, bas, cck, pyr (all compartments) and pyr: Adend3-Bdend (the LFP of Ketamine's paper)
#        synaptic currents: all of them in separe, and the sum in each compartment: Bdend, Soma, Adend1, Adend2, Adend3
#        lfp: Bdend, soma, Adend1, Adend2, Adend3
#        ica: Bdend, soma, Adend1, Adend2, Adend3 x 5

# outputs: psd_voltage_ca1/3.lzma  / psd_voltage_traces.lzma (not saved by now)
#          psd_synaptic_ca1/3.lzma / psd_synaptic_traces.lzma (not saved by now)
#          psd_lfp_ca1/3.lzma      / psd_lfp_traces.lzma (not saved by now)
#          psd_ica_ca1/3.lzma      / psd_ica_traces.lzma (not saved by now)

# without normalize: but for comparison reasons, i.e., between voltage and lfp, psd should be normalized
# units of PSD unit**2/Hz 
########################################################################################################

import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
from scipy.signal import welch 
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management
from inputs_parameters import *

def compute_psd(data_frame, variable = "soma_volt_mean", fs=1/0.002, nperseg=1024, noverlap=512, nfft=2048):
    iseed = np.unique(data_frame["iseed"])
    ibseed = np.unique(data_frame["ibseed"])

    data = dict.fromkeys(["psd","iseed","ibseed"])
    for key in data.keys():
        data[key] = []
    data_ = dict.fromkeys(["theta_power","gamma_power","theta_freq","gamma_freq","iseed","ibseed"])
    for key in data_.keys():
        data_[key] = []

    n = len(data_frame[variable][(data_frame["iseed"]==0) & (data_frame["ibseed"]==0)])
    time = np.linspace(0, n/fs, n)
    wtime = (time>2) & (time<30)

    psd_list = []
    for i in iseed:
        for j in ibseed: 
            sig = data_frame[(data_frame["iseed"]==i) & (data_frame["ibseed"]==j)][variable].values[wtime]
            f, psd = welch(sig, fs=fs, nperseg=nperseg,noverlap=noverlap, nfft=nfft)
            data["psd"].append(psd)
            data["iseed"].append([i]*len(psd))
            data["ibseed"].append([j]*len(psd))
            psd_list.append(psd)

            df = np.diff(f)[0]
            wtheta = (f>4) & (f<12)
            wgamma = (f>30) & (f<100)
            data_["theta_power"].append( np.sum(psd[wtheta])*df)
            data_["gamma_power"].append( np.sum(psd[wgamma])*df)
            data_["theta_freq"].append( f[wtheta][np.argmax(psd[wtheta])])
            data_["gamma_freq"].append( f[wgamma][np.argmax(psd[wgamma])])
            data_["iseed"].append(i)
            data_["ibseed"].append(j)
    
    data["psd"] = np.concatenate(data["psd"])
    data["iseed"] = np.concatenate(data["iseed"])
    data["ibseed"] = np.concatenate(data["ibseed"])
    data = pd.DataFrame(data)

    data1 = pd.DataFrame(data_)

    data2 = dict.fromkeys(["f","psd_mean","psd_sem"]) 
    data2["f"] = f 
    data2["psd_mean"] = np.mean(psd_list, axis=0)
    data2["psd_sem"] = np.std(psd_list, axis=0)/np.sqrt(len(psd_list))
    data2 = pd.DataFrame(data2)

    output = {} 
    output["realizations"] = data
    output["measurements"] = data1
    output["average"] = data2
    return output

def compute_psd_volt_ca3(folder):
    parent_dir = os.path.join(folder, "measurements")
    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048
    
    aux = ["Bdend_mean", "soma",  "Adend1","Adend2","Adend3","bas","olm"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem")
    data_frame = dict.fromkeys(columns) 

    data_dict = dict.fromkeys(aux)

    data = file_management.load_lzma(os.path.join(folder, "volt_pyr_ca3.lzma"))
    outputs = []
    outputs.append( compute_psd(data, variable = "Bdend_volt_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) ) 
    outputs.append( compute_psd(data, variable = "soma_volt_mean",  fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    outputs.append( compute_psd(data, variable = "Adend1_volt_mean",fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    outputs.append( compute_psd(data, variable = "Adend2_volt_mean",fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    outputs.append( compute_psd(data, variable = "Adend3_volt_mean",fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )

    data = file_management.load_lzma(os.path.join(folder, "volt_bas_ca3.lzma"))
    outputs.append( compute_psd(data, variable = "soma_volt_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )

    data = file_management.load_lzma(os.path.join(folder, "volt_olm_ca3.lzma"))
    outputs.append( compute_psd(data, variable = "soma_volt_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )

    for i,key in enumerate(aux):
        data_frame[key+"_mean"] = outputs[i]["average"]["psd_mean"]
        data_frame[key+"_sem"] = outputs[i]["average"]["psd_sem"]
        data_dict[key] = outputs[i]["measurements"]

    data_frame["freq"] = outputs[0]["average"]["f"]

    #alternative lfp
    data = file_management.load_lzma(os.path.join(folder, "volt_pyr_ca3.lzma"))
    data["alternative"] = data["Adend3_volt_mean"].values - data["Bdend_volt_mean"].values
    outputs.append( compute_psd(data, variable = "alternative", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    data_frame["alternative_mean"] = outputs[-1]["average"]["psd_mean"]
    data_frame["alternative_sem"] = outputs[-1]["average"]["psd_sem"]
    data_dict["alternative"] = outputs[-1]["measurements"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_voltage_ca3.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict, " psd_voltage_measurements_ca3.lzma", parent_dir=parent_dir)


def compute_psd_volt_ca1(folder):
    parent_dir = os.path.join(folder, "measurements")
    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048
    
    aux = ["Bdend_mean", "soma",  "Adend1","Adend2","Adend3","bas","olm","cck"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem")
    data_frame = dict.fromkeys(columns)
    data_dict = dict.fromkeys(aux)

    outputs = []
    data = file_management.load_lzma(os.path.join(folder, "volt_pyr_ca1.lzma"))
    outputs.append( compute_psd(data, variable = "Bdend_volt_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    outputs.append( compute_psd(data, variable = "soma_volt_mean",  fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    outputs.append( compute_psd(data, variable = "Adend1_volt_mean",fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    outputs.append( compute_psd(data, variable = "Adend2_volt_mean",fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    outputs.append( compute_psd(data, variable = "Adend3_volt_mean",fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )

    data = file_management.load_lzma(os.path.join(folder, "volt_bas_ca1.lzma"))
    outputs.append( compute_psd(data, variable = "soma_volt_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )

    data = file_management.load_lzma(os.path.join(folder, "volt_olm_ca1.lzma"))
    outputs.append( compute_psd(data, variable = "soma_volt_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )

    data = file_management.load_lzma(os.path.join(folder, "volt_cck_ca1.lzma"))
    outputs.append( compute_psd(data, variable = "soma_volt_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )

    for i,key in enumerate(aux):
        data_frame[key+"_mean"] = outputs[i]["average"]["psd_mean"]
        data_frame[key+"_sem"] = outputs[i]["average"]["psd_sem"]
        data_dict[key] = outputs[i]["measurements"]
    data_frame["freq"] = outputs[0]["average"]["f"]

    # alternative lfp
    data = file_management.load_lzma(os.path.join(folder, "volt_pyr_ca1.lzma"))
    data["alternative"] = data["Adend3_volt_mean"].values - data["Bdend_volt_mean"].values
    outputs.append( compute_psd(data, variable = "alternative", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
    data_frame["alternative_mean"] = outputs[-1]["average"]["psd_mean"]
    data_frame["alternative_sem"] = outputs[-1]["average"]["psd_sem"]
    data_dict["alternative"] = outputs[-1]["measurements"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_voltage_ca1.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict,  "psd_voltage_measurements_ca1.lzma", parent_dir=parent_dir)

def compute_psd_synaptic_currents_ca3(folder):
    parent_dir = os.path.join(folder, "measurements")
    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048

    aux = [ "iAdend3GABA_olm",     "iAdend3AMPA_ec2360", "iAdend3NMDA_ec2360", "iAdend3GABA_noise", "iAdend3AMPA_noise",
            "iAdend3AMPA_ec2180",  "iAdend3NMDA_ec2180", "iAdend1AMPA_dgreg",  "iAdend1NMDA_dgreg", "iAdend1AMPA_dgburst",
            "iAdend1NMDA_dgburst", "isomaGABA_bas",      "isomaAMPA_noise",    "isomaGABA_noise",   "iBdendAMPA_pyr",
            "iBdendNMDA_pyr"]

    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem")
    data_frame = dict.fromkeys(columns)
    data_dict = dict.fromkeys(aux)

    data = file_management.load_lzma(os.path.join(folder, "synapses_pyr_ca3.lzma"))
    outputs = []
    for key in aux:
        outputs.append( compute_psd(data, variable = key+"_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
        data_frame[key+"_mean"] = outputs[-1]["average"]["psd_mean"]
        data_frame[key+"_sem"]  = outputs[-1]["average"]["psd_sem"]
        data_dict[key] = outputs[-1]["measurements"]
    data_frame["freq"] = outputs[0]["average"]["f"]

    # sum of currents per compartment
    data["iAdend3"] = data["iAdend3GABA_olm_mean"] + data["iAdend3AMPA_ec2360_mean"] + data["iAdend3NMDA_ec2360_mean"] + data["iAdend3GABA_noise_mean"] + data["iAdend3AMPA_noise_mean"] + data["iAdend3AMPA_ec2180_mean"] + data["iAdend3NMDA_ec2180_mean"]
    data["iAdend2"] = data["iAdend3GABA_olm_mean"]*0
    data["iAdend1"] = data["iAdend1AMPA_dgreg_mean"] + data["iAdend1NMDA_dgreg_mean"] + data["iAdend1AMPA_dgburst_mean"] + data["iAdend1NMDA_dgburst_mean"]
    data["isoma"]   = data["isomaGABA_bas_mean"]  + data["isomaAMPA_noise_mean"] + data["isomaGABA_noise_mean"]
    data["iBdend"]  = data["iBdendAMPA_pyr_mean"] + data["iBdendNMDA_pyr_mean"]

    for key in ["iAdend3","iAdend2","iAdend1","isoma","iBdend"]:
        outputs.append( compute_psd(data, variable = key, fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
        data_frame[key+"_mean"] = outputs[-1]["average"]["psd_mean"]
        data_frame[key+"_sem"]  = outputs[-1]["average"]["psd_sem"]
        data_dict[key] = outputs[-1]["measurements"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_synaptic_ca3.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict,  "psd_synaptic_measurements_ca3.lzma", parent_dir=parent_dir)

def compute_psd_synaptic_currents_ca1(folder):
    parent_dir = os.path.join(folder, "measurements")
    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048
    aux = ["isomaAMPA_noise",    "isomaGABA_noise",    "iAdend3AMPA_noise",  "iAdend3GABA_noise",
          "iAdend3AMPA_ec3180", "iAdend3NMDA_ec3180", "iAdend3AMPA_ec3360", "iAdend3NMDA_ec3360", 
          "iAdend1AMPA_pyrCA3", "iAdend1NMDA_pyrCA3", "isomaGABA_cck",      "iAdend2GABA_cck",  
          "isomaGABA_bas",      "iAdend3GABA_olm",    "iBdendAMPA_pyr",     "iBdendNMDA_pyr"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem")
    data_frame = dict.fromkeys(columns)
    data_dict = dict.fromkeys(aux)

    data = file_management.load_lzma(os.path.join(folder, "synapses_pyr_ca1.lzma"))
    outputs = []
    for key in aux:
        outputs.append( compute_psd(data, variable = key+"_mean", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
        data_frame[key+"_mean"] = outputs[-1]["average"]["psd_mean"]
        data_frame[key+"_sem"]  = outputs[-1]["average"]["psd_sem"]
        data_dict[key] = outputs[-1]["measurements"]
    data_frame["freq"] = outputs[0]["average"]["f"]

    # sum of currents per compartment
    data["iAdend3"] = data["iAdend3GABA_olm_mean"]+data["iAdend3AMPA_noise_mean"] + data["iAdend3GABA_noise_mean"] + data["iAdend3AMPA_ec3180_mean"] + data["iAdend3NMDA_ec3180_mean"] + data["iAdend3AMPA_ec3360_mean"] + data["iAdend3NMDA_ec3360_mean"]
    data["iAdend2"] = data["iAdend2GABA_cck_mean"]
    data["iAdend1"] = data["iAdend1AMPA_pyrCA3_mean"] + data["iAdend1NMDA_pyrCA3_mean"]
    data["isoma"]   = data["isomaGABA_cck_mean"]  + data["isomaGABA_bas_mean"] + data["isomaAMPA_noise_mean"] + data["isomaGABA_noise_mean"]
    data["iBdend"]  = data["iBdendAMPA_pyr_mean"] + data["iBdendNMDA_pyr_mean"]

    for key in ["iAdend3","iAdend2","iAdend1","isoma","iBdend"]:
        outputs.append( compute_psd(data, variable = key, fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
        data_frame[key+"_mean"] = outputs[-1]["average"]["psd_mean"]
        data_frame[key+"_sem"]  = outputs[-1]["average"]["psd_sem"]
        data_dict[key] = outputs[-1]["measurements"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_synaptic_ca1.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict,  "psd_synaptic_measurements_ca1.lzma", parent_dir=parent_dir)

def compute_psd_lfp_ca3(folder):
    parent_dir = os.path.join(folder, "measurements")
    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048
    aux = ["Bdend", "soma",  "Adend1","Adend2","Adend3"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem")
    data_frame = dict.fromkeys(columns)
    data_dict = dict.fromkeys(aux)
    electrode = [-100,   10,   85,  235,  385]

    data = file_management.load_lzma(os.path.join(folder, "lfp_ca3.lzma"))
    outputs = []
    for key, elec in zip(aux,electrode):
        data_ = data[data["electrode"] == elec]
        outputs.append( compute_psd(data_, variable = "lfp", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
        data_frame[key+"_mean"] = outputs[-1]["average"]["psd_mean"]
        data_frame[key+"_sem"]  = outputs[-1]["average"]["psd_sem"]
        data_dict[key] = outputs[-1]["measurements"]

    data_frame["freq"] = outputs[0]["average"]["f"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_lfp_ca3.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict,  "psd_lfp_measurements_ca3.lzma", parent_dir=parent_dir)

def compute_psd_lfp_ca1(folder):
    parent_dir = os.path.join(folder, "measurements")
    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048    
    aux = ["Bdend", "soma",  "Adend1","Adend2","Adend3"]
    columns = []
    for key in aux:
        columns.append(key+"_mean")
        columns.append(key+"_sem")
    data_frame = dict.fromkeys(columns)
    data_dict = dict.fromkeys(aux)
    electrode = [-100,   10,   85,  235, 385]

    data = file_management.load_lzma(os.path.join(folder, "lfp_ca1.lzma"))
    outputs = []
    for key, elec in zip(aux,electrode):
        data_ = data[data["electrode"] == elec]
        outputs.append( compute_psd(data_, variable = "lfp", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
        data_frame[key+"_mean"] = outputs[-1]["average"]["psd_mean"]
        data_frame[key+"_sem"]  = outputs[-1]["average"]["psd_sem"]
        data_dict[key] = outputs[-1]["measurements"]
    data_frame["freq"] = outputs[0]["average"]["f"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_lfp_ca1.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict,  "psd_lfp_measurements_ca1.lzma", parent_dir=parent_dir)

def compute_psd_ica_ca3(folder):
    parent_dir = os.path.join(folder, "measurements")

    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048

    aux = ["Bdend", "soma",  "Adend1","Adend2","Adend3"]
    columns = []
    columns_ = []
    for key1 in aux:
        for key2 in aux:
            columns.append(key1+"_"+key2+"_mean")
            columns.append(key1+"_"+key2+"_sem")
            columns_.append(key1+"_"+key2)
    data_dict = dict.fromkeys(columns_)
    data_frame = dict.fromkeys(columns)
    data = file_management.load_lzma(os.path.join(folder, "ica_ca3.lzma"))
    electrode = [-100,   10,   85,  235,  385]

    outputs = []
    for key1, elec in zip(aux, electrode):
        for key2 in aux:
            w = (data["electrode"]==elec) & (data["component"]==key2)
            data_ = data[w]
            outputs.append( compute_psd(data_, variable = "ica", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
            data_frame[key1+"_"+key2+"_mean"] = outputs[-1]["average"]["psd_mean"]
            data_frame[key1+"_"+key2+"_sem"]  = outputs[-1]["average"]["psd_sem"]
            data_dict[key1+"_"+key2] = outputs[-1]["measurements"]
    data_frame["freq"] = outputs[0]["average"]["f"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_ica_ca3.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict,  "psd_ica_measurements_ca3.lzma", parent_dir=parent_dir)


def compute_psd_ica_ca1(folder):
    parent_dir = os.path.join(folder, "measurements")
    fs = 1/0.002
    nperseg = 1024
    noverlap = 512
    nfft = 2048
    aux = ["Bdend", "soma",  "Adend1","Adend2","Adend3"]
    columns, columns_ = [],[]
    for key1 in aux:
        for key2 in aux:
            columns.append(key1+"_"+key2+"_mean")
            columns.append(key1+"_"+key2+"_sem")
            columns_.append(key1+"_"+key2)
    data_dict = dict.fromkeys(columns_)
    data_frame = dict.fromkeys(columns)
    data = file_management.load_lzma(os.path.join(folder, "ica_ca1.lzma"))
    electrode = [-100,   10,   85,  235,  385]

    outputs = []
    for key1, elec in zip(aux, electrode):
        for key2 in aux:
            w = (data["electrode"]==elec) & (data["component"]==key2)
            data_ = data[w]
            outputs.append( compute_psd(data_, variable = "ica", fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft) )
            data_frame[key1+"_"+key2+"_mean"] = outputs[-1]["average"]["psd_mean"]
            data_frame[key1+"_"+key2+"_sem"]  = outputs[-1]["average"]["psd_sem"]
            data_dict[key1+"_"+key2] = outputs[-1]["measurements"]
    data_frame["freq"] = outputs[0]["average"]["f"]

    data_frame = pd.DataFrame(data_frame)
    file_management.save_lzma(data_frame, "psd_ica_ca1.lzma", parent_dir=parent_dir)
    file_management.save_lzma(data_dict,  "psd_ica_measurements_ca1.lzma", parent_dir=parent_dir)
