import numpy as np
import scipy.signal as signal
import pandas as pd
import scipy
import copy 
import sys
import os
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management
from signal_analysis import *

def compute_amplitude_coupling(data_frame,theta_reference,gamma_reference,fs = 1/0.002,f0_theta=8,df_theta=4.0, 
                               f0_gamma=60, df_gamma=40.0, norder=4, surrogate_test=True, nsurrogates=10):

    iseed = np.unique(data_frame["iseed"])
    ibseed = np.unique(data_frame["ibseed"])
    alpha = 0.05

    n = len(data_frame[theta_reference][(data_frame["iseed"]==0) & (data_frame["ibseed"]==0)])
    time = np.linspace(0, n/fs, n)
    wtime = (time>2) & (time<30)

    pbins = np.arange(0,2*np.pi,15*np.pi/180)
    nbins = len(pbins)

    amplitude_gamma_list = []
    amplitude_theta_list = []

    data = dict.fromkeys(["theta_amplitude","gamma_amplitude","iseed","ibseed"])
    for key in data.keys():
        data[key] = []

    np.random.seed(1234)
    amplitude_gamma_surrogates_list = []
    for i in iseed:
        for j in ibseed: 
            w = (data_frame["iseed"]==i) & (data_frame["ibseed"]==j) 

            sig1 = data_frame[w][theta_reference].values[wtime]
            sig2 = data_frame[w][gamma_reference].values[wtime]
            sig1 -= np.mean(sig1)
            sig2 -= np.mean(sig2)

            sig_theta, envelope_theta, phase_theta = bandpass_filter_and_hilbert_transform(sig1, fs=fs, f0=f0_theta, df=df_theta, norder=norder)
            sig_gamma, envelope_gamma, phase_gamma = bandpass_filter_and_hilbert_transform(sig2, fs=fs, f0=f0_gamma, df=df_gamma, norder=norder)

            amplitude1 = np.zeros(nbins-1)
            amplitude2 = np.zeros(nbins-1)
            for k in range(nbins-1):         
                pl = pbins[k]
                pr = pbins[k+1]                   
                indices=(phase_theta>=pl) & (phase_theta<pr)  
                amplitude1[k] = np.mean(envelope_gamma[indices])
                amplitude2[k] = np.mean(sig_theta[indices])

            bins = pbins[1:]-np.diff(pbins)[0]/2
            amplitude_gamma_list.append(amplitude1)
            amplitude_theta_list.append(amplitude2)
            data["theta_amplitude"].append(amplitude2)
            data["gamma_amplitude"].append(amplitude1)
            data["iseed"].append([i]*len(amplitude1))
            data["ibseed"].append([j]*len(amplitude1))

            # surrogate test, to test the amplitude of gamma is not random
            pi_surrogates = np.zeros((nsurrogates,nbins-1))
            if surrogate_test:
                surrogates = surrogate_generation(sig2, 10, method ="block boostrapping",nsplits=10)

                for l, surrogate in enumerate(surrogates): 
                    sig_gamma, envelope_gamma, phase_gamma = bandpass_filter_and_hilbert_transform(surrogate, fs=fs, f0=f0_gamma, df=df_gamma, norder=norder)
                    amplitude = np.zeros(len(pbins)-1)       
                    for k in range(len(pbins)-1):  
                        pl = pbins[k]
                        pr = pbins[k+1]                   
                        indices=(phase_theta>=pl) & (phase_theta<pr) 
                        amplitude[k] = np.mean(envelope_gamma[indices]) 
                    pi_surrogates[l] = amplitude
                amplitude_gamma_surrogates_list.append(pi_surrogates)
    
    data["theta_amplitude"] = np.concatenate(data["theta_amplitude"])
    data["gamma_amplitude"] = np.concatenate(data["gamma_amplitude"])
    data["iseed"] = np.concatenate(data["iseed"])
    data["ibseed"] = np.concatenate(data["ibseed"])
    data = pd.DataFrame(data)

    amplitude_gamma_mean = np.mean(amplitude_gamma_list,axis=0)
    amplitude_gamma_sem  = np.std(amplitude_gamma_list,axis=0)/np.sqrt(len(amplitude_gamma_list))
    amplitude_theta_mean = np.mean(amplitude_theta_list,axis=0)
    amplitude_theta_sem  = np.std(amplitude_theta_list,axis=0)/np.sqrt(len(amplitude_theta_list))

    amplitude_surrogates = np.concatenate(amplitude_gamma_surrogates_list,axis=0)

    p = np.zeros(np.shape(amplitude_surrogates)[-1])
    for i in range(np.shape(amplitude_surrogates)[-1]):
        # test signifcance 95& both sides
        mu, std = np.mean(amplitude_surrogates[:,i]), np.std(amplitude_surrogates[:,i])
        z = (amplitude_gamma_mean[i]-mu)/std
        z_critical = scipy.stats.norm.ppf(1-alpha/2)
        p[i] = (1-scipy.stats.norm.cdf(abs(z)))*2

    data2 = dict.fromkeys(["bins","theta_amplitude_mean","theta_amplitude_sem","gamma_amplitude_mean","gamma_amplitude_sem","pvalue"])
    data2["bins"] = bins
    data2["theta_amplitude_mean"] = amplitude_theta_mean
    data2["theta_amplitude_sem"]  = amplitude_theta_sem
    data2["gamma_amplitude_mean"] = amplitude_gamma_mean
    data2["gamma_amplitude_sem"]  = amplitude_gamma_sem 
    data2["pvalue"] = p
    data2 = pd.DataFrame(data2)

    output = {}
    output["realizations"] = data
    output["average"] = data2
    return output

def compute_amplitude_coupling_from_volt(folder, input_filename, output_filename):
    filename = os.path.join(folder, input_filename)
    data = file_management.load_lzma(filename)
    aux = ["Bdend","soma","Adend1","Adend2","Adend3"]

    columns = []
    for i,comp1 in enumerate(aux):
        for j, comp2 in enumerate(aux):
            columns.append(f"{comp1}_{comp2}")
    data_dict = dict.fromkeys(columns)

    outputs = []
    for i,comp1 in enumerate(aux):
        for j, comp2 in enumerate(aux):
            outputs.append(compute_amplitude_coupling(data,f"{comp1}_volt_mean",f"{comp2}_volt_mean"))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]

    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder, "measurements"))

def compute_amplitude_coupling_from_lfp(folder, input_filename, output_filename):
    filename = os.path.join(folder, input_filename)
    data = file_management.load_lzma(filename)
    # empty datafrmae
    data_new = pd.DataFrame()
    data_new = copy.deepcopy(data[data["electrode"]==-100])
    #remove column electrode
    data_new = data_new.drop(columns=["electrode"])
    # display(data_new)
    # rename column
    data_new = data_new.rename(columns={"lfp": "Bdend"})
    data_new["soma"]   = data[data["electrode"]==10]["lfp"].values
    data_new["Adend1"] = data[data["electrode"]==85]["lfp"].values
    data_new["Adend2"] = data[data["electrode"]==235]["lfp"].values
    data_new["Adend3"] = data[data["electrode"]==385]["lfp"].values

    aux = ["soma", "Bdend", "Adend1", "Adend2", "Adend3"]
    columns = []
    for i,comp1 in enumerate(aux):
        for j, comp2 in enumerate(aux):
            columns.append(f"{comp1}_{comp2}")
    data_dict = dict.fromkeys(columns)

    outputs = []
    for i,comp1 in enumerate(aux):
        for j, comp2 in enumerate(aux):
            outputs.append(compute_amplitude_coupling(data_new,comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]  
    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder, "measurements"))


def compute_amplitude_coupling_from_ica(folder, input_filename, output_filename):
    filename = os.path.join(folder, input_filename)
    data = file_management.load_lzma(filename)
    # empty datafrmae

    outputs = []
    data_new = pd.DataFrame()
    electrodes = [-100,10, 85, 235, 385]
    aux = ["soma", "Bdend", "Adend1", "Adend2", "Adend3"]
    data_new = copy.deepcopy(data[(data["electrode"]==-100) & (data["component"]=="Bdend")])
    data_new = data_new.drop(columns=["electrode", "component"])
    data_new = data_new.rename(columns={"ica": "Bdend_Bdend"})
    for elec, comp in zip(electrodes,aux):
        data_new[f"{comp}"] = data[(data["electrode"]==elec) & (data["component"]==comp)]["ica"].values #soma-soma, Bdend-Bdend, Adend1-Adend1, Adend2-Adend2, Adend3-Adend3

    columns = [] 
    for i,comp1 in enumerate(aux):
        for j, comp2 in enumerate(aux):
            columns.append(f"{comp1}_{comp2}")
    data_dict = dict.fromkeys(columns)

    for i,comp1 in enumerate(aux):
        for j, comp2 in enumerate(aux):
            outputs.append(compute_amplitude_coupling(data_new,comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]
    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder, "measurements"))


# two functions, one per ca3 and other per ca1 for not making too complicated functions
def compute_amplitude_coupling_from_synaptic_current_ca3(folder):
    ''' This function is different than the others: gamma-theta coupling it not computed across all posibilites, 
    instead, it is computed per each synaptic current and per the sum of synaptic current in the soma and dendrites.
    '''
    
    filename = os.path.join(folder, "synapses_pyr_ca3.lzma")
    data = file_management.load_lzma(filename)
    aux = [ "iAdend3GABA_olm",     "iAdend3AMPA_ec2360", "iAdend3NMDA_ec2360", "iAdend3GABA_noise", "iAdend3AMPA_noise",
            "iAdend3AMPA_ec2180",  "iAdend3NMDA_ec2180", "iAdend1AMPA_dgreg",  "iAdend1NMDA_dgreg", "iAdend1AMPA_dgburst",
            "iAdend1NMDA_dgburst", "isomaGABA_bas",      "isomaAMPA_noise",    "isomaGABA_noise",   "iBdendAMPA_pyr",
            "iBdendNMDA_pyr"]

    data_dict = {}
    outputs = []
    for i, comp in enumerate(aux):
        print(comp)
        outputs.append(compute_amplitude_coupling(data,f"{comp}_mean",f"{comp}_mean"))
        data_dict[f"{comp}"] = outputs[-1]["average"]

    # sum of currents per compartment
    data["iAdend3"] = data["iAdend3GABA_olm_mean"] + data["iAdend3AMPA_ec2360_mean"] + data["iAdend3NMDA_ec2360_mean"] + data["iAdend3GABA_noise_mean"] + data["iAdend3AMPA_noise_mean"] + data["iAdend3AMPA_ec2180_mean"] + data["iAdend3NMDA_ec2180_mean"]
    data["iAdend2"] = data["iAdend3GABA_olm_mean"]*0
    data["iAdend1"] = data["iAdend1AMPA_dgreg_mean"] + data["iAdend1NMDA_dgreg_mean"] + data["iAdend1AMPA_dgburst_mean"] + data["iAdend1NMDA_dgburst_mean"]
    data["isoma"]   = data["isomaGABA_bas_mean"]  + data["isomaAMPA_noise_mean"] + data["isomaGABA_noise_mean"]
    data["iBdend"]  = data["iBdendAMPA_pyr_mean"] + data["iBdendNMDA_pyr_mean"]
    
    for i, comp in enumerate(["iAdend3","iAdend2","iAdend1","isoma","iBdend"]):
        outputs.append(compute_amplitude_coupling(data,comp,comp))
        print(comp)
        data_dict[f"{comp}"] = outputs[-1]["average"]
    
    file_management.save_lzma(data_dict, "synapses_pyr_ca3_amplitude_coupling.lzma", parent_dir = os.path.join(folder, "measurements"))

def compute_amplitude_coupling_from_synaptic_current_ca1(folder):

    filename = os.path.join(folder, "synapses_pyr_ca1.lzma")
    data = file_management.load_lzma(filename)

    aux = ["isomaAMPA_noise",    "isomaGABA_noise",    "iAdend3AMPA_noise",  "iAdend3GABA_noise",
           "iAdend3AMPA_ec3180", "iAdend3NMDA_ec3180", "iAdend3AMPA_ec3360", "iAdend3NMDA_ec3360", 
           "iAdend1AMPA_pyrCA3", "iAdend1NMDA_pyrCA3", "isomaGABA_cck",      "iAdend2GABA_cck",  
           "isomaGABA_bas",      "iAdend3GABA_olm",    "iBdendAMPA_pyr",     "iBdendNMDA_pyr"]

    data_dict = {}
    outputs = []
    for i, comp in enumerate(aux):
        print(comp)
        outputs.append(compute_amplitude_coupling(data,f"{comp}_mean",f"{comp}_mean"))
        data_dict[f"{comp}"] = outputs[-1]["average"]
    # sum of currents per compartment
    # sum of currents per compartment
    data["iAdend3"] = data["iAdend3GABA_olm_mean"]+data["iAdend3AMPA_noise_mean"] + data["iAdend3GABA_noise_mean"] + data["iAdend3AMPA_ec3180_mean"] + data["iAdend3NMDA_ec3180_mean"] + data["iAdend3AMPA_ec3360_mean"] + data["iAdend3NMDA_ec3360_mean"]
    data["iAdend2"] = data["iAdend2GABA_cck_mean"]
    data["iAdend1"] = data["iAdend1AMPA_pyrCA3_mean"] + data["iAdend1NMDA_pyrCA3_mean"]
    data["isoma"]   = data["isomaGABA_cck_mean"]  + data["isomaGABA_bas_mean"] + data["isomaAMPA_noise_mean"] + data["isomaGABA_noise_mean"]
    data["iBdend"]  = data["iBdendAMPA_pyr_mean"] + data["iBdendNMDA_pyr_mean"]

    for i, comp in enumerate(["iAdend3","iAdend2","iAdend1","isoma","iBdend"]):
        print(comp)
        outputs.append(compute_amplitude_coupling(data,comp,comp))
        data_dict[f"{comp}"] = outputs[-1]["average"]
    
    file_management.save_lzma(data_dict, "synapses_pyr_ca1_amplitude_coupling.lzma", parent_dir = os.path.join(folder, "measurements"))
