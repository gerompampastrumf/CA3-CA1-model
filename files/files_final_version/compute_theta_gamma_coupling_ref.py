'''
As a diference to the other file called simililary, the idea of this function is to compute the gamma-theta coupling across sections
but keeping the phase reference of the soma pyr CA1, which is an standard. However this can be changed inside the functions.
'''

import numpy as np
import scipy.signal as signal
import pandas as pd
import scipy
from scipy import stats
import copy 
import sys
import os
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management
from signal_analysis import *

def compute_amplitude_coupling(data_frame,theta_reference,theta_component,gamma_component,fs = 1/0.002,f0_theta=8,df_theta=4.0, 
                               f0_gamma=60, df_gamma=40.0, norder=4,surrogate_test=True, nsurrogates=10):

    iseed = np.unique(data_frame["iseed"])
    ibseed = np.unique(data_frame["ibseed"])
    alpha = 0.05

    n = len(data_frame[theta_reference][(data_frame["iseed"]==0) & (data_frame["ibseed"]==0)])
    time = np.linspace(0, n/fs, n)
    wtime = (time>2) & (time<30)

    pbins = np.arange(0,2*np.pi+1e-3,5*np.pi/180) # binsize increase x3, from 15º to 5º in order to visualize thet-phase difference between soma and Adend3
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
            
            sig0 = data_frame[w][theta_reference].values[wtime]
            sig1 = data_frame[w][theta_component].values[wtime]
            sig2 = data_frame[w][gamma_component].values[wtime]
            sig0 -= np.mean(sig0)
            sig1 -= np.mean(sig1)
            sig2 -= np.mean(sig2)

            sig_ref, envelope_ref, phase_ref = bandpass_filter_and_hilbert_transform(sig0, fs=fs, f0=f0_theta, df=df_theta, norder=norder)

            sig_theta, envelope_theta, phase_theta = bandpass_filter_and_hilbert_transform(sig1, fs=fs, f0=f0_theta, df=df_theta, norder=norder)
            sig_gamma, envelope_gamma, phase_gamma = bandpass_filter_and_hilbert_transform(sig2, fs=fs, f0=f0_gamma, df=df_gamma, norder=norder)

            amplitude1 = np.zeros(nbins-1)
            amplitude2 = np.zeros(nbins-1)
            for k in range(nbins-1):         
                pl = pbins[k]
                pr = pbins[k+1]                   
                indices=(phase_ref>=pl) & (phase_ref<pr)  
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
                        indices=(phase_ref>=pl) & (phase_ref<pr) 
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
        # test signifcance 95& both sidesº
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

def compute_amplitude_coupling_from_volt(folder_ca3,folder_ca1,input_filename_ca3,input_filename_ca1, output_filename):
    filename_ca3 = os.path.join(folder_ca3, input_filename_ca3)
    filename_ca1 = os.path.join(folder_ca1, input_filename_ca1)
    data_ca3 = file_management.load_lzma(filename_ca3)
    data_ca1 = file_management.load_lzma(filename_ca1)

    data_new = pd.DataFrame()
    data_new["iseed"] = data_ca3["iseed"].values
    data_new["ibseed"] = data_ca3["ibseed"].values
    data_new["Bdend_ca3_volt_mean"]  = data_ca3["Bdend_volt_mean"].values
    data_new["soma_ca3_volt_mean"]   = data_ca3["soma_volt_mean"].values
    data_new["Adend1_ca3_volt_mean"] = data_ca3["Adend1_volt_mean"].values
    data_new["Adend2_ca3_volt_mean"] = data_ca3["Adend2_volt_mean"].values
    data_new["Adend3_ca3_volt_mean"] = data_ca3["Adend3_volt_mean"].values

    data_new["Bdend_ca1_volt_mean"]  = data_ca1["Bdend_volt_mean"].values
    data_new["soma_ca1_volt_mean"]   = data_ca1["soma_volt_mean"].values
    data_new["Adend1_ca1_volt_mean"] = data_ca1["Adend1_volt_mean"].values
    data_new["Adend2_ca1_volt_mean"] = data_ca1["Adend2_volt_mean"].values
    data_new["Adend3_ca1_volt_mean"] = data_ca1["Adend3_volt_mean"].values

    aux_ca3 = ["Bdend_ca3","soma_ca3", "Adend1_ca3", "Adend2_ca3", "Adend3_ca3"]
    aux_ca1 = ["Bdend_ca1","soma_ca1", "Adend1_ca1", "Adend2_ca1", "Adend3_ca1"]
    columns = []
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            columns.append(f"{comp1}_{comp2}")
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            columns.append(f"{comp1}_{comp2}")
    data_dict = dict.fromkeys(columns)

    outputs = []
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            outputs.append(compute_amplitude_coupling(data_new,"soma_ca1_volt_mean",f"{comp1}_volt_mean",f"{comp2}_volt_mean"))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            outputs.append(compute_amplitude_coupling(data_new,"soma_ca1_volt_mean",f"{comp1}_volt_mean",f"{comp2}_volt_mean"))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]

    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder_ca1, "measurements"))


def compute_amplitude_coupling_from_lfp(folder_ca3, folder_ca1, input_filename_ca3, input_filename_ca1, output_filename):
    # by default save in the ca1 folder 

    filename_ca3 = os.path.join(folder_ca3, input_filename_ca3)
    filename_ca1 = os.path.join(folder_ca1, input_filename_ca1)

    data_ca3 = file_management.load_lzma(filename_ca3)
    data_ca1 = file_management.load_lzma(filename_ca1)

    # empty datafrmae
    data_new = pd.DataFrame()
    data_new = copy.deepcopy(data_ca3[data_ca3["electrode"]==-100])
    #remove column electrode
    data_new = data_new.drop(columns=["electrode"])
    # display(data_new)
    # rename column
    data_new = data_new.rename(columns={"lfp": "Bdend_ca3"})
    data_new["soma_ca3"]   = data_ca3[data_ca3["electrode"]==10]["lfp"].values
    data_new["Adend1_ca3"] = data_ca3[data_ca3["electrode"]==85]["lfp"].values
    data_new["Adend2_ca3"] = data_ca3[data_ca3["electrode"]==235]["lfp"].values
    data_new["Adend3_ca3"] = data_ca3[data_ca3["electrode"]==385]["lfp"].values

    data_new["Bdend_ca1"] = data_ca1[data_ca1["electrode"]==-100]["lfp"].values
    data_new["soma_ca1"]   = data_ca1[data_ca1["electrode"]==10]["lfp"].values
    data_new["Adend1_ca1"] = data_ca1[data_ca1["electrode"]==85]["lfp"].values
    data_new["Adend2_ca1"] = data_ca1[data_ca1["electrode"]==235]["lfp"].values
    data_new["Adend3_ca1"] = data_ca1[data_ca1["electrode"]==385]["lfp"].values

    aux_ca3 = ["Bdend_ca3","soma_ca3", "Adend1_ca3", "Adend2_ca3", "Adend3_ca3"]
    aux_ca1 = ["Bdend_ca1","soma_ca1", "Adend1_ca1", "Adend2_ca1", "Adend3_ca1"]

    columns = []
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            columns.append(f"{comp1}_{comp2}")
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            columns.append(f"{comp1}_{comp2}")
    data_dict = dict.fromkeys(columns)

    outputs = []
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            outputs.append(compute_amplitude_coupling(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"] 
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            outputs.append(compute_amplitude_coupling(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]

    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder_ca1, "measurements"))


def compute_amplitude_coupling_from_ica(folder_ca3,folder_ca1, input_filename_ca3, input_filename_ca1, output_filename):
    input_filename_ca3  = os.path.join(folder_ca3, input_filename_ca3)
    input_filename_ca1  = os.path.join(folder_ca1, input_filename_ca1)
    data_ca3 = file_management.load_lzma(input_filename_ca3)
    data_ca1 = file_management.load_lzma(input_filename_ca1)

    outputs = []
    electrodes = [-100,10, 85, 235, 385]
    aux = ["Bdend", "soma", "Adend1", "Adend2", "Adend3"]
    aux_ca3 = ["Bdend_ca3","soma_ca3","Adend1_ca3", "Adend2_ca3", "Adend3_ca3"]
    aux_ca1 = ["Bdend_ca1","soma_ca1","Adend1_ca1", "Adend2_ca1", "Adend3_ca1"]
  
    data_new = {}
    for comp in aux:
        data_new[f"{comp}_ca3"] = []
        data_new[f"{comp}_ca1"] = []
    data_new["iseed"] = []
    data_new["ibseed"] = [] 
    
    iseed = np.sort(np.unique(data_ca3["iseed"].values))
    ibseed = np.sort(np.unique(data_ca3["ibseed"].values))

    try:
        ilist, jlist = [],[]
        for i in iseed:
            for j in ibseed:
                wca3 = (data_ca3["iseed"]==i) & (data_ca3["ibseed"]==j)
                wca1 = (data_ca1["iseed"]==i) & (data_ca1["ibseed"]==j)
                for elec, comp in zip(electrodes, aux):
                    w = wca3 & (data_ca3["electrode"]==elec) & (data_ca3["component"]==comp)
                    data_new[f"{comp}_ca3"].append( data_ca3[w]["ica"].values )
                    w = wca1 & (data_ca1["electrode"]==elec) & (data_ca1["component"]==comp)
                    data_new[f"{comp}_ca1"].append( data_ca1[w]["ica"].values )
                n = len(data_new[f"{comp}_ca3"][-1])
                ilist.append([i]*n)
                jlist.append([j]*n)
        for comp in aux: 
            data_new[f"{comp}_ca3"] = np.concatenate(data_new[f"{comp}_ca3"])
            data_new[f"{comp}_ca1"] = np.concatenate(data_new[f"{comp}_ca1"])
        data_new["iseed"] = np.concatenate(ilist)
        data_new["ibseed"] = np.concatenate(jlist)
        data_new = pd.DataFrame(data_new)
    except:
        print("Probably the data frames have different iseed and ibseed columns")
      
    columns = [] 
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            columns.append(f"{comp1}_{comp2}")
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            columns.append(f"{comp1}_{comp2}")
    data_dict = dict.fromkeys(columns)

    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            outputs.append(compute_amplitude_coupling(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            outputs.append(compute_amplitude_coupling(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]
    
    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder_ca1, "measurements"))




