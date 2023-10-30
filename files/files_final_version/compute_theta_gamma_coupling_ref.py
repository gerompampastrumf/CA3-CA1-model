import numpy as np
import scipy.signal as signal
import pandas as pd
import scipy
from scipy import stats
import copy 
import sys
import os
from scipy.interpolate import interp1d
from neurodsp.filt import filter_signal
from bycycle.cyclepoints import find_extrema, find_zerox
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management
from signal_analysis import *

def compute_amplitude_coupling(data_frame,theta_reference,theta_component,gamma_component,nbins=40,fs = 1/0.002,f0_theta=8,df_theta=2.0,f0_gamma=60, df_gamma=40.0, norder=4,surrogate_test=True, nsurrogates=10):
    iseed = np.unique(data_frame["iseed"])
    ibseed = np.unique(data_frame["ibseed"])
    alpha = 0.05

    n = len(data_frame[theta_reference][(data_frame["iseed"]==0) & (data_frame["ibseed"]==0)])
    time = np.linspace(0, n/fs, n)
    wtime = (time>2) & (time<30)
    
    pbins = np.linspace(0,2*np.pi,nbins+1)
    bins = pbins[1:]-np.diff(pbins)[0]/2
    # pbins = np.arange(0,2*np.pi+1e-3,5*np.pi/180) # igual esto hay que aumentarlo, 
    # nbins = len(pbins)

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
            sig0 -= np.nanmean(sig0)
            sig1 -= np.nanmean(sig1)
            sig2 -= np.nanmean(sig2)

            sig_ref, envelope_ref, phase_ref = bandpass_filter_and_hilbert_transform(sig0, fs=fs, f0=f0_theta, df=df_theta, norder=norder)
            sig_theta, envelope_theta, phase_theta = bandpass_filter_and_hilbert_transform(sig1, fs=fs, f0=f0_theta, df=df_theta, norder=norder)
            sig_gamma, envelope_gamma, phase_gamma = bandpass_filter_and_hilbert_transform(sig2, fs=fs, f0=f0_gamma, df=df_gamma, norder=norder)

            phase_ref = np.mod(phase_ref+np.pi,2*np.pi)
            phase_theta = np.mod(phase_theta+np.pi,2*np.pi)

            phase_diff = np.mod(phase_theta-phase_ref,2*np.pi) 
            phase_diff[phase_diff>np.pi] -= 2*np.pi
            counts, phase_bins = np.histogram(phase_diff, bins=21)
            phase_bins = phase_bins[1:]-np.diff(phase_bins)[0]/2
            phase_diff = phase_bins[np.argmax(counts)]
                
            phase_theta_adjusted = np.mod(phase_theta-phase_diff,2*np.pi)

            amplitude1 = np.zeros(nbins)
            amplitude2 = np.zeros(nbins)
            for k in range(nbins):         
                pl = pbins[k]
                pr = pbins[k+1]                   
                indices=(phase_theta_adjusted>=pl) & (phase_theta_adjusted<pr)  
                amplitude1[k] = np.nanmean(envelope_gamma[indices])
                amplitude2[k] = np.nanmean(sig_theta[indices])

            amplitude_gamma_list.append(amplitude1)
            amplitude_theta_list.append(amplitude2)
            data["theta_amplitude"].append(amplitude2)
            data["gamma_amplitude"].append(amplitude1)
            data["iseed"].append([i]*len(amplitude1))
            data["ibseed"].append([j]*len(amplitude1))

            # surrogate test, to test the amplitude of gamma is not random
            pi_surrogates = np.zeros((nsurrogates,nbins))
            if surrogate_test:
                surrogates = surrogate_generation(sig2, 10, method ="block boostrapping",nsplits=10)

                for l, surrogate in enumerate(surrogates): 
                    sig_gamma, envelope_gamma, phase_gamma = bandpass_filter_and_hilbert_transform(surrogate, fs=fs, f0=f0_gamma, df=df_gamma, norder=norder)
                    amplitude = np.zeros(nbins)       
                    for k in range(nbins):  
                        pl = pbins[k]
                        pr = pbins[k+1]                   
                        indices=(phase_theta_adjusted>=pl) & (phase_theta_adjusted<pr) 
                        amplitude[k] = np.nanmean(envelope_gamma[indices]) 
                    pi_surrogates[l] = amplitude
                amplitude_gamma_surrogates_list.append(pi_surrogates)
    
    data["theta_amplitude"] = np.concatenate(data["theta_amplitude"])
    data["gamma_amplitude"] = np.concatenate(data["gamma_amplitude"])
    data["iseed"] = np.concatenate(data["iseed"])
    data["ibseed"] = np.concatenate(data["ibseed"])
    data = pd.DataFrame(data)

    amplitude_gamma_mean = np.nanmean(amplitude_gamma_list,axis=0)
    amplitude_gamma_sem  = np.nanstd(amplitude_gamma_list,axis=0)/np.sqrt(len(amplitude_gamma_list))
    amplitude_theta_mean = np.nanmean(amplitude_theta_list,axis=0)
    amplitude_theta_sem  = np.nanstd(amplitude_theta_list,axis=0)/np.sqrt(len(amplitude_theta_list))

    amplitude_surrogates = np.concatenate(amplitude_gamma_surrogates_list,axis=0)

    p = np.zeros(np.shape(amplitude_surrogates)[-1])
    for i in range(np.shape(amplitude_surrogates)[-1]):
        # test signifcance 95& both sidesº
        mu, std = np.nanmean(amplitude_surrogates[:,i]), np.nanstd(amplitude_surrogates[:,i])
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

def compute_amplitude_coupling_bycycle(data_frame,theta_reference,theta_component,gamma_component,nbins=40,fs = 1/0.002,f0_theta=8,df_theta=2.0, 
                               f0_gamma=60, df_gamma=40.0, norder=4,surrogate_test=True, nsurrogates=10):
    ''' Alternative method where phase conversion of the theta component is done in the temporal domain, using the bycycle package.
    Here the bins are non-constant,'''
    
    f_theta = (f0_theta-df_theta, f0_theta+df_theta)
    f_lowpass = 20
    n_seconds_filter = .5
    n_seconds_theta = .75

    iseed = np.unique(data_frame["iseed"])
    ibseed = np.unique(data_frame["ibseed"])
    alpha = 0.05

    n = len(data_frame[theta_reference][(data_frame["iseed"]==0) & (data_frame["ibseed"]==0)])
    time = np.linspace(0, n/fs, n)
    wtime = (time>2) & (time<30)
    time = time[wtime]
    
    data = dict.fromkeys(["theta_amplitude","gamma_amplitude","iseed","ibseed"])
    for key in data.keys():
        data[key] = []

    pbins = np.linspace(0, 2*np.pi, nbins+1)
    bins = pbins[1:]-np.diff(pbins)[0]/2

    np.random.seed(1234)
    amplitude_gamma_list = []
    amplitude_theta_list = []
    amplitude_gamma_surrogates_list = []
    for i in iseed:
        for j in ibseed: 
            w = (data_frame["iseed"]==i) & (data_frame["ibseed"]==j)    
            
            sig0 = data_frame[w][theta_reference].values[wtime]
            sig1 = data_frame[w][theta_component].values[wtime]
            sig2 = data_frame[w][gamma_component].values[wtime]
            sig0 -= np.nanmean(sig0)
            sig1 -= np.nanmean(sig1)
            sig2 -= np.nanmean(sig2)

            # band pass filtering gamma reference
            # sig_theta, envelope_theta, phase_theta = bandpass_filter_and_hilbert_transform(sig1, fs=fs, f0=f0_theta, df=df_theta, norder=norder)
            sig_gamma, envelope_gamma, phase_gamma = bandpass_filter_and_hilbert_transform(sig2, fs=fs, f0=f0_gamma, df=df_gamma, norder=norder)
            ###############################################
            sig0_low = filter_signal(sig0, fs, 'lowpass', f_lowpass, n_seconds=n_seconds_filter, remove_edges=False)
            sig0_narrow = filter_signal(sig0, fs, 'bandpass', f_theta, n_seconds=n_seconds_theta, remove_edges=False)
            peaks, troughs = find_extrema(sig0_low, fs, f_theta, filter_kwargs={'n_seconds':n_seconds_theta})
            rises, decays  = find_zerox(sig0_narrow, peaks, troughs)

            troughs_min = np.min(troughs)
            troughs_max = np.max(troughs)

            while rises[0] < troughs_min:
                rises = rises[1:]
            while rises[-1] > troughs_max:
                rises = rises[:-1]

            while peaks[0] < troughs_min:
                peaks = peaks[1:]
            while peaks[-1] > troughs_max:
                peaks = peaks[:-1]

            while decays[0] <troughs_min:
                decays = decays[1:]
            while decays[-1] > troughs_max:
                decays = decays[:-1]

            time_ref_binned,phase_theta_ref_binned = [],[]
            for ii in range(len(troughs)-1):
                time_ref_binned.append(np.linspace(time[troughs[ii]], time[rises[ii]],     nbins//4, endpoint=False))
                time_ref_binned.append(np.linspace(time[rises[ii]],   time[peaks[ii]],     nbins//4, endpoint=False))
                time_ref_binned.append(np.linspace(time[peaks[ii]],   time[decays[ii]],    nbins//4, endpoint=False))
                time_ref_binned.append(np.linspace(time[decays[ii]],  time[troughs[ii+1]], nbins//4, endpoint=False))
                phase_theta_ref_binned.append(np.linspace(0, 2*np.pi, nbins, endpoint=False))

            time_ref_binned = np.concatenate(time_ref_binned)
            phase_theta_ref_binned = np.concatenate(phase_theta_ref_binned)
    
            ####################################################################################
            sig1_low = filter_signal(sig1, fs, 'lowpass', f_lowpass, n_seconds=n_seconds_filter, remove_edges=False)
            sig1_narrow = filter_signal(sig1, fs, 'bandpass', f_theta, n_seconds=n_seconds_theta, remove_edges=False)
            peaks, troughs = find_extrema(sig1_low, fs, f_theta, filter_kwargs={'n_seconds':n_seconds_theta})
            rises, decays  = find_zerox(sig1_narrow, peaks, troughs)
            
            troughs_min = np.min(troughs)
            troughs_max = np.max(troughs)

            while rises[0] < troughs_min:
                rises = rises[1:]
            while rises[-1] > troughs_max:
                rises = rises[:-1]

            while peaks[0] < troughs_min:
                peaks = peaks[1:]
            while peaks[-1] > troughs_max:
                peaks = peaks[:-1]

            while decays[0] <troughs_min:
                decays = decays[1:]
            while decays[-1] > troughs_max:
                decays = decays[:-1]

            time_sig1_binned,phase_theta_sig1_binned = [],[]
            for ii in range(len(troughs)-1):
                time_sig1_binned.append(np.linspace(time[troughs[ii]], time[rises[ii]],     nbins//4, endpoint=False))
                time_sig1_binned.append(np.linspace(time[rises[ii]],   time[peaks[ii]],     nbins//4, endpoint=False))
                time_sig1_binned.append(np.linspace(time[peaks[ii]],   time[decays[ii]],    nbins//4, endpoint=False))
                time_sig1_binned.append(np.linspace(time[decays[ii]],  time[troughs[ii+1]], nbins//4, endpoint=False))
                phase_theta_sig1_binned.append(np.linspace(0, 2*np.pi, nbins, endpoint=False))

            time_sig1_binned = np.concatenate(time_sig1_binned)
            phase_theta_sig1_binned = np.concatenate(phase_theta_sig1_binned)

            # adjusting the phases 
            interpolator = interp1d(time_ref_binned, phase_theta_ref_binned, kind='linear', fill_value='extrapolate')
            interpolated_phase_theta_ref_binned = interpolator(time_sig1_binned)

            phase_diff = np.mod(phase_theta_sig1_binned-interpolated_phase_theta_ref_binned,2*np.pi)
            phase_diff[phase_diff>np.pi] -= 2*np.pi
            counts, phase_bins = np.histogram(phase_diff, bins=21)
            phase_bins = phase_bins[1:]-np.diff(phase_bins)[0]/2
            phase_diff = phase_bins[np.argmax(counts)]

            phase_theta_sig1_adjusted = np.mod(phase_theta_sig1_binned-phase_diff,2*np.pi)
            # band pass filtering theta component

            # binning the envelope 
            envelope_gamma_binned = np.ones(len(time_sig1_binned))*np.nan
            sig_theta_binned = np.ones(len(time_sig1_binned))*np.nan
            for ii in range(len(time_sig1_binned)):
                if ii == len(time_sig1_binned)-1:
                    bin_start = time_sig1_binned[-1]
                    bin_end = time_sig1_binned[-1]+0.002
                else:
                    bin_start = time_sig1_binned[ii]
                    bin_end = time_sig1_binned[ii+1]
                
                # Extract values that fall within the current bin
                sig_theta_binned[ii] = np.nan if not np.any(np.where((bin_start <= time) & (time < bin_end))) else np.nanmean(sig1_low[(bin_start <= time) & (time < bin_end)])
                envelope_gamma_binned[ii] = np.nan if not np.any(np.where((bin_start <= time) & (time < bin_end))) else np.nanmean(envelope_gamma[(bin_start <= time) & (time < bin_end)])
                # envelope_gamma_binned[ii] = np.nanmean(envelope_gamma[(bin_start <= time) & (time < bin_end)])
                # sig_theta_binned[ii] = np.nanmean(sig1_low[(bin_start <= time) & (time < bin_end)])
            
            amplitude1 = np.zeros(nbins)
            amplitude2 = np.zeros(nbins)
            for k in range(nbins):
                pl = pbins[k]
                pr = pbins[k+1]                   
                indices=(phase_theta_sig1_adjusted>=pl) & (phase_theta_sig1_adjusted<pr)  
                amplitude1[k] = np.nanmean(envelope_gamma_binned[indices])
                amplitude2[k] = np.nanmean(sig_theta_binned[indices])
                
            amplitude_gamma_list.append(amplitude1)
            amplitude_theta_list.append(amplitude2)
            data["theta_amplitude"].append(amplitude2)
            data["gamma_amplitude"].append(amplitude1)
            data["iseed"].append([i]*len(amplitude1))
            data["ibseed"].append([j]*len(amplitude1))

            # surrogate test, to test the amplitude of gamma is not random
            pi_surrogates = np.zeros((nsurrogates,nbins))
            if surrogate_test:
                surrogates = surrogate_generation(sig2, 10, method ="block boostrapping",nsplits=10)

                for l, surrogate in enumerate(surrogates): 
                    sig_gamma, envelope_gamma, phase_gamma = bandpass_filter_and_hilbert_transform(surrogate, fs=fs, f0=f0_gamma, df=df_gamma, norder=norder)

                    # binning the envelope 
                    envelope_gamma_binned = np.ones(len(time_sig1_binned))*np.nan
                    for ii in range(len(time_sig1_binned)):
                        if ii == len(time_sig1_binned)-1:
                            bin_start = time_sig1_binned[-1]
                            bin_end = time_sig1_binned[-1]+0.002
                        else:
                            bin_start = time_sig1_binned[ii]
                            bin_end = time_sig1_binned[ii+1]
                        envelope_gamma_binned[ii] = np.nan if not np.any(np.where((bin_start <= time) & (time < bin_end))) else np.nanmean(envelope_gamma[(bin_start <= time) & (time < bin_end)])

                    amplitude = np.zeros(nbins)       
                    for k in range(nbins):
                        pl = pbins[k]
                        pr = pbins[k+1]                   
                        indices=(phase_theta_sig1_adjusted>=pl) & (phase_theta_sig1_adjusted<pr)
                        amplitude[k] = np.nanmean(envelope_gamma_binned[indices]) 
                    pi_surrogates[l] = amplitude
                amplitude_gamma_surrogates_list.append(pi_surrogates)
    
    data["theta_amplitude"] = np.concatenate(data["theta_amplitude"])
    data["gamma_amplitude"] = np.concatenate(data["gamma_amplitude"])
    data["iseed"] = np.concatenate(data["iseed"])
    data["ibseed"] = np.concatenate(data["ibseed"])
    data = pd.DataFrame(data)

    amplitude_gamma_mean = np.nanmean(amplitude_gamma_list,axis=0)
    amplitude_gamma_sem  = np.nanstd(amplitude_gamma_list,axis=0)/np.sqrt(len(amplitude_gamma_list))
    amplitude_theta_mean = np.nanmean(amplitude_theta_list,axis=0)
    amplitude_theta_sem  = np.nanstd(amplitude_theta_list,axis=0)/np.sqrt(len(amplitude_theta_list))

    amplitude_surrogates = np.concatenate(amplitude_gamma_surrogates_list,axis=0)

    p = np.zeros(np.shape(amplitude_surrogates)[-1])
    for i in range(np.shape(amplitude_surrogates)[-1]):
        # test signifcance 95& both sides
        mu, std = np.nanmean(amplitude_surrogates[:,i]), np.nanstd(amplitude_surrogates[:,i])
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

#######################################################################################################################

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

    outputs = []
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            outputs.append(compute_amplitude_coupling_bycycle(data_new,"soma_ca1_volt_mean",f"{comp1}_volt_mean",f"{comp2}_volt_mean"))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            outputs.append(compute_amplitude_coupling_bycycle(data_new,"soma_ca1_volt_mean",f"{comp1}_volt_mean",f"{comp2}_volt_mean"))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]
    
    output_filename = output_filename.replace(".lzma","_bycycle.lzma")
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

    outputs = []
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            outputs.append(compute_amplitude_coupling_bycycle(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"] 
    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            outputs.append(compute_amplitude_coupling_bycycle(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]

    output_filename = output_filename.replace(".lzma","_bycycle.lzma")
    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder_ca1, "measurements"))


def compute_amplitude_coupling_from_ica(folder_ca3,folder_ca1, input_filename_ca3, input_filename_ca1, output_filename):
    input_filename_ca3  = os.path.join(folder_ca3, input_filename_ca3)
    input_filename_ca1  = os.path.join(folder_ca1, input_filename_ca1)
    data_ca3 = file_management.load_lzma(input_filename_ca3)
    data_ca1 = file_management.load_lzma(input_filename_ca1)

    # empty datafrmae
    electrodes = [-100, 10, 85, 235, 385]
    aux     = ["Bdend","soma","Adend1", "Adend2", "Adend3"]
    aux_ca3 = ["Bdend_ca3","soma_ca3","Adend1_ca3", "Adend2_ca3", "Adend3_ca3"]
    aux_ca1 = ["Bdend_ca1","soma_ca1","Adend1_ca1", "Adend2_ca1", "Adend3_ca1"]

    data_new = {}
    for comp in aux:
        data_new[f"{comp}_ca3"] = []
        data_new[f"{comp}_ca1"] = []
    data_new["iseed"] = []
    data_new["ibseed"] = [] 
    
    # try: 
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

    outputs = []
    for i,comp1 in enumerate(aux_ca3):
        for j, comp2 in enumerate(aux_ca3):
            outputs.append(compute_amplitude_coupling_bycycle(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]

    for i,comp1 in enumerate(aux_ca1):
        for j, comp2 in enumerate(aux_ca1):
            outputs.append(compute_amplitude_coupling_bycycle(data_new,"soma_ca1",comp1,comp2))
            data_dict[f"{comp1}_{comp2}"] = outputs[-1]["average"]

    output_filename = output_filename.replace(".lzma","_bycycle.lzma")
    file_management.save_lzma(data_dict, output_filename, parent_dir = os.path.join(folder_ca1, "measurements"))

