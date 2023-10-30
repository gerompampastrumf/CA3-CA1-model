import numpy as np 
import pandas as pd
import os 
import sys
from neurodsp.filt import filter_signal
from bycycle.cyclepoints import find_extrema, find_zerox
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version')
import file_management 
import signal_analysis as sa

def spiking_pattern_theta_cycle(time_series_list, spikes_rate_list, nbins=40,fs = 1/0.002, f0_theta=8, df_theta=2,norder=4,plot=False, ax=None):
    pbins = np.linspace(0,2*np.pi,nbins+1) # increased the bins to visualized phase difference between the thetas.
    bins = pbins[1:]-np.diff(pbins)[0]/2
    amp_list = []
    for time_series, spikes_rate in zip(time_series_list, spikes_rate_list):
        # plt.plot(spikes_rate)
        sig_theta, envelope_theta, phase_theta = sa.bandpass_filter_and_hilbert_transform(time_series, fs=fs, f0=f0_theta, df=df_theta, norder=norder)

        phase_theta = np.mod(phase_theta+np.pi,2*np.pi)

        amp = np.zeros(nbins)
        amp_theta = np.zeros(nbins)
        for k in range(nbins):         
            pl = pbins[k]
            pr = pbins[k+1]                   
            indices=(phase_theta>=pl) & (phase_theta<pr)  
            amp[k] = np.nanmean(spikes_rate[indices])
            amp_theta[k] = np.nanmean(sig_theta[indices])
        # amp = np.array(list(amp)+[amp[0]])
        amp_list.append(amp)
    
    amp_error = np.nanstd(amp_list,axis=0)/np.sqrt(len(amp_list))
    amp_mean = np.nanmean(amp_list,axis=0)

    # fig, ax = plt.subplots(figsize=(20,5))
    if plot:
        ax.step(bins*180/np.pi,amp_mean,where="mid")
        ax.fill_between(bins*180/np.pi, amp_mean,0,alpha=0.5,step="mid")
        ax.twinx().plot(bins*180/np.pi, amp_theta, color="black")

    return bins, amp_mean, amp_error


def spiking_pattern_theta_cycle_bycycle(time_series_list, spikes_list, time, nbins=40,fs = 1/0.002, f0_theta=8, df_theta=2,norder=4,plot=False, ax=None):

    f_lowpass = 20
    n_seconds_filter = 0.5
    n_seconds_theta = 0.75
    f_theta = ( f0_theta-df_theta, f0_theta+df_theta)

    amp_list = []
    pbins = np.linspace(0,2*np.pi,nbins+1) # increased the bins to visualized phase difference between the thetas.
    bins = pbins[1:]-np.diff(pbins)[0]/2
    for time_series, spikes in zip(time_series_list, spikes_list):
        # plt.plot(spikes_rate)
        # sig_theta, envelope_theta, phase_theta = sa.bandpass_filter_and_hilbert_transform(time_series, fs=fs, f0=f0_theta, df=df_theta, norder=norder)
        # phase_theta = np.mod(phase_theta+np.pi,2*np.pi)
        ##################################################################

        # Lowpass filter (theta reference)
        sig1_low = filter_signal(time_series, fs, 'lowpass', f_lowpass, n_seconds=n_seconds_filter, remove_edges=False)
        # Narrowband filter signal (theta reference)
        sig1_narrow = filter_signal(time_series, fs, 'bandpass', f_theta, n_seconds=n_seconds_theta, remove_edges=False)
        # Find peaks and troughs (this function also does the above)
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

        time_binned,phase_theta_binned = [],[]
        for ii in range(len(troughs)-1):
            time_binned.append(np.linspace(time[troughs[ii]], time[rises[ii]],     nbins//4, endpoint=False))
            time_binned.append(np.linspace(time[rises[ii]],   time[peaks[ii]],     nbins//4, endpoint=False))
            time_binned.append(np.linspace(time[peaks[ii]],   time[decays[ii]],    nbins//4, endpoint=False))
            time_binned.append(np.linspace(time[decays[ii]],  time[troughs[ii+1]], nbins//4, endpoint=False))
            phase_theta_binned.append(np.linspace(0, 2*np.pi, nbins,endpoint=False))

        time_binned = np.concatenate(time_binned)
        phase_theta_binned = np.concatenate(phase_theta_binned)

        # envelope_gamma_binned = np.ones(len(time_binned))*np.nan
        sig_theta_binned = np.ones(len(time_binned))*np.nan
        for ii in range(len(time_binned)):
            if ii == len(time_binned)-1:
                bin_start = time_binned[-1]
                bin_end = time_binned[-1]+0.002
            else:
                bin_start = time_binned[ii]
                bin_end = time_binned[ii+1]
            sig_theta_binned[ii] = np.nan if not np.any(np.where((bin_start <= time) & (time < bin_end))) else np.nanmean(sig1_low[(bin_start <= time) & (time < bin_end)])
        
        spikes_binned = np.histogram(spikes, bins=time_binned.tolist()+[time_binned[-1]+0.002], density=True)[0]
        # this can generate a warning due to it can give dbin = 0 in some cases.! 
        
        ##########################################################################
        amp = np.zeros(nbins)
        amp_theta = np.zeros(nbins)
        for k in range(nbins):         
            pl = pbins[k]
            pr = pbins[k+1]                   
            indices=(phase_theta_binned>=pl) & (phase_theta_binned<pr)  
            amp[k] = np.nanmean(spikes_binned[indices])
            amp_theta[k] = np.nanmean(sig_theta_binned[indices])
        amp_list.append(amp)
    
    amp_error = np.nanstd(amp_list,axis=0)/np.sqrt(len(amp_list))
    amp_mean = np.nanmean(amp_list,axis=0)

    # bins = pbins[1:]-np.diff(pbins)[0]/2
    # fig, ax = plt.subplots(figsize=(20,5))
    if plot:
        ax.step(bins*180/np.pi,amp_mean,where="mid")
        ax.fill_between(bins*180/np.pi, amp_mean,0,alpha=0.5,step="mid")
        ax.twinx().plot(bins*180/np.pi, amp_theta, color="black")

    return bins, amp_mean, amp_error


def compute_spiking_patterns_with_dynamic_reference(folder_ca3,folder_ca1):
    ''' Three different ways: 

    1) CA3 dynamic reference: volt
    2) CA3 dynmica reference: lfp
    3) CA3 dynamic reference: ica
    4) CA3 dynamic reference: lfp alt
    5) CA1 dynamic reference: volt
    6) CA1 dynamic reference: lfp
    7) CA1 dynamic reference: ica
    8) CA1 dynamic reference: lfp alt

    With both hilbert and alternative method of signal binning 
    '''

    folder = dict.fromkeys(["ca3","ca1"])
    folder["ca3"] = folder_ca3
    folder["ca1"] = folder_ca1

    spikes_files = dict.fromkeys(["ca3","ca1"])
    ref_files    = dict.fromkeys(["ca3","ca1"])

    spikes_files["ca3"] = ["external_inputs_pyr_ca3.lzma","spikes_bas_ca3.lzma","spikes_olm_ca3.lzma"]
    spikes_files["ca1"] = ["spikes_pyr_ca1.lzma","spikes_bas_ca1.lzma","spikes_olm_ca1.lzma","spikes_cck_ca1.lzma"]

    ref_files["ca3"] = ["volt_pyr_ca3.lzma", "lfp_ca3.lzma", "ica_ca3.lzma"]
    ref_files["ca1"] = ["volt_pyr_ca1.lzma", "lfp_ca1.lzma", "ica_ca1.lzma"]

    spiking_patterns = dict.fromkeys(["ca3","ca1"])
    spiking_patterns["ca3"] = dict.fromkeys(["volt_soma","lfp_soma","ica_soma","lfp_alt"])
    spiking_patterns["ca1"] = dict.fromkeys(["volt_soma","lfp_soma","ica_soma","lfp_alt"])

    spiking_patterns_bycycle = dict.fromkeys(["ca3","ca1"])
    spiking_patterns_bycycle["ca3"] = dict.fromkeys(["volt_soma","lfp_soma","ica_soma","lfp_alt"])
    spiking_patterns_bycycle["ca1"] = dict.fromkeys(["volt_soma","lfp_soma","ica_soma","lfp_alt"])
    
    for r1 in ["ca3","ca1"]:

        spiking_patterns[r1]["volt_soma"] = {}
        spiking_patterns[r1]["lfp_soma"]  = {}
        spiking_patterns[r1]["ica_soma"]  = {}
        spiking_patterns[r1]["lfp_alt"]   = {}

        spiking_patterns_bycycle[r1]["volt_soma"] = {}
        spiking_patterns_bycycle[r1]["lfp_soma"]  = {}
        spiking_patterns_bycycle[r1]["ica_soma"]  = {}
        spiking_patterns_bycycle[r1]["lfp_alt"]   = {}

        data1 = file_management.load_lzma(os.path.join(folder[r1], ref_files[r1][0]))
        data2 = file_management.load_lzma(os.path.join(folder[r1], ref_files[r1][1]))
        data3 = file_management.load_lzma(os.path.join(folder[r1], ref_files[r1][2]))

        lfp_soma_aux = data2[(data2["iseed"] == 0) & (data2["ibseed"] == 0) & (data2["electrode"]==10)]["lfp"].values
        n = len(lfp_soma_aux)
        time = np.linspace(0,n*0.002,n)
        wt  = (time>=2) & (time<=30)
        time = time[wt]

        volt_soma_list   = []
        lfp_soma_list    = []
        ica_soma_list    = []
        lfp_alt_list     = []

        for i in range(10):
            for j in range(10):

                w1 = (data1["iseed"] == i) & (data1["ibseed"] == j)
                w2 = (data2["iseed"] == i) & (data2["ibseed"] == j) & (data2["electrode"]==10)
                w3 = (data3["iseed"] == i) & (data3["ibseed"] == j) & (data3["electrode"]==10) & (data3["component"] == "soma")

                volt_soma_list.append( data1[w1]["soma_volt_mean"].values[:-1][wt])
                lfp_soma_list.append( data2[w2]["lfp"].values[wt])
                ica_soma_list.append( data3[w3]["ica"].values[wt])
                lfp_alt_list.append( data1[w1]["Adend3_volt_mean"].values[:-1][wt]-data1[w1]["Bdend_volt_mean"].values[:-1][wt])

        for r2 in ["ca3","ca1"]:
            for file in spikes_files[r2]:
                label = file.split("_")[1]
                print(folder[r2], file)
                data0 = file_management.load_lzma(os.path.join(folder[r2], file))
                
                spikes_rate_list = []
                spikes_list = [] # in the bycycle method the binning is done once the ref is binned
                for i in range(10):
                    for j in range(10):
                            
                        w0 = (data0["iseed"] == i) & (data0["ibseed"] == j)
                        spikes = data0[w0]["tvec"].values/1e3
                        spikes = spikes[(spikes>=2) & (spikes<=30)]
                        spikes_rate_list.append(np.histogram(spikes, bins=time.tolist()+[time[-1]+0.002], density=True)[0])
                        spikes_list.append(spikes)
                
                # regular method 
                bins, amp_mean, amp_error = spiking_pattern_theta_cycle(volt_soma_list, spikes_rate_list)
                spiking_patterns[r1]["volt_soma"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns[r1]["volt_soma"][label+f"_{r2}_error"] = amp_error

                bins, amp_mean, amp_error = spiking_pattern_theta_cycle(lfp_soma_list, spikes_rate_list)
                spiking_patterns[r1]["lfp_soma"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns[r1]["lfp_soma"][label+f"_{r2}_error"] = amp_error

                bins, amp_mean, amp_error = spiking_pattern_theta_cycle(ica_soma_list, spikes_rate_list)
                spiking_patterns[r1]["ica_soma"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns[r1]["ica_soma"][label+f"_{r2}_error"] = amp_error

                bins, amp_mean, amp_error = spiking_pattern_theta_cycle(lfp_alt_list, spikes_rate_list)
                spiking_patterns[r1]["lfp_alt"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns[r1]["lfp_alt"][label+f"_{r2}_error"] = amp_error 

                # bycycle method
                bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(volt_soma_list, spikes_list, time)
                spiking_patterns_bycycle[r1]["volt_soma"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns_bycycle[r1]["volt_soma"][label+f"_{r2}_error"] = amp_error

                bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(lfp_soma_list, spikes_list, time)
                spiking_patterns_bycycle[r1]["lfp_soma"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns_bycycle[r1]["lfp_soma"][label+f"_{r2}_error"] = amp_error

                bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(ica_soma_list, spikes_list, time)
                spiking_patterns_bycycle[r1]["ica_soma"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns_bycycle[r1]["ica_soma"][label+f"_{r2}_error"] = amp_error

                bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(lfp_alt_list, spikes_list, time)
                spiking_patterns_bycycle[r1]["lfp_alt"][label+f"_{r2}_mean"]  = amp_mean
                spiking_patterns_bycycle[r1]["lfp_alt"][label+f"_{r2}_error"] = amp_error
                
        
        spiking_patterns[r1]["volt_soma"]["bins"] = bins
        spiking_patterns[r1]["lfp_soma"]["bins"]  = bins
        spiking_patterns[r1]["ica_soma"]["bins"]  = bins
        spiking_patterns[r1]["lfp_alt"]["bins"]   = bins

        spiking_patterns_bycycle[r1]["volt_soma"]["bins"] = bins
        spiking_patterns_bycycle[r1]["lfp_soma"]["bins"]  = bins
        spiking_patterns_bycycle[r1]["ica_soma"]["bins"]  = bins
        spiking_patterns_bycycle[r1]["lfp_alt"]["bins"]   = bins


        file_management.save_lzma(spiking_patterns[r1],f"spiking_patterns_ref_{r1}.lzma",parent_dir=os.path.join(folder_ca1,"measurements"))
        file_management.save_lzma(spiking_patterns_bycycle[r1],f"spiking_patterns_ref_{r1}_bycycle.lzma",parent_dir=os.path.join(folder_ca1,"measurements"))

def compute_spiking_patterns_with_dynamic_reference_external_inputs(folder_ca3,folder_ca1):
    ''' Three different ways: 
    1) CA3 dynamic reference: volt
    2) CA3 dynmica reference: lfp
    3) CA3 dynamic reference: ica
    4) CA3 dynamic reference: lfp alt
    5) CA1 dynamic reference: volt
    6) CA1 dynamic reference: lfp
    7) CA1 dynamic reference: ica
    8) CA1 dynamic reference: lfp alt
    '''

    folder = {"ca3": folder_ca3, "ca1":folder_ca1}

    spikes_files = ["external_inputs_sep_180.lzma","external_inputs_sep_360.lzma","external_inputs_ec2_180.lzma",
                    "external_inputs_ec3_360.lzma","external_inputs_dg_regular.lzma","external_inputs_dg_burst.lzma"]
    
    aux = ["sep_180","sep_360","ec2_180","ec3_360","dg_regular","dg_burst"]
    spiking_patterns = dict.fromkeys(["ca3","ca1"])
    spiking_patterns_bycycle = dict.fromkeys(["ca3","ca1"])
    for r1 in ["ca3","ca1"]:
        spiking_patterns[r1] = {}
        spiking_patterns_bycycle[r1] = {} 

    ref_files    = dict.fromkeys(["ca3","ca1"])
    ref_files["ca3"] = ["volt_pyr_ca3.lzma", "lfp_ca3.lzma", "ica_ca3.lzma"]
    ref_files["ca1"] = ["volt_pyr_ca1.lzma", "lfp_ca1.lzma", "ica_ca1.lzma"]

    for r1 in ["ca3","ca1"]:

        spiking_patterns[r1]["volt_soma"] = {}
        spiking_patterns[r1]["lfp_soma"]  = {}
        spiking_patterns[r1]["ica_soma"]  = {}
        spiking_patterns[r1]["lfp_alt"]   = {}

        spiking_patterns_bycycle[r1]["volt_soma"] = {}
        spiking_patterns_bycycle[r1]["lfp_soma"]  = {}
        spiking_patterns_bycycle[r1]["ica_soma"]  = {}
        spiking_patterns_bycycle[r1]["lfp_alt"]   = {}

        data1 = file_management.load_lzma(os.path.join(folder[r1], ref_files[r1][0]))
        data2 = file_management.load_lzma(os.path.join(folder[r1], ref_files[r1][1]))
        data3 = file_management.load_lzma(os.path.join(folder[r1], ref_files[r1][2]))

        lfp_soma_aux = data2[(data2["iseed"] == 0) & (data2["ibseed"] == 0) & (data2["electrode"]==10)]["lfp"].values
        n = len(lfp_soma_aux)
        time = np.linspace(0,n*0.002,n)
        wt  = (time>=2) & (time<=30)
        time = time[wt]

        volt_soma_list   = []
        lfp_soma_list    = []
        ica_soma_list    = []
        lfp_alt_list     = []

        for i in range(10):
            for j in range(10):

                w1 = (data1["iseed"] == i) & (data1["ibseed"] == j)
                w2 = (data2["iseed"] == i) & (data2["ibseed"] == j) & (data2["electrode"]==10)
                w3 = (data3["iseed"] == i) & (data3["ibseed"] == j) & (data3["electrode"]==10) & (data3["component"] == "soma")

                volt_soma_list.append( data1[w1]["soma_volt_mean"].values[:-1][wt])
                lfp_soma_list.append( data2[w2]["lfp"].values[wt])
                ica_soma_list.append( data3[w3]["ica"].values[wt])
                lfp_alt_list.append( data1[w1]["Adend3_volt_mean"].values[:-1][wt]-data1[w1]["Bdend_volt_mean"].values[:-1][wt])

        for file,label in zip(spikes_files,aux):
            print(file,label)

            data0 = file_management.load_lzma(os.path.join(folder["ca3"], file))
            spikes_rate_list = []
            spikes_list = [] # in the bycycle method the binning is done once the ref is binned

            for i in range(10):
                for j in range(10):
                    w0 = (data0["iseed"] == i) #& (data0["ibseed"] == j)
                    spikes = data0[w0]["tvec"].values/1e3
                    spikes = spikes[(spikes>=2) & (spikes<=30)]
                    spikes_rate_list.append(np.histogram(spikes, bins=time.tolist()+[time[-1]+0.002], density=True)[0])
                    spikes_list.append(spikes)
            
            # regular method 
            bins, amp_mean, amp_error = spiking_pattern_theta_cycle(volt_soma_list, spikes_rate_list)
            spiking_patterns[r1]["volt_soma"][label+f"_mean"]  = amp_mean
            spiking_patterns[r1]["volt_soma"][label+f"_error"] = amp_error

            bins, amp_mean, amp_error = spiking_pattern_theta_cycle(lfp_soma_list, spikes_rate_list)
            spiking_patterns[r1]["lfp_soma"][label+f"_mean"]  = amp_mean
            spiking_patterns[r1]["lfp_soma"][label+f"_error"] = amp_error

            bins, amp_mean, amp_error = spiking_pattern_theta_cycle(ica_soma_list, spikes_rate_list)
            spiking_patterns[r1]["ica_soma"][label+f"_mean"]  = amp_mean
            spiking_patterns[r1]["ica_soma"][label+f"_error"] = amp_error

            bins, amp_mean, amp_error = spiking_pattern_theta_cycle(lfp_alt_list, spikes_rate_list)
            spiking_patterns[r1]["lfp_alt"][label+f"_mean"]  = amp_mean
            spiking_patterns[r1]["lfp_alt"][label+f"_error"] = amp_error 

            # bycycle method
            bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(volt_soma_list, spikes_list, time)
            spiking_patterns_bycycle[r1]["volt_soma"][label+f"_mean"]  = amp_mean
            spiking_patterns_bycycle[r1]["volt_soma"][label+f"_error"] = amp_error

            bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(lfp_soma_list, spikes_list, time)
            spiking_patterns_bycycle[r1]["lfp_soma"][label+f"_mean"]  = amp_mean
            spiking_patterns_bycycle[r1]["lfp_soma"][label+f"_error"] = amp_error

            bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(ica_soma_list, spikes_list, time)
            spiking_patterns_bycycle[r1]["ica_soma"][label+f"_mean"]  = amp_mean
            spiking_patterns_bycycle[r1]["ica_soma"][label+f"_error"] = amp_error

            bins, amp_mean, amp_error = spiking_pattern_theta_cycle_bycycle(lfp_alt_list, spikes_list, time)
            spiking_patterns_bycycle[r1]["lfp_alt"][label+f"_mean"]  = amp_mean
            spiking_patterns_bycycle[r1]["lfp_alt"][label+f"_error"] = amp_error
            
        spiking_patterns[r1]["volt_soma"]["bins"] = bins
        spiking_patterns[r1]["lfp_soma"]["bins"]  = bins
        spiking_patterns[r1]["ica_soma"]["bins"]  = bins
        spiking_patterns[r1]["lfp_alt"]["bins"]   = bins

        spiking_patterns_bycycle[r1]["volt_soma"]["bins"] = bins
        spiking_patterns_bycycle[r1]["lfp_soma"]["bins"]  = bins
        spiking_patterns_bycycle[r1]["ica_soma"]["bins"]  = bins
        spiking_patterns_bycycle[r1]["lfp_alt"]["bins"]   = bins


        file_management.save_lzma(spiking_patterns[r1],f"spiking_patterns_external_inputs_ref_{r1}.lzma",parent_dir=os.path.join(folder["ca1"],"measurements"))
        file_management.save_lzma(spiking_patterns_bycycle[r1],f"spiking_patterns_external_inputs_ref_{r1}_bycycle.lzma",parent_dir=os.path.join(folder["ca1"],"measurements"))

