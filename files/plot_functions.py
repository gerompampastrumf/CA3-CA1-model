import numpy as np
from scipy import signal as sp
import sys 
import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import glob
import shutil
import scipy
from matplotlib.colorbar import Colorbar
from spectrum import *
from sklearn.neighbors import KernelDensity

#sys.path.append("/home/jaime/Desktop/hippocampus/files/")
import file_management
import matplotlib as mpl
mpl.rcParams["figure.max_open_warning"]=0
mpl.rcParams['pcolor.shading']="auto"

# #####################################################################################ç
# Analysis functions 
#######################################################################################
def spline_fit(x,y,s=5, ndim=100):
    ''' spline fit with degree k=3 (defaualt). s is the number of knots '''
    from scipy.interpolate import UnivariateSpline
    us = UnivariateSpline(x, y, s=s)
    xs = np.linspace(0, x[-1], ndim)
    return xs, us(xs)

def smoothing(x,ymean,ysem, s):
    ''' spine fit of the the average of a time series with its error'''
    xs, ymean_fit = spline_fit(x, ymean, s)
    xs, ysup_fit  = spline_fit(x, ymean+ysem, s)
    xs, yinf_fit  = spline_fit(x, ymean-ysem, s)
    return xs, ymean_fit, ysup_fit, yinf_fit 

def get_spiking_phase_distribution(xlist, ylist, phase_reference, simulation_time=4000, nbins=20, dt=0.1): 
    time = np.arange(0,simulation_time-dt,dt)
    nsims = len(xlist) 

    starting_cycles = np.where( np.abs( np.diff( phase_reference )) > 1.9*np.pi)[0]
    ncycles =  len(starting_cycles)-1
    nspikes = np.zeros((nsims*ncycles, nbins))
    nspikes_std = np.zeros((nsims*ncycles, nbins))
    
    for k, (x, y) in enumerate(zip(xlist, ylist)): 

        for ii in range(ncycles):
            n1,n2 = starting_cycles[ii]+1, starting_cycles[ii+1]
            t1,t2 = time[n1],time[n2]
            phi1, phi2 = phase_reference[n1], phase_reference[n2]

            time_aux = np.linspace(t1,t2,nbins+1)
            phi_aux  = np.linspace(phi1,phi2, nbins+1)
            dphi = np.diff(phi_aux[:2])
            phi_seq = phi_aux[:-1]+dphi/2.0

            for jj in range(len(time_aux[:-1])):
                tj1,tj2 = time_aux[jj],time_aux[jj+1]
                w = np.logical_and(x >= tj1, x < tj2)
                nspikes[(k+1)*ii,jj] = len( x[w] )

    spikes_mean = np.mean(nspikes, axis=0)
    spikes_std  = np.std(nspikes, axis=0)
    spikes_sem  = np.std(nspikes, axis=0)/np.sqrt(ncycles*nsims)

    spike_phase_peak = phi_seq[np.argmax(spikes_mean)]*180/np.pi

    phiseq = np.array( phi_seq.tolist() + (phi_seq+2*np.pi).tolist() )*180/np.pi
    spike_phase_mean = np.array( spikes_mean.tolist()*2)
    spike_phase_std  = np.array( spikes_std.tolist()*2 )
    spike_phase_sem  = np.array( spikes_sem.tolist()*2 )
    return phiseq, spike_phase_mean, spike_phase_std, spike_phase_sem, spike_phase_peak

def get_nfft(m):
    nfft = 2
    k = 0
    while m > nfft:
        k+=1
        nfft = 2**k
    return nfft

def bandpass_filter(xs, norder, f_range, fs=1e4): 
    sos = sp.butter(N=norder, Wn=f_range, btype="bandpass", fs=fs, output="sos")
    return sp.sosfilt(sos,xs)

def get_mean_activity_per_cycle(xlist, ylist, theta_cycle=125):
    '''
    xlist: different realizations list, spikes
    ylist: different realizations list, neuron number
    '''
    n = 8
    initial_time = 50 + n*theta_cycle
    
    xaux = np.concatenate(xlist)
    final_time = xaux.max()-np.mod(xaux.max(),theta_cycle)
    
    nspikes = []
    active_neurons = []
    active_same_neurons = []
    for x, y in zip(xlist, ylist):
        t1, t2 = initial_time, initial_time+theta_cycle
        while t2<=final_time:
            w = np.logical_and(x>=t1, x<t2)
            nspikes.append( len(x[w]) ) 
            active_neurons.append( len(np.unique(y[w])) )
            
            yaux = Counter(y[w]) 
            yaux = np.array(list(yaux.values()))
            active_same_neurons.append( len(yaux[np.where(yaux>1)]) )
            
            t1+=theta_cycle
            t2+=theta_cycle
    
    ns = len(nspikes)
    nspikes_std = np.round( np.std(nspikes)/np.sqrt(ns),2)
    nspikes_mean = np.round( np.mean(nspikes), 2)
    active_neurons_std = np.round( np.std(active_neurons)/np.sqrt(ns), 2)
    active_neurons_mean = np.round( np.mean(active_neurons), 2)
    active_same_neurons_std = np.round( np.std(active_same_neurons)/np.sqrt(ns), 2)
    active_same_neurons_mean = np.round( np.mean(active_same_neurons), 2) 
    
    return nspikes_mean, nspikes_std, active_neurons_mean, active_neurons_std, active_same_neurons_mean, active_same_neurons_std

def get_voltage_diff_per_cycle(t,x):
    theta_cycle = 125
    initial_time = 50 + 5*theta_cycle
    t1, t2 = initial_time, initial_time+theta_cycle 
    final_time = t[-1]-np.mod(t[-1], theta_cycle)
    
    vdiff = []
    vmin, vmax = [],[]
    while t2<=final_time:
        w = np.logical_and(t>=t1, t<t2)
        
        vmin.append( np.min( x[w] ) ) 
        vmax.append( np.max( x[w] ) )

        vdiff.append( vmax[-1]-vmin[-1] )    

        t1+=theta_cycle
        t2+=theta_cycle

    vdiff_mean = np.mean(vdiff)
    vdiff_sem  = np.std(vdiff)/np.sqrt(len(vdiff))
    vmin_mean  = np.mean(vmin)                                     
    vmin_sem   = np.std(vmin)/np.sqrt(len(vdiff))
    vmax_mean  = np.mean(vmax)                                     
    vmax_sem   = np.std(vmax)/np.sqrt(len(vdiff))
    return vmin_mean, vmin_sem, vmax_mean, vmax_sem, vdiff_mean, vdiff_sem

# ###############################################################################
# DataFrames generation
# ###############################################################################
def get_lfp_and_ica_DataFrames(folder, save_dir, label, net="ca1"): 
    '''In principle only working for agosto folder'''
    
    files = os.listdir(folder)
    file_names = []
    for file in files: 
        file_names.append(file.split(".")[0].split('_')[0])
    
    lfp_df = pd.DataFrame()
    ica_df = pd.DataFrame()
    
    if files:
        file_names = np.unique(file_names)
        aux = files[0].split(".")[0].split('_')
        all_labels = ["k","i","j"]
        nl = len(aux)-3
        if nl == 3: 
            files_label = label+"_double_scan"
        else: 
            files_label = label
        
        files_ = {}
        if "lfp" in file_names: 
            files_["lfp"] = []
            for file in files: 
                aux = file.split(".")[0].split('_')[0]
                aux2 = len(file.split(".")[0].split('_'))
                if (aux == "lfp") & (aux2 == nl+3): 
                    files_["lfp"].append(file)

            for file in files_["lfp"]:
                aux = file.split(".")[0].split('_')[-nl:]
                data = file_management.load_lzma(os.path.join(folder, file))[net]
                aux_df = pd.DataFrame()
                for i,xl in enumerate(data):
                    aux_df[str(i)] = xl[::10].astype("float32")
                for l, label in enumerate(all_labels[:nl]):
                    aux_df[label] = int(aux[l])
                lfp_df = pd.concat([lfp_df, aux_df])

            for l, label in enumerate(all_labels[:nl]):
                lfp_df[label] = lfp_df[label].astype("int8")
            lfp_df.sort_values(by=all_labels[:nl])
                
            print("LFP DataFrame created")
            print("------------------------")
            file_management.save_lzma( lfp_df, "lfp_df_"+files_label+".lzma", parent_dir=save_dir)
            print("LFP DataFrame saved")
            print("------------------------")
            for file in files_["lfp"]: 
                os.remove( os.path.join(folder, file ))
            print("All ica files data removed")
            print("------------------------")
            
        if "ica" in file_names: 
            files_["ica"] = []
            for file in files: 
                aux = file.split(".")[0].split('_')[0]
                aux2 = len(file.split(".")[0].split('_'))
                if (aux == "ica") & (aux2 == nl+3):
                    files_["ica"].append(file)

            for file in files_["ica"]:
                aux = file.split(".")[0].split('_')[-nl:]
                data = file_management.load_lzma(os.path.join(folder, file))[net]
                aux_df = pd.DataFrame()
                for i,xl in enumerate(data):
                    aux_df[str(i)] = xl[::10].astype("float32")
                for l, label in enumerate(all_labels[:nl]):
                    aux_df[label] = int(aux[l])
                ica_df = pd.concat([ica_df, aux_df])

            for l, label in enumerate(all_labels[:nl]):
                ica_df[label] = ica_df[label].astype("int8")
            ica_df.sort_values(by=all_labels[:nl])

            print("Ica DataFrame created")
            print("------------------------")
            file_management.save_lzma( ica_df, "ica_df_"+files_label+".lzma", parent_dir=save_dir)
            print("Ica DataFrame saved")
            print("------------------------")
            for file in files_["ica"]: 
                os.remove( os.path.join(folder, file ))
            print("All ica files data removed")
            print("------------------------")
            
    else:
        print("No files") 
        shutil.rmtree(folder)
        print("folder deteled")
            
    return lfp_df, ica_df            
            
def get_DataFrames_general(folder, save_dir, label, net="ca1"):
    ''' DataFrames are created based on the filename and then those files are removed, 
        in case there are no files, folder is removed'''

    files = os.listdir(folder)
    file_names = []
    for file in files: 
        file_names.append(file.split(".")[0].split('_')[0])
    
    voltage_df = pd.DataFrame()
    spikes_df = pd.DataFrame()
    
    if files:
        file_names = np.unique(file_names)
        aux = files[0].split(".")[0].split('_')
        all_labels = ["k","i","j"]
        nl = len(aux)-3
        if nl == 3: 
            files_label = label+"_double_scan"
        else: 
            files_label = label
        print("file_names:", file_names) 
        print("------------------------")
        files_ = {}
        if "spikes" in file_names: 
            spikes_df = pd.DataFrame()
            files_["spikes"] = []
            for file in files: 
                aux = file.split(".")[0].split('_')[0]
                aux2 = len(file.split(".")[0].split('_'))
                if (aux == "spikes") & (aux2 == nl+3): 
                    files_["spikes"].append(file)

            for file in files_["spikes"]:
                aux = file.split(".")[0].split('_')[-nl:]
                data = file_management.load_lzma(os.path.join(folder, file))[net]
                cells = data.keys() 
                for cell in cells: 
                    aux_df = pd.DataFrame()
                    aux_df["time"] = data[cell]["tvec"]
                    aux_df["id"] = data[cell]["idvec"]
                    aux_df["cell"] = cell
                    for l, label in enumerate(all_labels[:nl]):
                        aux_df[label] = int(aux[l])
                    spikes_df = pd.concat([spikes_df, aux_df])

            for l, label in enumerate(all_labels[:nl]):
                spikes_df[label] = spikes_df[label].astype("int16")
            spikes_df["id"] = spikes_df["id"].astype("int16")
            spikes_df["time"] = spikes_df["time"].astype("float32")
            spikes_df = spikes_df.sort_values(by=all_labels[:nl])
            
            print("Spikes DataFrame created")
            print("------------------------")
            file_management.save_lzma( spikes_df, "spikes_df_"+files_label+".lzma", parent_dir=save_dir)
            print("Spikes DataFrame saved")
            print("------------------------")
            #for file in files_["spikes"]: 
                #os.remove( os.path.join(folder, file ))
            print("All spiking files data removed")

        if "volt" in file_names: 
            voltage_df = pd.DataFrame()
            files_["volt"] = []
            for file in files: 
                aux = file.split(".")[0].split('_')[0]
                aux2 = len(file.split(".")[0].split('_'))
                if (aux == "volt") & (aux2 == nl+3):
                    files_["volt"].append(file)

            for file in files_["volt"]:
                aux = file.split(".")[0].split('_')[-nl:]
                data = file_management.load_lzma(os.path.join(folder, file))[net]
                cells = data.keys() 
                for cell in cells: 
                    aux_df = pd.DataFrame()
                    aux_df["volt"] = data[cell][0][::10]
                    aux_df["cell"] = cell
                    for l, label in enumerate(all_labels[:nl]):
                        aux_df[label] = int(aux[l])
                    voltage_df = pd.concat([voltage_df, aux_df])

            #for l, label in enumerate(all_labels[:nl]):
            #    voltage_df[label] = voltage_df[label].astype("int16")
            voltage_df["volt"] = voltage_df["volt"].astype("float32")
            voltage_df = voltage_df.sort_values(by=all_labels[:nl])

            print("Voltage DataFrame created")
            print("------------------------")
            file_management.save_lzma( voltage_df, "voltage_df_"+files_label+".lzma", parent_dir=save_dir)
            print("Voltage DataFrame saved")
            print("------------------------")
            #for file in files_["volt"]: 
                #os.remove( os.path.join(folder, file ))
            print("All voltage files data removed")
            print("------------------------")
    else:
        print("No files") 
        shutil.rmtree(folder)
        print("folder deteled")
            
    return voltage_df, spikes_df

def get_DataFrames(folder, save_dir, label):
    ''' Specific when only one parameter is scanned'''
    ntrials = 3
    n = 20
    extra_label_title = "only_ca1_"
    cells = ["pyr","bas","olm", "cck"]

    spikes_df = pd.DataFrame()
    voltage_df = pd.DataFrame() 

    for k in range(ntrials):
        for i in range(n): 
            argvs = f"{k}_{i}.lzma"
            title = "spikes_"+extra_label_title+argvs
            file = os.path.join(folder, title)
            if os.path.exists(file):
                data = file_management.load_lzma( file )["ca1"]
                for cell in cells: 
                    aux = pd.DataFrame() 
                    aux["time"] = data[cell]["tvec"]
                    aux["id"] = data[cell]["idvec"]
                    aux["cell"] = cell
                    aux["k"] = k
                    aux["i"] = i 
                    spikes_df = pd.concat([spikes_df, aux])
            else: 
                for cell in cells:
                    aux = pd.DataFrame() 
                    aux["time"] = np.nan
                    aux["id"] = np.nan
                    aux["cell"] = cell
                    aux["k"] = k
                    aux["i"] = i 
                    spikes_df = pd.concat([spikes_df, aux])
                
    print("Spikes DataFrame created")

    for k in range(ntrials):
        for i in range(n): 
            argvs = f"{k}_{i}.lzma"
            title = "volt_"+extra_label_title+argvs
            file = os.path.join(folder, title)
            if os.path.exists(file):
                data = file_management.load_lzma( file )["ca1"]
                for cell in cells: 
                    aux = pd.DataFrame() 
                    aux["volt"] = data[cell][0][::10]
                    #aux["volt_std"] = data[cell][1][::10]
                    aux["cell"] = cell 
                    aux["k"] = k 
                    aux["i"] = i
                    voltage_df = pd.concat([voltage_df,aux])
            else: 
                for cell in cells: 
                    aux = pd.DataFrame() 
                    aux["volt"] = np.nan
                    #aux["volt_std"] = data[cell][1][::10]
                    aux["cell"] = cell 
                    aux["k"] = k 
                    aux["i"] = i
                    voltage_df = pd.concat([voltage_df,aux])
                
    print("Voltage DataFrame created")

    voltage_df["k"] = voltage_df["k"].astype("int8")
    voltage_df["i"] = voltage_df["i"].astype("int8")
    voltage_df["volt"] = voltage_df["volt"].astype("float32")

    spikes_df["k"] = spikes_df["k"].astype("int8")
    spikes_df["i"] = spikes_df["i"].astype("int8")
    spikes_df["id"] = spikes_df["id"].astype("int16")
    spikes_df["time"] = spikes_df["time"].astype("float32")

    file_management.save_lzma(voltage_df, "voltage_df_"+label+".lzma",parent_dir = save_dir )
    file_management.save_lzma(spikes_df,  "spikes_df_"+label+".lzma", parent_dir = save_dir )
    print("DataFrames saved")
    
    return voltage_df, spikes_df

def get_data_per_cycle(runtime, voltage_df, spikes_df, save_dir, label):

    n = np.max(voltage_df["i"])+1
    ntrials = np.max(voltage_df["k"])+1
    cells = voltage_df["cell"].unique()
    w = (voltage_df["cell"]=="pyr") & (voltage_df["k"]==0) & (voltage_df["i"]==0)
    nn = len( voltage_df["volt"][w] ) 
    time = np.linspace(0,runtime, nn)

    keys = ["nspikes_mean", "nspikes_sem",
           "active_mean", "active_sem",
           "active_same_mean", "active_same_sem",
           "vmax_mean", "vmax_sem",
           "vmin_mean", "vmin_sem",
           "vdiff_mean", "vdiff_sem",
           "cell","k", "i"] 

    aux = dict.fromkeys(keys)
    for key in keys: 
        aux[key] = []

    for k in range(ntrials): 
        for i in range(n): 
            for cell in cells: 
                w = (spikes_df["cell"]==cell) & (spikes_df["k"]==k) & (spikes_df["i"]==i) & (spikes_df["time"]>500.0)
                xlist = np.array( spikes_df["time"][w] )
                ylist = np.array( spikes_df["id"][w] )

                w = (voltage_df["cell"]==cell) & (voltage_df["k"]==k) & (voltage_df["i"]==i)
                volt  = np.array( voltage_df["volt"][w] )

                if len( xlist) > 0:
                    output = get_mean_activity_per_cycle([xlist], [ylist], theta_cycle=125)
                    aux["nspikes_mean"].append(output[0])
                    aux["nspikes_sem"].append(output[1])
                    aux["active_mean"].append(output[2])
                    aux["active_sem"].append(output[3]) 
                    aux["active_same_mean"].append(output[4])
                    aux["active_same_sem"].append(output[5])
                else: 
                    aux["nspikes_mean"].append(np.nan)
                    aux["nspikes_sem"].append(np.nan)
                    aux["active_mean"].append(np.nan )
                    aux["active_sem"].append(np.nan)
                    aux["active_same_mean"].append(np.nan )
                    aux["active_same_sem"].append(np.nan)

                if len(volt) > 0:
                    output = get_voltage_diff_per_cycle( time, volt )
                    aux["vmin_mean"].append(output[0])
                    aux["vmin_sem"].append(output[1])
                    aux["vmax_mean"].append(output[2])
                    aux["vmax_sem"].append(output[3])
                    aux["vdiff_mean"].append(output[4])
                    aux["vdiff_sem"].append(output[5])
                else:
                    aux["vmin_mean"].append(np.nan)
                    aux["vmin_sem"].append(np.nan)
                    aux["vmax_mean"].append(np.nan)
                    aux["vmax_sem"].append(np.nan)
                    aux["vdiff_mean"].append(np.nan)
                    aux["vdiff_sem"].append(np.nan)

                aux["cell"].append( cell )  
                aux["k"].append(k) 
                aux["i"].append(i) 

    data_df = pd.DataFrame( aux )
    file_management.save_lzma( data_df, "data_df_"+label+".lzma",parent_dir = save_dir )
    print("DataFrame created and saved")

def get_data_per_cycle_double_scan(runtime, voltage_df, spikes_df, save_dir, label):
    ''' for double scan simulations ''' 
    n1 = np.max(voltage_df["i"])+1
    n2 = np.max(voltage_df["j"])+1 
    ntrials = np.max(voltage_df["k"])+1
    cells = voltage_df["cell"].unique()
    w = (voltage_df["cell"]=="pyr") & (voltage_df["k"]==0) & (voltage_df["i"]==0) & (voltage_df["j"] == 0)
    nn = len( voltage_df["volt"][w] ) 
    time = np.linspace(0,runtime, nn)

    keys = ["nspikes_mean", "nspikes_sem",
           "active_mean", "active_sem",
           "active_same_mean", "active_same_sem",
           "vmax_mean", "vmax_sem",
           "vmin_mean", "vmin_sem",
           "vdiff_mean", "vdiff_sem",
           "cell","k","i", "j"] 

    aux = dict.fromkeys(keys)
    for key in keys: 
        aux[key] = []

    for k in range(ntrials): 
        for i in range(n1): 
            for j in range(n2):
                for cell in cells: 
                    w = (spikes_df["cell"]==cell) & (spikes_df["k"]==k) & (spikes_df["i"]==i) & (spikes_df["j"]==j) & (spikes_df["time"]>550.0)
                    xlist = np.array( spikes_df["time"][w] )
                    ylist = np.array( spikes_df["id"][w] )

                    w = (voltage_df["cell"]==cell) & (voltage_df["k"]==k) & (voltage_df["j"]==j) & (voltage_df["i"]==i)
                    volt  = np.array( voltage_df["volt"][w] )

                    if len( xlist) > 0:
                        output = get_mean_activity_per_cycle([xlist], [ylist], theta_cycle=125)
                        aux["nspikes_mean"].append(output[0])
                        aux["nspikes_sem"].append(output[1])
                        aux["active_mean"].append(output[2])
                        aux["active_sem"].append(output[3]) 
                        aux["active_same_mean"].append(output[4])
                        aux["active_same_sem"].append(output[5])
                    else: 
                        aux["nspikes_mean"].append(np.nan)
                        aux["nspikes_sem"].append(np.nan)
                        aux["active_mean"].append(np.nan )
                        aux["active_sem"].append(np.nan)
                        aux["active_same_mean"].append(np.nan )
                        aux["active_same_sem"].append(np.nan)

                    if len(volt) > 0:
                        output = get_voltage_diff_per_cycle( time, volt )
                        aux["vmin_mean"].append(output[0])
                        aux["vmin_sem"].append(output[1])
                        aux["vmax_mean"].append(output[2])
                        aux["vmax_sem"].append(output[3])
                        aux["vdiff_mean"].append(output[4])
                        aux["vdiff_sem"].append(output[5])
                    else:
                        aux["vmin_mean"].append(np.nan)
                        aux["vmin_sem"].append(np.nan)
                        aux["vmax_mean"].append(np.nan)
                        aux["vmax_sem"].append(np.nan)
                        aux["vdiff_mean"].append(np.nan)
                        aux["vdiff_sem"].append(np.nan)

                    aux["cell"].append( cell )  
                    aux["k"].append(k) 
                    aux["i"].append(i) 
                    aux["j"].append(j)

        data_df = pd.DataFrame( aux )
        file_management.save_lzma( data_df, "data_df_"+label+"_double_scan.lzma",parent_dir = save_dir )
        print("DataFrame created and saved")

# ###############################################################################
# Ploting functions 
# ###############################################################################

def get_external_inputs_times():
    external_inputs_times = { }
    external_inputs_times["ec2_180"]    = np.round( np.arange(50+125/2.0, 8000, 125), 1)
    external_inputs_times["ec3_360"]    = np.round( np.arange(50, 8000, 125), 1 )
    external_inputs_times["sep_180"]    = np.round( np.arange(50+125/2.0, 8000, 125), 1)
    external_inputs_times["sep_360"]    = np.round( np.arange(50, 8000, 125), 1 ) 
    external_inputs_times["dg_regular"] = np.round( np.arange(50+0.375*125, 8000, 125), 1)
    external_inputs_times["dg_burst"]   = np.round( np.arange(50+0.125*125, 8000, 125), 1)

    external_inputs = { }
    external_inputs_folder = "/home/jaime/Desktop/hippocampus/external_inputs_python/10"
    external_inputs["pyr_ca3"] = file_management.load_lzma( os.path.join( external_inputs_folder, "external_inputs_pyr_ca3_0_0.lzma"))
    external_inputs_times["pyr_ca3"] = np.round( np.arange(50+80/360*125, 8000, 125 ), 1) # theoretical
    external_inputs_times["pyr_ca1"] = np.round( np.arange(50+200/360*125, 8000, 125), 1) 

    external_inputs_folder = "/home/jaime/Desktop/hippocampus/external_inputs_python/"
    external_inputs["ec3_360"] = file_management.load_lzma( os.path.join( external_inputs_folder, "external_inputs_ec3_360_0.lzma"))

    external_inputs_time_phase = {}
    for lb in ["pyr_ca3", "ec3_360"]: 
        external_inputs_time_phase[lb] = np.zeros( (len(external_inputs_times[lb]),2))
        for i, h in enumerate( external_inputs_times[lb]): 
            aux = external_inputs[lb]["tvec"]
            w = ( aux >= h-65) & ( aux < h+65)
            external_inputs_time_phase[lb][i,0] = np.median(aux[w]) 
            external_inputs_time_phase[lb][i,1] = np.std(aux[w]) 
    return external_inputs_times, external_inputs, external_inputs_time_phase 

def get_external_inputs_raster_plots( external_inputs_time, external_inputs, external_inputs_time_phase ):
    
    fig, axs = plt.subplot_mosaic([['pyr_ca3'], ['ec3_360']],constrained_layout=True, figsize=(15,5))
    for lb in ["pyr_ca3","ec3_360"]:
        axs[lb].plot( external_inputs[lb]["tvec"], external_inputs[lb]["idvec"],'o', markersize=1)
        for (h,dh) in external_inputs_time_phase[lb]: 
            axs[lb].axvline(h,color="tab:red")
            axs[lb].axvspan(h-dh, h+dh, alpha=0.4, color='tab:red')
        axs[lb].set_xlim([4000,8000])
        axs[lb].set_title(lb)
        axs[lb].grid(True)

def get_cycle_analysis_plots(data_df):
    
    n = np.max(data_df["i"])+1
    ntrials = np.max(data_df["k"])+1
    cells = data_df["cell"].unique() 
    
    x = np.arange(n)
    mosaic_list = [cells]
    fig, ax = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,3))
    for cell in cells:
        w = ( data_df["cell"] == cell ) & ( data_df["k"] == 0 )
        y, yerr = data_df["nspikes_mean"][w], data_df["nspikes_sem"][w]
        ax[cell].errorbar(x = x, y = y, yerr=yerr,marker='o') 
        ax[cell].grid(True)
        ax[cell].set_xlabel('# scan')
        ax[cell].set_ylabel('# Number of spikes')
        ax[cell].set_title(cell)

    fig, ax = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,3))
    for cell in cells:
        w = ( data_df["cell"] == cell ) & ( data_df["k"] == 0 )
        y1, y1err = data_df["vmax_mean"][w], data_df["vmax_sem"][w]
        y2, y2err = data_df["vmin_mean"][w], data_df["vmin_sem"][w]
        ax[cell].errorbar(x = x, y = y1, yerr = y1err, marker='o') 
        ax[cell].errorbar(x = x, y = y2, yerr = y2err, marker='o') 
        ax[cell].grid(True)
        ax[cell].set_xlabel('# scan')
        ax[cell].set_ylabel('Vmax,Vmin per cycle [mV]')
        ax[cell].set_title(cell)

    fig, ax = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,3))
    for cell in cells:
        w = ( data_df["cell"] == cell ) & ( data_df["k"] == 0 )
        y, y1err = data_df["vdiff_mean"][w], data_df["vdiff_sem"][w]

        ax[cell].errorbar(x = x, y = y, yerr = yerr, marker='o') 
        ax[cell].grid(True)
        ax[cell].set_xlabel('# scan')
        ax[cell].set_ylabel('Vmax-Vmin per cycle [mV]')
        ax[cell].set_title(cell)

    fig, ax = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,3))
    for cell in cells:
        w = ( data_df["cell"] == cell ) & ( data_df["k"] == 0 )
        x,y = data_df["vdiff_mean"][w], data_df["nspikes_mean"][w]
        ax[cell].plot(x,y,'o-') 
        ax[cell].grid(True)
        ax[cell].set_xlabel('Vmax-Vmin per cycle [mV]')
        ax[cell].set_ylabel('# Number of spikes')
        ax[cell].set_title(cell)
        
def get_cycle_analysis_plots_double_scan(data_df):
    n1 = np.max(data_df["i"])+1
    n2 = np.max(data_df["j"])+1
    ntrials = np.max(data_df["k"])+1
    cells = list(data_df["cell"].unique())
    ncells = len(cells) 

    x1 = np.arange(n1)
    x2 = np.arange(n2)
    bar_list = [ cell+"_bar" for cell in cells ]
    mosaic_list = [ bar_list, cells] 
    widths= [1 for i in range( ncells)]
    heights = [0.05, 1] 

    label_list = [ ["nspikes_mean", "Number of spikes: mean"], 
                   ["nspikes_sem", "Number of spikes: error"], 
                   ["vmax_mean", "Vmax per cycle"], 
                   ["vmin_mean","Vmin per cycle"], 
                   ["vdiff_mean", "Vmax - Vmin per cycle"] ]

    for label,title in label_list:
        fig, ax = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,4), 
                                    gridspec_kw={"width_ratios": widths, "height_ratios": heights})
        for cell in cells:
            w = ( data_df["cell"] == cell ) & ( data_df["k"] == 0 ) 
            df = data_df[w].pivot('i','j',label)
            x, y, z= df.columns.values, df.index.values, df.values

            plc = ax[cell].pcolormesh(x,y,z) 
            ax[cell].set_xlabel('#2 scan',fontsize=15)
            ax[cell].set_title( cell, fontsize=15)
            cb = Colorbar(ax = ax[cell+"_bar"], mappable = plc, orientation = 'horizontal', ticklocation = 'top')
            st = fig.suptitle(title, fontsize=20)
            st.set_y(0.95)
        ax[cells[0]].set_ylabel('#1 scan',fontsize=15)

def get_spiking_and_voltage_plots( runtime, voltage_df, spikes_df, external_inputs, external_inputs_time_phase, external_inputs_time,
                                 t1=2000,t2=4000):
    n = np.max(voltage_df["i"])+1
    ntrials = np.max( voltage_df["k"])+1
    cells = voltage_df["cell"].unique()

    w = (voltage_df["cell"]=="pyr") & (voltage_df["k"]==0) & (voltage_df["i"]==0)
    nn = len( voltage_df["volt"][w] ) 
    time = np.linspace(0,runtime, nn)

    colors = ["tab:blue","tab:orange","tab:green", "tab:red"]
    colors_external_inputs = ["tab:red", "tab:orange", "tab:purple","tab:pink"]

    mosaic_labels = ["spikes", "voltage"]
    mosaic_list = [ ["spikes"], ["voltage"]] 

    for i in range(0,n,2):
        fig, axs = plt.subplot_mosaic(mosaic_list,constrained_layout=True, figsize=(20,8))

        # first plot
        for ii, cell in enumerate(cells): 
            w = (spikes_df["cell"]==cell) & (spikes_df["k"]==0) & (spikes_df["i"]==i)
            stime, idd = spikes_df["time"][w], spikes_df["id"][w]
            wt = (stime >= t1) & (stime<=t2)
            axs["spikes"].plot(stime[wt], idd[wt],'o', markersize=2,label=cell,color=colors[ii])
        axs["spikes"].grid(True)
        axs["spikes"].legend(loc='upper right')
        axs["spikes"].set_title( str(i)+" 0",fontsize=20)

        # second plot
        vmin, vmax = [],[]
        for ii, cell in enumerate(cells): 
            w = (voltage_df["cell"]==cell) & (voltage_df["k"]==0) & (voltage_df["i"]==i)
            y = voltage_df["volt"][w]
            wt = (time>=t1) & (time<=t2) 
            axs["voltage"].plot(time[wt], y[wt], label=cell, color=colors[ii])
            vmin.append( np.min(y) )
            vmax.append( np.max(y) )

        vmin = np.min(vmin)
        vmax = np.max(vmax)
        dv = vmax-vmin

        x,y = external_inputs["pyr_ca3"]["tvec"], external_inputs["pyr_ca3"]["idvec"]
        y = (y-y.min())/(y-y.min()).max()*dv+vmin

        axs["voltage"].plot(x,y,'o', markersize=1) 
        axs["voltage"].grid(True)
        axs["voltage"].legend(loc='upper right')
        axs["voltage"].set_title( str(i)+" 0",fontsize=20)

        # common to both plots:
        for ax_lb in mosaic_labels:

            for lb in ["pyr_ca3","ec3_360"]:
                for (h,dh) in external_inputs_time_phase[lb]: 
                    axs[ax_lb].axvline(h,color="black", linestyle="--")
                    if lb == "pyr_ca3":
                        axs[ax_lb].axvspan(h-dh, h+dh, alpha=0.5, color='gray')

            for yh in external_inputs_time["pyr_ca1"]: 
                axs[ax_lb].axvline(yh, color='red', linewidth=4, linestyle="--")
                axs[ax_lb].set_xlim([t1,t2])

            #axs[ax_lb].set_xlim([t1,t2])
            
def get_volt_traces_plots(runtime, voltage_df, external_inputs, external_inputs_time_phase, 
                          external_inputs_time,t1=2000,t2=4000):
    
    n = np.max(voltage_df["i"])+1
    ntrials = np.max( voltage_df["k"])+1
    cells = voltage_df["cell"].unique()
    
    w = (voltage_df["cell"]=="pyr") & (voltage_df["k"]==0) & (voltage_df["i"]==0)
    nn = len( voltage_df["volt"][w] ) 
    time = np.linspace(0,runtime, nn)
    
    vmap = []
    # first plot 
    plt.figure(figsize=(20,6))
    for i in range(n):
        w = (voltage_df["cell"] == "pyr") & (voltage_df["k"]==0) & (voltage_df["i"]==i)
        y = voltage_df["volt"][w]
        vmap.append(y)
        plt.plot(time, -15*i+y,label="pyr", color="tab:blue", alpha=0.5+0.02*i)
    for (h,dh) in external_inputs_time_phase["pyr_ca3"]:
        plt.axvline(h, color="black")
    plt.xlim([t1,t2])
    plt.grid(True)
    
    # second plot
    vmap = np.array(vmap)
    plt.figure(figsize=(20,6),dpi=100)
    plt.pcolormesh(time, np.arange(n),vmap,shading="auto" )
    for (h,dh) in external_inputs_time_phase["pyr_ca3"]:
        plt.axvline(h, color="black")
    for h in external_inputs_time["pyr_ca1"]:
        plt.axvline(h, color='red', linewidth=4, linestyle="--")
    plt.xlim([t1,t2])
    plt.colorbar()
    plt.xlabel("Time [ms]",fontsize=15)
    plt.ylabel("Delay [ms]",fontsize=15)
    
def get_spiking_distribution_plots( runtime, spikes_df, show_histogram = True, show_spiking_cck_pyr = True, show_spiking_per_cell = True, sigma=0.5, ylabel="# scan", title=""): 
    n = np.max(spikes_df["i"])+1
    ntrials = np.max( spikes_df["k"])+1
    cells = spikes_df["cell"].unique()

    w = (spikes_df["cell"]=="pyr") & (spikes_df["k"]==0) & (spikes_df["i"]==0)

    binsize = 5 
    t0=50
    time_binned = np.arange(t0,runtime,binsize)
    nbins = len(time_binned)

    rate_cycles = dict.fromkeys(cells)
    rate_mean = dict.fromkeys(cells)
    rate_sem = dict.fromkeys(cells)
    ncycles = int((runtime-50)/125)
    
    phase = np.linspace(0,360,25)

    for cell in cells:
        rate_mean[cell]   = [ [] for k in range(n) ]
        rate_sem[cell]    = [ [] for k in range(n) ]
        rate_cycles[cell] = [ [] for k in range(n) ]

        for i in range(n):
            w = (spikes_df["cell"] == cell) & (spikes_df["k"]==0) & (spikes_df["i"]==i)
            y, x = np.histogram( np.array(spikes_df["time"][w]),bins=time_binned)
            x = x[1:]-np.diff(x)[0]/2

            rate = []
            for l in range(ncycles):
                t1 = t0+l*125 
                t2 = t0+(l+1)*125
                rate.append(y[(x>=t1) & (x<t2)])

            rate = np.array(rate)
            rate_cycles[cell][i] = np.concatenate(rate)
            rate_mean[cell][i] = np.mean(rate[2:-1],axis=0)
            rate_sem[cell][i]  = np.std(rate[2:-1],axis=0)/np.sqrt(ncycles-3)

        rate_mean[cell] = np.array(rate_mean[cell]) 
        rate_sem[cell]  = np.array(rate_sem[cell])

    for cell in cells:
        rate_cycles[cell] = np.array( rate_cycles[cell] )

    # first plot
    if show_histogram:
        plt.figure(figsize=(20,6),dpi=100)
        plt.pcolormesh(time_binned[time_binned<=t2], np.arange(n), rate_cycles["pyr"],shading="auto" )
        plt.colorbar()
        plt.xlabel("Time [ms]",fontsize=15)
        plt.ylabel("# scan",fontsize=15)
        plt.title("Spiking histogram")

    # second plot 
    if show_spiking_cck_pyr:
        fig, ax = plt.subplots(figsize=(6,5),dpi=100) 
        for i,(r,rs) in enumerate( zip(rate_mean["pyr"][:-1], rate_sem["pyr"][:-1])):
            ax.plot(phase,r,alpha=np.tanh(i/20), color="tab:blue")
        ax.grid(True)
        ax.set_ylabel("# spikes pyr", color="tab:blue")
        ax.set_xlabel("Theta phase")
        ax.set_title(" Phase relation")
        ax2 = ax.twinx()
        for i, (r, rs) in enumerate( zip(rate_mean["cck"], rate_sem["cck"])) :
            ax2.plot(phase,r, alpha=0.05*i, color="tab:orange")
        ax2.set_ylabel("# spikes cck", color="tab:orange")

    # thir plot
    if show_spiking_per_cell:
        cells = list(spikes_df["cell"].unique())
        ncells = len(cells) 
        bar_list = [ cell+"_bar" for cell in cells ]
        mosaic_list = [ bar_list, cells] 
        widths= [1 for i in range(ncells) ]
        heights = [0.05, 1] 

        fig, axs = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,4), 
                                     gridspec_kw={"width_ratios": widths, "height_ratios": heights})
        for cell in cells: 
            data = scipy.ndimage.gaussian_filter( rate_mean[cell] , sigma = sigma)
            pcl = axs[cell].pcolormesh( phase, np.arange(n), data, shading="auto")
            axs[cell].set_title(cell,fontsize=15)
            axs[cell].set_xlabel("Theta phase",fontsize=15)
            cb = Colorbar(ax = axs[cell+"_bar"], mappable = pcl, orientation = 'horizontal', ticklocation = 'top')
        axs[cells[0]].set_ylabel(ylabel,fontsize=15)
        st = fig.suptitle("Spiking pattern "+title, fontsize = 20)
        st.set_y(0.95)
    return time_binned, rate_cycles


def get_spiking_distribution_plots_alternative( runtime, spikes_df, ylabel="# scan", title=""): 
    n = np.max(spikes_df["i"])+1
    ntrials = np.max( spikes_df["k"])+1
    cells = spikes_df["cell"].unique()

    w = (spikes_df["cell"]=="pyr") & (spikes_df["k"]==0) & (spikes_df["i"]==0)

    binsize = 5 
    t0 = 50
    theta_rythm = 125
    time_binned = np.arange(t0,runtime,binsize)
    nbins = len(time_binned)

    rate_cycles = dict.fromkeys(cells)
    rate_mean = dict.fromkeys(cells)
    rate_sem = dict.fromkeys(cells)
    rate_mean_fit = dict.fromkeys(cells)
    rate_sup_fit = dict.fromkeys(cells)
    rate_inf_fit = dict.fromkeys(cells)
    
    ncycles = int((runtime-t0)/theta_rythm)
    phase = np.linspace(0,360,25)

    for cell in cells:
        rate_mean[cell]   = [ [] for k in range(n) ]
        rate_sem[cell]    = [ [] for k in range(n) ]
        rate_cycles[cell] = [ [] for k in range(n) ]
        rate_mean_fit[cell] = [ [] for k in range(n) ]
        rate_sup_fit[cell] = [ [] for k in range(n) ]
        rate_inf_fit[cell] = [ [] for k in range(n) ]

        for i in range(n):
            w = (spikes_df["cell"] == cell) & (spikes_df["k"]==0) & (spikes_df["i"]==i)
            
            y, x = np.histogram( np.array(spikes_df["time"][w]),bins=time_binned)
            x = x[1:]-np.diff(x)[0]/2
            
            ny = len(y)
            ncycles = int( ny/(theta_rythm/binsize) )
            nnew = int(ncycles*theta_rythm/binsize)
            yy = np.reshape( y[:nnew], (ncycles, int(theta_rythm/binsize)) )
            phase = np.linspace(0, 360, int(theta_rythm/binsize)) 
            
            rate_mean[cell][i] = np.mean(yy, axis=0)
            rate_sem[cell][i] = np.std(yy, axis=0)/np.sqrt(ncycles) 
            
            phase,ym,y1,y2 = smoothing( phase, rate_mean[cell][i], rate_sem[cell][i], s=15)
            rate_mean_fit[cell][i] = ym
            rate_sup_fit[cell][i] = y1
            rate_inf_fit[cell][i] = y2 
            
        rate_mean[cell] = np.array(rate_mean[cell]) 
        rate_sem[cell]  = np.array(rate_sem[cell])
        rate_mean_fit[cell] = np.array(rate_mean_fit[cell])
        rate_sup_fit[cell] = np.array(rate_sup_fit[cell])
        rate_inf_fit[cell] = np.array(rate_inf_fit[cell])
    
    cells = list(voltage_df["cell"].unique())
    ncells = len(cells) 
    bar_list = [ cell+"_bar" for cell in cells ]
    mosaic_list = [ bar_list, cells] 
    widths= [1 for i in range(ncells) ]
    heights = [0.05, 1] 
    
    fig, axs = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,4), 
                                 gridspec_kw={"width_ratios": widths, "height_ratios": heights})
    for cell in cells: 
        data = scipy.ndimage.gaussian_filter( rate_mean_fit[cell] , sigma = sigma)
        pcl = axs[cell].pcolormesh( phase, np.arange(n), data, shading="auto")
        axs[cell].set_title(cell,fontsize=15)
        axs[cell].set_xlabel("Theta phase",fontsize=15)
        cb = Colorbar(ax = axs[cell+"_bar"], mappable = pcl, orientation = 'horizontal', ticklocation = 'top')
    axs[cells[0]].set_ylabel(ylabel,fontsize=15) 
    st = fig.suptitle("Spiking pattern "+title, fontsize = 20)
    st.set_y(0.95)
    
    return phase, rate_mean, rate_mean_fit


def get_spiking_distribution_kde_plots( runtime, spikes_df, ylabel="# scan", title=""): 
    n = np.max(spikes_df["i"])+1
    ntrials = np.max( spikes_df["k"])+1
    cells = spikes_df["cell"].unique()

    w = (spikes_df["cell"]=="pyr") & (spikes_df["k"]==0) & (spikes_df["i"]==0)
    
    t0 = 50
    theta_rythm = 125
    rate_mean = dict.fromkeys(cells)
    phase = np.linspace(0,360,101)
    for cell in cells:
        rate_mean[cell]   = [ [] for k in range(n) ]

        for i in range(n):
            w = (spikes_df["cell"] == cell) & (spikes_df["k"]==0) & (spikes_df["i"]==i)
            
            y = np.array(spikes_df["time"][w])-t0
            y = np.mod(y[y>4*theta_rythm],theta_rythm)/theta_rythm*360
            
            kde = KernelDensity(bandwidth=3.0, kernel='gaussian')
            kde.fit(y[:, None])
            # score_samples returns the log of the probability density
            logprob = kde.score_samples(phase[:, None])           
            rate_mean[cell][i] = np.exp(logprob)
            
        rate_mean[cell] = np.array(rate_mean[cell]) 
    
    cells = list(spikes_df["cell"].unique())
    ncells = len(cells) 
    bar_list = [ cell+"_bar" for cell in cells ]
    mosaic_list = [ bar_list, cells] 
    widths= [1 for i in range(ncells) ]
    heights = [0.05, 1] 
    
    fig, axs = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,4), 
                                  gridspec_kw={"width_ratios": widths, "height_ratios": heights})
    for cell in cells: 
        pcl = axs[cell].pcolormesh( phase, np.arange(n), rate_mean[cell], shading="auto")
        axs[cell].set_title(cell,fontsize=15)
        axs[cell].set_xlabel("Theta phase",fontsize=15)
        cb = Colorbar(ax = axs[cell+"_bar"], mappable = pcl, orientation = 'horizontal', ticklocation = 'top')
    axs[cells[0]].set_ylabel(ylabel,fontsize=15) 
    st = fig.suptitle("Spiking pattern "+title, fontsize = 20)
    st.set_y(0.95)
    
    return phase, rate_mean

def get_mebrane_potential_distribution_plots( runtime, voltage_df ,sigma=0.5, ylabel="# scan", title=""): 
    n = np.max(voltage_df["i"])+1
    ntrials = np.max( voltage_df["k"])+1
    cells = voltage_df["cell"].unique()

    t0 = 50
    theta_rythm = 125
    ncycles = int((runtime-t0)/theta_rythm)
    w = (voltage_df["cell"] == "pyr") & (voltage_df["k"]==0) & (voltage_df["i"]==0)
    v0 = np.array(voltage_df["volt"])[t0:]
    nv = len(v0)
    nvnew = int( theta_rythm*ncycles ) 
    phase = np.linspace(0,360,theta_rythm)
    vmean = dict.fromkeys(cells)
    vsem  = dict.fromkeys(cells)
                          
    for cell in cells:
        vmean[cell] = [ [] for i in range(n) ]
        vsem[cell]  = [ [] for i in range(n) ] 

        for i in range(n):
            w = (voltage_df["cell"] == cell) & (voltage_df["k"]==0) & (voltage_df["i"]==i)
            v = np.array(voltage_df["volt"][w])[t0:]
            vv = np.reshape( v[:nvnew], (ncycles, theta_rythm))
            
            vmean[cell][i] = np.mean(vv, axis=0)
            vsem[cell][i] = np.std(vv, axis=0)/np.sqrt(ncycles) 
                
        vmean[cell] = np.array(vmean[cell]) 
        vsem[cell]  = np.array(vsem[cell])
                          
    cells = list(voltage_df["cell"].unique())
    ncells = len(cells) 
    bar_list = [ cell+"_bar" for cell in cells ]
    mosaic_list = [ bar_list, cells] 
    widths= [1 for i in range(ncells) ]
    heights = [0.05, 1] 
    
    fig, axs = plt.subplot_mosaic(mosaic_list, constrained_layout=True, figsize=(12,4), 
                                 gridspec_kw={"width_ratios": widths, "height_ratios": heights})
    for cell in cells: 
        data = scipy.ndimage.gaussian_filter( vmean[cell] , sigma = sigma)
        pcl = axs[cell].pcolormesh( phase, np.arange(n), data, shading="auto")
        axs[cell].set_title(cell,fontsize=15)
        axs[cell].set_xlabel("Theta phase",fontsize=15)
        cb = Colorbar(ax = axs[cell+"_bar"], mappable = pcl, orientation = 'horizontal', ticklocation = 'top')
    
    axs[cells[0]].set_ylabel(ylabel,fontsize=15) 
    st = fig.suptitle("Membrane potential "+title, fontsize = 20)
    st.set_y(0.95)
    return phase, vmean, vsem

############################################################################################################
# Spectral analysis 
############################################################################################################

def multitapering(x,y, t0=50,NW = 2.5, k=4):
    '''
    Multitaper analysis
    '''
    lowcut  = [25.0,  6]
    highcut = [100.0, 10]
    dt = np.diff(x)
    
    Power_theta = []
    Power_gamma = []

    nsamples = len(y)
    # classical FFT
    yf = np.fft.fft(y-y.mean())
    yf = abs(yf)**2 / nsamples * dt
    fr = np.fft.freqfft(n=nsamples, d=1.0/dt)[:nsamples//2]
    df = np.diff(xf)[0]
    #NW=2.5 #k=4
    [tapers, eigen] = dpss(nsamples, NW, k)
    Sk_complex, weights, eigenvalues=pmtm(y-y.mean(), e=eigen, v=tapers, NFFT=nsamples, show=False)
    Sk = abs(Sk_complex)**2
    Sk = np.mean(Sk*np.transpose(weights), axis=0)*dt

    # comparing total power
    theta_band = (xf>=lowcut[0]) & (xf<=highcut[0])
    gamma_band = (xf>=lowcut[1]) & (xf<=highcut[1])
    Sk = Sk[0:nsamples//2]
    '''
    Power_theta.append(np.mean(Sk[theta_band])*df)
    Power_beta.append(np.mean(Sk[beta_band])*df)
    Power_gamma.append(np.mean(Sk[gamma_band])*df)

    theta_max.append( xf[theta_band][np.argmax(Sk[theta_band])] )
    beta_max.append(  xf[beta_band][np.argmax(Sk[beta_band])] )
    gamma_max.append( xf[gamma_band][np.argmax(Sk[gamma_band])] )

    theta_max_amp.append( np.max(Sk[theta_band]) )
    beta_max_amp.append(  np.max(Sk[beta_band]) )
    gamma_max_amp.append( np.max(Sk[gamma_band]) )
    
    freq_spectrum.append(xf)
    Power_spectrum.append(Sk)
    power = [Power_theta,  Power_beta,   Power_gamma]
    fmax  = [theta_max,    beta_max,     gamma_max]
    amax  = [beta_max_amp, beta_max_amp, gamma_max_amp]
    '''
    #return freq_spectrum, Power_spectrum, power, fmax, amax
    return fr, Sk, yf