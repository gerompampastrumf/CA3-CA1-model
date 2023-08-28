import os 
import sys 
import compute_psd as cp  
import compute_phases_from_spikes as cpfs 
import compute_theta_gamma_coupling as ctgc

# compute data of the first table:
folders = ["baseline","baseline_nobasketCA3","baseline_noolmCA3","baseline_noinhibitionCA3",
           "baseline_nobasketCA1","baseline_noolmCA1","baseline_nocckCA1",
           "baseline_nobasketnocckCA1","baseline_noolmnocckCA1","baseline_nobasketnoolmCA1",
           "baseline_noinhibitionCA1",]

path = '/home/jaime/Desktop/hippocampus/external_inputs/'

path_ca3 = os.path.join(path, folders[int(sys.argv[1])])
path_ca1 = path_ca3#os.path.join(path,'baseline') # they will be different in a specific cases.. where actually ca3 would be already computed

print("Computing PSDs CA3 ...")
try:
    print("... from voltages")
    cp.compute_psd_volt_ca3(path_ca3)
    print("... from LFPs")
    cp.compute_psd_lfp_ca3(path_ca3)
    print("... from ICAs")
    cp.compute_psd_ica_ca3(path_ca3)
    print("... from synaptic currents")
    cp.compute_psd_synaptic_currents_ca3(path_ca3)
except FileNotFoundError: 
    print("... no CA3 files in this folder")

print("Computing PSDs CA1 ...")
try:
    print("... from voltages")
    cp.compute_psd_volt_ca1(path_ca1)
    print("... from LFPs")
    cp.compute_psd_lfp_ca1(path_ca1)
    print("... from ICAs")
    cp.compute_psd_ica_ca1(path_ca1)
    print("... from synaptic currents")
    cp.compute_psd_synaptic_currents_ca1(path_ca1)
    print(" ... done")
    print(" ")
except FileNotFoundError: 
    print("... no CA1 files in this folder")

# Compute phases
print("Computing phases CA3 ...")
try:
    cpfs.compute_phases_from_spikes_ca3(path_ca3)
except FileNotFoundError: 
    print("... no CA3 files in this folder")
print("Computing phases CA1 ...")
try:
    cpfs.compute_phases_from_spikes_ca1(path_ca1)
    print(" ... done")
except FileNotFoundError: 
    print("... no CA1 files in this folder")
print(" ")

# Compute theta-gamma coupling
print("Computing theta-gamma coupling CA3 ...")
try:
    print("... from voltages")
    ctgc.compute_amplitude_coupling_from_volt(path_ca3, input_filename="volt_pyr_ca3.lzma",output_filename="volt_pyr_ca3_amplitude_coupling.lzma")
    print("... from LFPs")
    ctgc.compute_amplitude_coupling_from_lfp(path_ca3, input_filename="lfp_ca3.lzma", output_filename="lfp_ca3_amplitude_coupling.lzma")
    print("... from ICAs")
    ctgc.compute_amplitude_coupling_from_ica(path_ca3, input_filename="ica_ca3.lzma", output_filename="ica_ca3_amplitude_coupling.lzma")
    print("... from synaptic currents")
    ctgc.compute_amplitude_coupling_from_synaptic_current_ca3(path_ca3)
except FileNotFoundError: 
    print("... no CA3 files in this folder")

print("Computing theta-gamma coupling CA1 ...")
try:
    print("... from voltages")
    ctgc.compute_amplitude_coupling_from_volt(path_ca1,input_filename="volt_pyr_ca1.lzma",output_filename="volt_pyr_ca1_amplitude_coupling.lzma")
    print("... from LFPs")
    ctgc.compute_amplitude_coupling_from_lfp(path_ca1, input_filename="lfp_ca1.lzma", output_filename="lfp_ca1_amplitude_coupling.lzma")
    print("... from ICAs")
    ctgc.compute_amplitude_coupling_from_ica(path_ca1, input_filename="ica_ca1.lzma", output_filename="ica_ca1_amplitude_coupling.lzma")
    print("... from synaptic currents")
    ctgc.compute_amplitude_coupling_from_synaptic_current_ca1(path_ca1)
    print(" ... done")
    print(" ")
except FileNotFoundError:
    print("... no CA1 files in this folder")

