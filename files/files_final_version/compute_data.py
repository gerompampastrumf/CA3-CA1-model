import os 
import sys 
import compute_psd as cp  
import compute_phases_from_spikes as cpfs 
import compute_theta_gamma_coupling as ctgc


path = '/home/jaime/Desktop/hippocampus/external_inputs/'

path_ca3 = os.path.join(path,'baseline')
path_ca1 = path_ca3#os.path.join(path,'baseline') # they will be different in a specific cases.. where actually ca3 would be already computed

print("Computing PSDs CA3 ...")
print("... from voltages")
cp.compute_psd_volt_ca3(path_ca3)
print("... from LFPs")
cp.compute_psd_lfp_ca3(path_ca3)
print("... from ICAs")
cp.compute_psd_ica_ca3(path_ca3)
print("... from synaptic currents")
cp.compute_psd_synaptic_currents_ca3(path_ca3)

print("Computing PSDs CA1 ...")
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

# Compute phases
# print("Computing phases CA3 ...")
# cpfs.compute_phases_from_spikes_ca3(path_ca3)
# print("Computing phases CA1 ...")
# cpfs.compute_phases_from_spikes_ca1(path_ca1)
# print(" ... done")
# print(" ")

# Compute theta-gamma coupling
print("Computing theta-gamma coupling CA3 ...")
# print("... from voltages")
# ctgc.compute_amplitude_coupling_from_volt(path_ca3, input_filename="volt_pyr_ca3.lzma",output_filename="volt_pyr_ca3_amplitude_coupling.lzma")
# print("... from LFPs")
# ctgc.compute_amplitude_coupling_from_lfp(path_ca3, input_filename="lfp_ca3.lzma", output_filename="lfp_ca3_amplitude_coupling.lzma")
# print("... from ICAs")
# ctgc.compute_amplitude_coupling_from_ica(path_ca3, input_filename="ica_ca3.lzma", output_filename="ica_ca3_amplitude_coupling.lzma")
print("... from synaptic currents")
# ctgc.compute_amplitude_coupling_from_synaptic_current_ca3(path_ca3)

print("Computing theta-gamma coupling CA1 ...")
# print("... from voltages")
# ctgc.compute_amplitude_coupling_from_volt(path_ca1,input_filename="volt_pyr_ca1.lzma",output_filename="volt_pyr_ca1_amplitude_coupling.lzma")
# print("... from LFPs")
# ctgc.compute_amplitude_coupling_from_lfp(path_ca1, input_filename="lfp_ca1.lzma", output_filename="lfp_ca1_amplitude_coupling.lzma")
# print("... from ICAs")
# ctgc.compute_amplitude_coupling_from_ica(path_ca1, input_filename="ica_ca1.lzma", output_filename="ica_ca1_amplitude_coupling.lzma")
# print("... from synaptic currents")
# ctgc.compute_amplitude_coupling_from_synaptic_current_ca1(path_ca1)
print(" ... done")
print(" ")
