'''
Date 11/05/2023

'''
import warnings
warnings.filterwarnings("ignore")
from neuron import h
import numpy as np
import sys
import os
from neuron.units import ms, mV
import time as tm
sys.path.append('/home/dimitrios/Neurons/CA1model/final_files/')
#from neurons import *
from network_hippocampus_ca3 import *
from measurements import *
from LFPsimpy import LfpElectrode
from parameters import *
import file_management
from functions import *

'''###########################################################################
                    Parameters
###########################################################################'''
tick = tm.time()

# inputs values from bash
inputs_argvs = [] 
column_labels = []
for argv in sys.argv[1:]:
    inputs_argvs.append(int(argv))

if len(inputs_argvs) > 2:
    iseed  = inputs_argvs[0] 
    ibseed = inputs_argvs[1] 
    # must be added the other variables in case they are used
    # ... = inputs[2]
    # ... = inputs[3]
    # ...
    number_of_argvs = len(sys.argv[1:])
    argvs=""
    for i in range(number_of_argvs-1):
        argvs += sys.argv[1:][i]+'_'
    argvs+=sys.argv[-1]

    column_labels.append("iseed")
    column_labels.append("ibseed")
    # column_labels ...

elif len(inputs_argvs) == 2: 
    # in case no inputs are given
    iseed  = inputs_argvs[0] # external inputs basal seed
    # control the seed for the external inputs generation (online generation, not implemented)
    # and select the trial of the external input (offline generation, implemented)
    ibseed = inputs_argvs[1] # background noise basal seed
    # control the seed for the background noise generation (online generation)
        # and would also select the trial of the background noise (offline generation, not implemented)
    number_of_argvs = len(sys.argv[1:])
    argvs=""
    for i in range(number_of_argvs-1):
        argvs += sys.argv[1:][i]+'_'
    argvs+=sys.argv[-1]

    column_labels.append("iseed")
    column_labels.append("ibseed")
else: 
    print("provide at least two inputs")
    sys.exit()

h.dt = 0.1 # time resolution of the simulation
h.celsius = 34.0 # temperature of the simulation
simulation_time = 30000.0 # total time of the simulation
time_resolution = 2.0 # time resolution of the measurements

DoMakeNoise = True # Control the external noise
MakeNetStim = True # Online generation of the background noise 
DoMakeExternalInputs = True # Control the external inputs
MakeCellStim         = False # Online generation of the external inputs

create_ca1_network  = False # Create the CA1 network as well
record_all_synapses = True # Record all the synapses
record_lfp          = True # Record the LFP and transmembrane currents

# synapses to be recorded
synlist = [ "Adend3GABA_olm","Adend3AMPA_ec2360","Adend3NMDA_ec2360", "Adend3GABA_noise", "Adend3AMPA_noise","Adend3AMPA_ec2180","Adend3NMDA_ec2180",
            "Adend1AMPA_dgreg", "Adend1NMDA_dgreg","Adend1AMPA_dgburst","Adend1NMDA_dgburst",
            "somaGABA_bas", "somaAMPA_noise", "somaGABA_noise",
            "BdendAMPA_pyr","BdendNMDA_pyr"]

# membrane potential to be recorded
record_mp = {"Bdend": False, "soma": True, "Adend1": False, "Adend2": False, "Adend3": False}

# what data to save (be consistent to what you record)
save_data_volt   = True
save_data_spikes = True 
# these are the ones that require consistency what previously declared
save_data_syn    = True
save_data_lfp    = True
save_data_ica    = True

# folders for saving
# output of ca3 will be stored in the same folder as the external inputs
current_folder = os.getcwd()
inputs_folder = '/home/dimitrios/Neurons/CA1model/baselineCA3/external_inputs'
save_folder   = '/home/dimitrios/Neurons/CA1model/baselineCA3/'
#save_folder = os.path.join(current_folder,"test_data")

# file_name = __file__
# file_name = file_name.replace(current_folder,"") 
# file_name = file_name.split('_')[-1][:-3] 
# save_dir = os.path.join(current_folder, file_name)

def unify_data(variable, cells, nr=1):
    f''' unify data from all cells in the network: membrane voltage, synaptic currents, spikes, etc '''
    local_data = { cell.id: list( np.round( cell.__dict__[variable].to_python(),nr)) for cell in cells }
    all_data = pc.py_alltoall( [local_data] + [None] * (pc.nhost() - 1) )
    return all_data

'''###########################################################################
                   Parameters optimized (to optimize)
###########################################################################'''
dg_reg_pyr_sc = np.linspace(0,2.5,11)[6]
dg_bur_pyr_sc = np.linspace(0,1,11)[7] 
olm_pyr_scale = np.linspace(0,1,11)[4]
pyr_bas_scale = np.linspace(0,2,11)[2] 
bas_pyr_scale = np.linspace(0,2,11)[4]
bas_bas_scale = np.linspace(0,2,11)[1]

weights_neurons_neurons["bas_ca3_to_bas_ca3"] = [[ 4.05e-3 *bas_bas_scale]]
weights_neurons_neurons["pyr_ca3_to_bas_ca3"] = [[ 2.341e-3*pyr_bas_scale, 1.5*1.38e-3*pyr_bas_scale]]
weights_neurons_neurons["bas_ca3_to_pyr_ca3"] = [[ 0.3175e-3 *bas_pyr_scale]]
weights_neurons_neurons["olm_ca3_to_pyr_ca3"] = [[ 57.6e-3*olm_pyr_scale]]  # [0.8*72.0e-3] #*0.9

weights_noise_neurons["bas_ca3"]["somaAMPA_noise"]   = 0.005e-3
weights_noise_neurons["pyr_ca3"]["somaAMPA_noise"]   = 0.0125e-3
weights_noise_neurons["pyr_ca3"]["Adend3AMPA_noise"] = 0.0125e-3

weights_inputs_neurons["sep_180_to_bas_ca3"] = [[6.4e-4]]
weights_inputs_neurons["sep_360_to_olm_ca3"] = [[3.2e-4]]

weights_inputs_neurons["dg_regular_to_pyr_ca3"] = [[ 1.44e-4*dg_reg_pyr_sc, 0.6e-4*dg_reg_pyr_sc ]]
weights_inputs_neurons["dg_regular_to_bas_ca3"] = [[ 0.9e-4,  0.5e-4 ]]
weights_inputs_neurons["dg_burst_to_pyr_ca3"]   = [[ 1.44e-4*dg_bur_pyr_sc, 0.6e-4*dg_bur_pyr_sc ]]
weights_inputs_neurons["dg_burst_to_bas_ca3"]   = [[ 0.9e-4,  0.5e-4 ]]

# function for online changes 
# h.load_file("stdrun.hoc")
# def fi():
#     for i in range(0,int(simulation_time),100):
#         h.cvode.event(i, "print " + str(i))

# def change_bursting():
#     w_new = net.burst_basal_ncl_[0].weight[0]
#     for nc in net.burst_var_ncl_:
#         nc.weight[0] = w_new
#     for nc in net.burst_basal_ncl_:
#         nc.weight[0] = 0.0

'''############################################################################
                                The network
#############################################################################'''

net = Network( weights_inputs_neurons  = weights_inputs_neurons,
               nsyns_inputs_neurons    = nsyns_inputs_neurons,
               syn_inputs_neurons      = syn_inputs_neurons,
               weights_neurons_neurons = weights_neurons_neurons,
               nsyns_neurons_neurons   = nsyns_neurons_neurons,
               syn_neurons_neurons     = syn_neurons_neurons,
               weights_noise_neurons   = weights_noise_neurons,
               delay_neurons_neurons   = delay_neurons_neurons,
               delay_inputs_neurons    = delay_inputs_neurons,
               iseed=3421*(2*iseed+1),  # seed for external inputs, this parameter will not be used (only for "online" generation)
               bseed=27*(2*ibseed+1),   # seed for background inputs    
               create_ca1_network = create_ca1_network,
               DoMakeNoise = DoMakeNoise,
               DoMakeExternalInputs = DoMakeExternalInputs,
               MakeCellStim = MakeCellStim,
               MakeNetStim  = MakeNetStim,
               inputs_folder = inputs_folder,
               external_inputs_iseed = iseed, # !!!
               n_pyr_ca3 = 800,
               n_bas_ca3 = 100,
               n_olm_ca3 = 30,
               n_pyr_ca1 = 800,
               n_bas_ca1 = 100,
               n_olm_ca1 = 30,
               n_cck_ca1 = 70,
               nseg_pyr  = 3,
               record_mp = record_mp,
               resolution = time_resolution)
               
net.set_noise_inputs(simulation_time) # set background noise and external input
if net.MakeNetStim:
    net.init_NetStims()  # init rngs of background
if net.MakeCellStim:
    if not net.MakeCellStim:
        net.init_CellStims() # init rngs of external inputs

'''###########################################################################
                     external inputs vecstims
############################################################################'''

inputs, vecstims = [],[]
net.inputs_list = []
net.vsl_ = []

if net.DoMakeExternalInputs:
    print("external inputs with vecstims")
    print(net.external_inputs_data.keys())
    for key in net.external_inputs_data.keys():
        if key.endswith("ca3"):
            print(key)
            idv      = net.external_inputs_data[key]["idv"]
            spikes   = net.external_inputs_data[key]["spikes"]
            trg      = net.external_inputs_data[key]["population"]
            syn_list = net.external_inputs_data[key]["synapses"]
            delays   = net.external_inputs_data[key]["delays"]
            w_list   = net.external_inputs_data[key]["weights"]
            conn     = net.external_inputs_data[key]["connectivity"]

            for k,conn_ in enumerate(conn):
                for post_id, all_pre in enumerate(conn_):
                    net.inputs_list.append([ ])
                    for j, pre_id in enumerate(all_pre):
                        idvs = np.where(idv==pre_id)[0]
                        net.inputs_list[-1].append( np.sort(spikes[idvs]))
                        inputs.append(h.Vector( np.sort(spikes[idvs]) ))
                        vecstims.append(h.VecStim())
                        vecstims[-1].play(inputs[-1])

                        spike_list = vecstims[-1]
                        net.vsl_.append(spike_list)
                        for syn,w in zip(syn_list[k], w_list[k]):
                            net.ncl_.append(h.NetCon(spike_list, trg.cell[post_id].__dict__[syn].syn, 0, delays[k][j,post_id], w))

        if net.create_ca1_network & key.endswith("ca1"):
            print(key)
            idv      = net.external_inputs_data[key]["idv"]
            spikes   = net.external_inputs_data[key]["spikes"]
            trg      = net.external_inputs_data[key]["population"]
            syn_list = net.external_inputs_data[key]["synapses"]
            delays   = net.external_inputs_data[key]["delays"]
            w_list   = net.external_inputs_data[key]["weights"]
            conn     = net.external_inputs_data[key]["connectivity"]

            for post_id, all_pre in enumerate(conn):
                net.inputs_list.append([ ])
                for j, pre_id in enumerate(all_pre):
                    idvs = np.where(idv==pre_id)[0]
                    net.inputs_list[-1].append( np.sort(spikes[idvs]) )
                    inputs.append(h.Vector( np.sort(spikes[idvs])) )
                    vecstims.append(h.VecStim())
                    vecstims[-1].play(inputs[-1])

                    spike_list = vecstims[-1]
                    net.vsl_.append(spike_list)
                    for syn,w in zip(syn_list, w_list):
                        net.ncl_.append(h.NetCon(spike_list, trg.cell[post_id].__dict__[syn].syn, 0, delays[j,post_id], w))

    tock = tm.time()
    diff = tock-tick
    horas = int(diff/3600)
    diff  = diff-horas*3600
    mint  = int(diff/60)
    diff  = diff-mint*60
    seg   = int(diff)

    print('Vecstim done after:', horas,'h', mint, 'min', seg, 's')
''' #######################################################################'''

pc = h.ParallelContext()
pc.set_maxstep(100*ms)
t = h.Vector().record(h._ref_t)
h.celsius = 34.0
# h.steps_per_ms = 1/h.dt
h.dt = 0.1

if record_all_synapses:
    # this record only the specific synapses in the ca3 pyramidal cells
    # there is a function called record_all_synapses() in the net class that record all synapses in the network
    # but the synlist should be provided first to the neurons.
    for cell in net.pyr_ca3.cell:
        cell.syn_list = synlist
        cell.record_synapses()

#Record LFP
if record_lfp:
    # electrode_y_coords = np.arange(0,850,50)
    zcoords = [-100,10,85,235, 385]
    electrode = {}
    electrode["ca3"] = []

    for i, z in enumerate(zcoords): 
        electrode["ca3"].append( LfpElectrode(x=25.0, y=25.0, z=z, sampling_period=time_resolution, neuron_type = "Pyramidal CA3"))
    
    if net.create_ca1_network:
        electrode["ca1"] = []
        for i,z in enumerate(zcoords):
            # this x position should be corrected, however we're not going to simulate the whole network
            electrode["ca1"].append( LfpElectrode(x=1e3+25.0, y=25.0, z=z, sampling_period=time_resolution, neuron_type="Pyramidal CA1"))
      
# h.tstop = simulation_time
# h.stdinit()
h.init()
h.finitialize()
h.fcurrent()
h.frecord_init()
h.tstop=simulation_time
h.dt = 0.1
h.celsius = 34
pc.psolve(simulation_time*ms) # simulation running starts

tock = tm.time()
diff = tock-tick
horas = int(diff/3600)
diff  = diff-horas*3600
mint  = int(diff/60)
diff  = diff-mint*60
seg   = int(diff)

print("Simulation done after", horas,'h', mint, 'min', seg, 's')
print("-------------------------------------------")
print("Preparing the data to be saved") 

t = np.array(t.to_python())
all_spikes, all_volt, lfp, ica = {},{},{},{}

# only saving the data if the simulation was run in the master node

if pc.id() == 0:

    if save_data_volt:
        #############################################
        print("Unifying the membrane potential data")
        #############################################
        data_volt = {}
        all_volt  = dict.fromkeys(["pyr_ca3","bas_ca3","olm_ca3"])
        keys = [] 
        for key in record_mp.keys():
            if record_mp[key] == True:
                keys.append(key)
        if keys: 
            all_volt["pyr_ca3"] = dict.fromkeys(keys)
            for key in keys:        
                all_volt["pyr_ca3"][key] = unify_data(f"{key}_volt", cells=net.__dict__["pyr_ca3"].cell,nr=4)
            
            data_volt["pyr_ca3"] = process_volt_data(all_volt["pyr_ca3"])

            for i in range(number_of_argvs):
                value = inputs_argvs[i]
                data_volt["pyr_ca3"][column_labels[i]] = [value]*len(data_volt["pyr_ca3"])
            
            title = f"volt_pyr_ca3_{argvs}.lzma"
            file_management.save_lzma( data_volt["pyr_ca3"], title, parent_dir = save_folder)

            if record_mp["soma"] == True:
                all_volt["bas_ca3"] = dict.fromkeys(["soma"])
                all_volt["olm_ca3"] = dict.fromkeys(["soma"])
                all_volt["bas_ca3"]["soma"] = unify_data("soma_volt", cells=net.__dict__["bas_ca3"].cell,nr=4)
                all_volt["olm_ca3"]["soma"] = unify_data("soma_volt", cells=net.__dict__["olm_ca3"].cell,nr=4)

                data_volt["bas_ca3"] = process_volt_data(all_volt["bas_ca3"])
                data_volt["olm_ca3"] = process_volt_data(all_volt["olm_ca3"])

                for i in range(number_of_argvs): # adding the input values
                    value = inputs_argvs[i]
                    data_volt["bas_ca3"][column_labels[i]] = [value]*len(data_volt["bas_ca3"])
                    data_volt["olm_ca3"][column_labels[i]] = [value]*len(data_volt["olm_ca3"])

                title = f"volt_bas_ca3_{argvs}.lzma"
                file_management.save_lzma( data_volt["bas_ca3"], title, parent_dir = save_folder)
                title = f"volt_olm_ca3_{argvs}.lzma"
                file_management.save_lzma( data_volt["olm_ca3"], title, parent_dir = save_folder)

            if net.create_ca1_network:
                if keys: 
                    all_volt["pyr_ca1"] = dict.fromkeys(keys)
                    for key in keys:        
                        all_volt["pyr_ca1"][key] = unify_data(f"{key}_volt", cells=net.__dict__["pyr_ca3"].cell,nr=4)

                data_volt["pyr_ca1"] = process_volt_data(all_volt["pyr_ca1"])
                for i in range(number_of_argvs): # add the input values
                    value = inputs_argvs[i]
                    data_volt["pyr_ca1"][column_labels[i]] = [value]*len(data_volt["pyr_ca1"])

                title = f"volt_pyramidal_ca1_{argvs}.lzma"
                file_management.save_lzma( data_volt["pyr_ca1"], title, parent_dir = save_folder)

                if record_mp["soma"] == True:
                    all_volt["bas_ca1"] = dict.fromkeys(["soma"])
                    all_volt["olm_ca1"] = dict.fromkeys(["soma"])
                    all_volt["bas_ca1"]["soma"] = unify_data("soma_volt", cells=net.__dict__["bas_ca1"].cell,nr=4)
                    all_volt["olm_ca1"]["soma"] = unify_data("soma_volt", cells=net.__dict__["olm_ca1"].cell,nr=4)
                    all_volt["cck_ca1"]["soma"] = unify_data("soma_volt", cells=net.__dict__["cck_ca1"].cell,nr=4)
                    
                    data_volt["bas_ca1"] = process_volt_data(all_volt["bas_ca1"])
                    data_volt["olm_ca1"] = process_volt_data(all_volt["olm_ca1"])
                    data_volt["cck_ca1"] = process_volt_data(all_volt["cck_ca1"])

                    for i in range(number_of_argvs): # adding the input values
                        value = inputs_argvs[i]
                        data_volt["bas_ca1"][column_labels[i]] = [value]*len(data_volt["bas_ca1"])
                        data_volt["olm_ca1"][column_labels[i]] = [value]*len(data_volt["olm_ca1"])
                        data_volt["cck_ca1"][column_labels[i]] = [value]*len(data_volt["cck_ca1"])

                    title = f"volt_bas_ca1_{argvs}.lzma"
                    file_management.save_lzma( data_volt["bas_ca1"], title, parent_dir = save_folder)
                    title = f"volt_olm_ca1_{argvs}.lzma"
                    file_management.save_lzma( data_volt["olm_ca1"], title, parent_dir = save_folder)
                    title = f"volt_cck_ca1_{argvs}.lzma"
                    file_management.save_lzma( data_volt["cck_ca1"], title, parent_dir = save_folder)

    if save_data_spikes:
        ################################
        print("Unifying the spike data")
        ################################
        all_spikes = dict.fromkeys(["pyr_ca3","bas_ca3","olm_ca3"])
        data_spikes = {} 

        all_spikes["pyr_ca3"] = unify_data("spike_times", net.__dict__["pyr_ca3"].cell)
        all_spikes["bas_ca3"] = unify_data("spike_times", net.__dict__["bas_ca3"].cell)
        all_spikes["olm_ca3"] = unify_data("spike_times", net.__dict__["olm_ca3"].cell)

        data_spikes["pyr_ca3"] = process_spike_data( all_spikes["pyr_ca3"] )
        data_spikes["bas_ca3"] = process_spike_data( all_spikes["bas_ca3"] )
        data_spikes["olm_ca3"] = process_spike_data( all_spikes["olm_ca3"] )

        for i in range(number_of_argvs): # adding the input values
            value = inputs_argvs[i]
            data_spikes["pyr_ca3"][column_labels[i]] = [value]*len(data_spikes["pyr_ca3"])
            data_spikes["bas_ca3"][column_labels[i]] = [value]*len(data_spikes["bas_ca3"])
            data_spikes["olm_ca3"][column_labels[i]] = [value]*len(data_spikes["olm_ca3"])
        
        title = f"external_inputs_pyr_ca3_{argvs}.lzma"
        file_management.save_lzma( data_spikes["pyr_ca3"], title, parent_dir = save_folder)
        title = f"spikes_bas_ca3_{argvs}.lzma"
        file_management.save_lzma( data_spikes["bas_ca3"], title, parent_dir = save_folder)
        title = f"spikes_olm_ca3_{argvs}.lzma"
        file_management.save_lzma( data_spikes["olm_ca3"], title, parent_dir = save_folder)

        if net.create_ca1_network: 
            all_spikes["pyr_ca1"] = unify_data("spike_times", net.__dict__["pyr_ca1"].cell)
            all_spikes["bas_ca1"] = unify_data("spike_times", net.__dict__["bas_ca1"].cell)
            all_spikes["olm_ca1"] = unify_data("spike_times", net.__dict__["olm_ca1"].cell)
            all_spikes["cck_ca1"] = unify_data("spike_times", net.__dict__["cck_ca1"].cell)

            data_spikes["pyr_ca1"] = process_spike_data( all_spikes["pyr_ca1"] )
            data_spikes["bas_ca1"] = process_spike_data( all_spikes["bas_ca1"] )
            data_spikes["olm_ca1"] = process_spike_data( all_spikes["olm_ca1"] )
            data_spikes["cck_ca1"] = process_spike_data( all_spikes["cck_ca1"] )

            for i in range(number_of_argvs): # adding the input values
                value = inputs_argvs[i]
                data_spikes["pyr_ca1"][column_labels[i]] = [value]*len(data_spikes["pyr_ca1"])
                data_spikes["bas_ca1"][column_labels[i]] = [value]*len(data_spikes["bas_ca1"])
                data_spikes["olm_ca1"][column_labels[i]] = [value]*len(data_spikes["olm_ca1"])
                data_spikes["cck_ca1"][column_labels[i]] = [value]*len(data_spikes["cck_ca1"])

            title = f"spikes_pyr_ca1_{argvs}.lzma"
            file_management.save_lzma( data_spikes["pyr_ca1"], title, parent_dir = save_folder)
            title = f"spikes_bas_ca1_{argvs}.lzma"
            file_management.save_lzma( data_spikes["bas_ca1"], title, parent_dir = save_folder)
            title = f"spikes_olm_ca1_{argvs}.lzma"
            file_management.save_lzma( data_spikes["olm_ca1"], title, parent_dir = save_folder)
            title = f"spikes_cck_ca1_{argvs}.lzma"

    if save_data_syn:
        if record_all_synapses:
            ################################
            print("Unifying the synaptic currents")
            ################################

            all_syn = {}
            all_syn["ca3"] = {}
            for key in synlist:
                all_syn["ca3"][key] = unify_data(f"i{key}", cells=net.__dict__["pyr_ca3"].cell,nr=4)
            data_syn = {} 
            data_syn["ca3"] = {}
            data_syn["ca3"] = process_syn_data( all_syn["ca3"] )

            for i in range(number_of_argvs): # adding the inpOut values
                value = inputs_argvs[i]
                data_syn["ca3"][column_labels[i]] = [value]*len(data_syn["ca3"])
            
            title = f"synapses_pyr_ca3_{argvs}.lzma"
            file_management.save_lzma( data_syn["ca3"], title, parent_dir = save_folder)

            if net.create_ca1_network: 
                all_syn["ca1"] = {}
                for key in synlist:
                    all_syn["ca1"][key] = unify_data(f"i{key}", cells=net.__dict__["pyr_ca1"].cell,nr=4)
                data_syn["ca1"] = {}
                data_syn["ca1"] = process_syn_data( all_syn["ca1"] )

                for i in range(len(sys.argv)+1): # adding the input values
                    value = inputs_argvs[i]
                    data_syn["ca1"][column_labels[i]] = [value]*len(data_syn["ca1"])
                
                title = f"synapses_pyr_ca1_{argvs}.lzma"
                file_management.save_lzma( data_syn["ca1"], title, parent_dir = save_folder )

    if save_data_lfp:

        # LFPsimpy is already implemented with the parallelization framework
        if record_lfp:
            data_lfp = {}
            data_lfp["ca3"] =  {"lfp":[],"electrode":[]}
            
            for i, lfp_ in enumerate(electrode["ca3"]):
                data_lfp["ca3"]["lfp"].extend( np.array(lfp_.values) )
                data_lfp["ca3"]["electrode"].extend( [zcoords[i]]*len(lfp_.values) ) # z-component inseatd of names, in case we put more
            
            data_lfp["ca3"] = pd.DataFrame(data_lfp["ca3"])
            for i in range(number_of_argvs):
                value = inputs_argvs[i]
                # data_lfp["ca3"][f"input{i+1}"] = []
                data_lfp["ca3"][column_labels[i]] = [value]*len(data_lfp["ca3"]) 
            
            title = f"lfp_ca3_{argvs}.lzma"
            file_management.save_lzma(data_lfp["ca3"], title, parent_dir=save_folder)

            if net.create_ca1_network:
                data_lfp["ca1"] = {"lfp":[],"electrode":[]}
                for i, lfp_ in enumerate(electrode["ca1"]):
                    data_lfp["ca1"]["lfp"].extend( np.array(lfp_.values) )
                    data_lfp["ca1"]["electrode"].extend( [zcoords[i]]*len(lfp_.values) )
                
                data_lfp["ca1"] = pd.DataFrame(data_lfp["ca1"])
                for i in range(number_of_argvs):
                    value = inputs_argvs[i]
                    # data_lfp["ca1"][f"input{i+1}"] = []
                    data_lfp["ca1"][column_labels[i]] = [value]*len(data_lfp["ca1"]) 
                
                title = f"lfp_ca1_{argvs}.lzma"
                file_management.save_lzma(data_lfp["ca1"], title, parent_dir=save_folder)
    
    if save_data_ica:
        if record_lfp:
            # LFPsimpy is already implemented with the parallelization framework
            data_ica = {}
            data_ica["ca3"] = {"ica":[],"electrode":[],"component":[]}
            
            for i, lfp_ in enumerate(electrode["ca3"]):
                for j,sublfp_ in enumerate(np.array(lfp_.values_per_section).T):
                    data_ica["ca3"]["ica"].extend( np.array(sublfp_) )
                    data_ica["ca3"]["electrode"].extend( [zcoords[i]]*len(sublfp_) )
                    data_ica["ca3"]["component"].extend( [zcoords[j]]*len(sublfp_) )

            data_ica["ca3"] = pd.DataFrame(data_ica["ca3"])
            for i in range(number_of_argvs):
                value = inputs_argvs[i]
                data_ica["ca3"][column_labels[i]] = [value]*len(data_ica["ca3"]) 
            
            title = f"ica_ca3_{argvs}.lzma"
            file_management.save_lzma(data_ica["ca3"], title, parent_dir=save_folder)

            if net.create_ca1_network:
                data_lfp["ca1"] = {"ica":[],"electrode":[],"component":[]}
                for i, lfp_ in enumerate(electrode["ca1"]):
                    for j,sublfp_ in enumerate(np.array(lfp_.values_per_section).T):
                        data_ica["ca1"]["ica"].extend( np.array(sublfp_) )
                        data_ica["ca1"]["electrode"].extend( [zcoords[i]]*len(sublfp_) )
                        data_ica["ca1"]["component"].extend( [["Bdend","soma","Adend1","Adend2","Adend3"][j]]*len(sublfp_) )

                data_ica["ca1"] = pd.DataFDrame(data_ica["ca1"])
                for i in range(number_of_argvs):
                    value = inputs_argvs[i]
                    data_ica["ca1"][column_labels[i]] =  [value]*len(data_ica["ca1"]) 

                title = f"ica_ca1_{argvs}.lzma"
                file_management.save_lzma(data_ica["ca1"], title, parent_dir=save_folder)
