'''
Date 09/02/2023

Modifico la entrada de los inputs de ECIII-360. 
Con la idea de que no haya que añadir ningun "fake delay" 

De primeras voy a utilizar el metodo mas sucio para ir lo mas rapido posible

Este fichero se trata de una copia del fichero folder54 de la carpeta de enero (nuredduna), 
donde escaneaba el efecto de la inhibición de las cck sobre las pyramidales. En este se observaba 
cómo la intensidad de dicha inhibición retrasaba la fase de de las pyramidales. 
Como simple referencia escojo un valor intermedio de dicha inhibición. 

De momento no voya añadir la sigma que había considerado en la proyección CA3-CA1
'''

import warnings
warnings.filterwarnings("ignore")
from neuron import h
import numpy as np
import sys
import os
from neuron.units import ms, mV
import time as tm
#sys.path.append('/home/dimitrios/Neurons/CA1model/files/')
#from neurons import *
from network_hippocampus_ca1 import *
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
    input1 = inputs_argvs[2] # first variable to scan
    input2 = inputs_argvs[3] # second variable to scan

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
    column_labels.append("input1")
    column_labels.append("input2")

    # column_labels ...

elif len(inputs_argvs) == 2: 
    # in case no inputs are given
    iseed  = inputs_argvs[0] # external inputs basal seed
    # control the seed for the external inputs generation (online generation, not implemented)
    # and select the trial of the external input (offline generation, implemented)
    ibseed = inputs_argvs[1] # background noise basal seed
    # control the seed for the background noise generation (online generation)
        # and would also select the trial of the background noise (offline generation, not implemented)

    # default in this example
    input1 = 0
    input2 = 0

    number_of_argvs = len(sys.argv[1:])
    argvs=""
    for i in range(number_of_argvs-1):
        argvs += sys.argv[1:][i]+'_'
    argvs+=sys.argv[-1]

    column_labels.append("iseed")
    column_labels.append("ibseed")
    column_labels.append("input1")
    column_labels.append("input2")

else: 
    print("provide at least two inputs")
    sys.exit()

h.dt = 0.1 # time resolution of the simulation
h.celsius = 34.0 # temperature of the simulation
simulation_time = 300.0 # total time of the simulation
time_resolution = 2.0 # time resolution of the measurements

DoMakeNoise = True # Control the external noise
MakeNetStim = True # Online generation of the background noise 
DoMakeExternalInputs = True # Control the external inputs
MakeCellStim         = False # Online generation of the external inputs (False means loading the external inputs from a file)

record_all_synapses = False # Record all the synapses
record_lfp          = False # Record the LFP and transmembrane currents

# synapses to be recorded
synlist = [ "Adend3GABA_olm","Adend3AMPA_ec3360","Adend3NMDA_ec3360", "Adend3GABA_noise", "Adend3AMPA_noise",
            "Adend2GABA_cck",
            "Adend1AMPA_pyrCA3", "Adend1NMDA_pyrCA3",
            "somaGABA_bas", "somaGABA_cck", "somaAMPA_noise", "somaGABA_noise"]

# membrane potential to be recorded
record_mp = {"Bdend": False, "soma": False, "Adend1": False, "Adend2": False, "Adend3": False}

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
inputs_folder = '/home/jaime/Desktop/hippocampus/external_inputs/baseline/'
save_folder = os.path.join(current_folder, "test_data")
#save_folder   = '/home/jaime/Desktop/hippocampus/external_inputs/baseline/'
# file_name = __file__
# file_name = file_name.replace(current_folder,"") 
# file_name = file_name.split('_')[-1][:-3] 
# save_dir = os.path.join(current_folder, file_name)

def unify_data(variable, cells, nr=1):
    f''' unify data from all cells in the network: membrane voltage, synaptic currents, spikes, etc '''
    local_data = { cell.id: list( np.round( cell.__dict__[variable].to_python(),nr)) for cell in cells }
    all_data = pc.py_alltoall( [local_data] + [None] * (pc.nhost() - 1) )
    return all_data

'''#######################################################################
                  Parameters to optimize in CA1
            (All parameters written and explicited below)
            i) neurons-neurons
               a) weights
               b) nsyns (number of connections)
               c) type of synapsis
            ii) inputs-neurons
               a) weights
               b) nsyns (number of connections)
               c) type of synapsis
###########################################################################'''

'''###########################################################################
                    Neurons to neurons projections
###########################################################################'''
# weights
gain2 = np.linspace(0,2,21)
gain = np.linspace(0.0,1.5,21)# esta va con la conexion de abajo 
weights_neurons_neurons["pyr_ca1_to_bas_ca1"] = [[ gain[5]*2.341e-3, gain[5]*1.5*1.38e-3]]
weights_neurons_neurons["pyr_ca1_to_olm_ca1"] = [[ 0.969e-3, 0.7e-3]] # [0.36e-3+scale*0.7e-3, 0.7e-3]
weight_seq = np.linspace(0,2e-3,21)
weights_neurons_neurons["pyr_ca1_to_cck_ca1"] = [[ 0.0*weight_seq[input2],0]] #[[0.10*weight_seq[5], 0.1*weight_seq[5]]] #[4.05e-3], [ 4.05e-3 ]] #[[cck_inhibition[input1]],[cck_inhibition[input1]]]#[
gain = np.linspace(0.5,1,21)
weights_neurons_neurons["bas_ca1_to_bas_ca1"] = [[ gain[10]*4.05e-3 ]] # initial value [[4.5e-3]]
weight = np.linspace(0,4e-3,21)
weights_neurons_neurons["bas_ca1_to_pyr_ca1"] = [[ 0.5*0.5*weight[10] ]] # initial [[0.576e-3]]
weights_neurons_neurons["olm_ca1_to_pyr_ca1"] = [[ 57.6e-3 ]]  # [0.8*72.0e-3] #*0.9

cck_inhibition = 8*np.linspace(0,1,21)*4.05e-4
gain = np.linspace(0,4,21)
weights_neurons_neurons["cck_ca1_to_pyr_ca1"] =  [[ cck_inhibition[10]], [cck_inhibition[10]]] #[[0.0],[0.0]] #[4.05e-3], [ 4.05e-3 ]] ##

weight = np.linspace(0,1e-3,21)
weights_neurons_neurons["cck_ca1_to_cck_ca1"] = [[ weight[10] ]]
# nsyns (number of projections)
ncons = np.arange(20,101,4)
nsyns_neurons_neurons["pyr_ca1_to_bas_ca1"] = [90] # [100] before the reduction of ca3
nsyns_neurons_neurons["pyr_ca1_to_olm_ca1"] = [10]
nsyns_neurons_neurons["pyr_ca1_to_cck_ca1"] = [90] # to be optimized
nsyns_neurons_neurons["bas_ca1_to_bas_ca1"] = [30] # [30] [60] before the reduction of ca3
nsyns_neurons_neurons["bas_ca1_to_pyr_ca1"] = [20] # [42] before the reduction of ca3
nsyns_neurons_neurons["olm_ca1_to_pyr_ca1"] = [10]
nsyns_neurons_neurons["cck_ca1_to_pyr_ca1"] = [30, 10]  #to optimizate keeping 70-30 of the total projections
nsyns_neurons_neurons["cck_ca1_to_cck_ca1"] = [30]

# syn (synapses type)
syn_neurons_neurons["pyr_ca1_to_bas_ca1"] = [["somaAMPA_pyr", "somaNMDA_pyr"]]
syn_neurons_neurons["pyr_ca1_to_olm_ca1"] = [["somaAMPA_pyr", "somaNMDA_pyr"]]
syn_neurons_neurons["pyr_ca1_to_cck_ca1"] = [["somaAMPA_pyr", "somaNMDA_pyr"]]
syn_neurons_neurons["bas_ca1_to_bas_ca1"] = [["somaGABA_bas"]]
syn_neurons_neurons["bas_ca1_to_pyr_ca1"] = [["somaGABA_bas"]]
syn_neurons_neurons["olm_ca1_to_pyr_ca1"] = [["Adend3GABA_olm"]] # previously  [["Adend3GABAf"]]
syn_neurons_neurons["cck_ca1_to_pyr_ca1"] = [["Adend2GABA_cck"], ["somaGABA_cck"]]
syn_neurons_neurons["cck_ca1_to_cck_ca1"] = [["somaGABA_cck"]]

# delay (synpatic delay)
delay_neurons_neurons["cck_ca1_to_pyr_ca1"] = [2.0,1.0]
delay_neurons_neurons["pyr_ca1_to_cck_ca1"] = [2.0,1.0]
delay_neurons_neurons["pyr_ca1_to_olm_ca1"] = [2.0,1.0]
delay_neurons_neurons["olm_ca1_to_pyr_ca1"] = [2.0,1.0]
delay_neurons_neurons["pyr_ca1_to_bas_ca1"] = [2.0,1.0]
delay_neurons_neurons["bas_ca1_to_pyr_ca1"] = [2.0,1.0]
delay_neurons_neurons["bas_ca1_to_bas_ca1"] = [2.0,1.0]
delay_neurons_neurons["cck_ca1_to_cck_ca1"] = [2.0,1.0]

'''###########################################################################
                    External inputs to neurons projections
###########################################################################'''
# weights
gain = np.linspace(1,10,21)
weights_inputs_neurons["sep_180_to_bas_ca1"] = [[ 6.4e-4 ]]
weights_inputs_neurons["sep_360_to_olm_ca1"] = [[ 3.2e-4 ]]
weights_inputs_neurons["sep_180_to_cck_ca1"] = [[ 3.2e-4 ]]
# those below null (by now)
weights_inputs_neurons["ec3_180_to_pyr_ca1"] = [[ 0.0, 0.0 ]]
#weights_inputs_neurons["ec3_180_to_bas_ca1"] = [[0.0,0.0]]
#weights_inputs_neurons["ec3_180_to_cck_ca1"] = [[0.0,0.0]]
gain = np.linspace(0,1,21)
weights_inputs_neurons["ec3_360_to_cck_ca1"] = [[ 4*5*3.3e-4*gain[12], 1.8e-4*gain[12]]] #[[3.3e-4, ]]
gain = np.linspace(0,1,21)
weights_inputs_neurons["ec3_360_to_bas_ca1"] = [[ 0.5*3.3e-4, 0.5*1.8e-4]]
weights_inputs_neurons["ec3_360_to_pyr_ca1"] = [[ 3.3e-4, 1.8e-4]] # [[2*3.3e-4, 2*1.8e-4]]#[[2*3.3e-4, 2*1.8e-4]]
weight = np.linspace(0,2,21)
weights_inputs_neurons["ec3_180_to_pyr_ca1"] = [[ 0.0*weight[8]*3.3e-4, 0.0*weight[8]*1.8e-4 ]] # [[2*3.3e-4, 2*1.8e-4]]#[[2*3.3e-4, 2*1.8e-4]]

# Schaffer colateral
weight_seq = np.linspace(0.75,1.25,21)*0.001045
gseq_pyr = np.linspace(0.3e-3, 0.7e-3, 21)
gseq_cck = np.linspace(1e-4,1e-3,21) # scale factor
delay_seq = np.linspace(2,8,21)

weight_seq = np.linspace(0,0.23e-3,21)
gain = np.linspace(0.25,0.5,21)
weights_inputs_neurons["pyr_ca3_to_pyr_ca1"] =  [[ gain[10]*0.75*0.32e-3, 4*0.25*0.5*0.23e-3 ]]#[[ 0.25*0.32e-3, 0.25*0.5*0.23e-3]]#weight_seq[input1] ]] # [[ gain[input1]*1.44e-4, gain[input2]*0.6e-4 ]] #[[1.5*weight_seq[5], 1.5*weight_seq[5]/2.0 ]] #[[ 1.5*weight_seq[4], 1.5*weight_seq[4]/2.0 ]]  # to be optimize

weight_seq = np.linspace(0,0.4e-3,21)
gain = np.linspace(0,2,21)
weights_inputs_neurons["pyr_ca3_to_bas_ca1"] = [ [0.8*0.5*3.3e-4, 4*0.8*0.5*1.8e-4] ] #[[ 0.9e-4,  0.5e-4 ]]
gain = np.linspace(0,2,21)
weights_inputs_neurons["pyr_ca3_to_cck_ca1"] = [ [ 0.25*0.16e-3, 0.25*0.014e-3 ] ] # to be optimize

ncon = np.linspace(25,125,21).astype(int)
nsyns_inputs_neurons["pyr_ca3_to_pyr_ca1"] = [ 100 ]  # theoretically 160 to be optimized
nsyns_inputs_neurons["pyr_ca3_to_bas_ca1"] = [ 100 ]  # theoretically 160 to be optimized
nsyns_inputs_neurons["pyr_ca3_to_cck_ca1"] = [ 100 ]  # theoretically 160, to be optimized

delay_seq = np.arange(2,23,1)
delay_external = 10.0
delay_CA3 = 5.0
delay_CA3_seq = np.linspace(2,12,21)
sigma_seq = np.linspace(0,20,21)
delay_inputs_neurons["pyr_ca3_to_pyr_ca1"] = [ delay_CA3_seq[input1], 10.0 ]# [delay_seq[input2]]
delay_inputs_neurons["pyr_ca3_to_bas_ca1"] = [ delay_CA3_seq[input1], 1.0 ]
delay_inputs_neurons["pyr_ca3_to_cck_ca1"] = [ delay_CA3_seq[input1], 1.0 ]

delay_inputs_neurons["sep_360_to_olm_ca1"] = [ delay_external,1.0 ]
delay_inputs_neurons["sep_180_to_bas_ca1"] = [ delay_external,1.0 ]
delay_inputs_neurons["sep_180_to_cck_ca1"] = [ delay_external,1.0 ]

delay_inputs_neurons["ec3_360_to_cck_ca1"] = [ delay_external,1.0 ]
delay_inputs_neurons["ec3_360_to_pyr_ca1"] = [ delay_external,1.0 ]
delay_inputs_neurons["ec3_360_to_bas_ca1"] = [ delay_external,1.0 ]

# not used yet
delay_inputs_neurons["ec3_180_to_pyr_ca1"] = [ delay_external,1.0 ]
delay_inputs_neurons["ec3_180_to_cck_ca1"] = [ delay_external,1.0 ]
delay_inputs_neurons["ec3_180_to_bas_ca1"] = [ delay_external,1.0 ]

# nsyns (number of projections)
nsyns_inputs_neurons["sep_180_to_bas_ca1"] = [10]
nsyns_inputs_neurons["sep_360_to_olm_ca1"] = [10]
nsyns_inputs_neurons["sep_180_to_cck_ca1"] = [10]
nsyns_inputs_neurons["ec3_180_to_pyr_ca1"] = [25]
# nsyns_inputs_neurons["ec3_180_to_bas_ca1"] = [25]
# nsyns_inputs_neurons["ec3_180_to_cck_ca1"] = [25]
nsyns_inputs_neurons["ec3_360_to_pyr_ca1"] = [25]
nsyns_inputs_neurons["ec3_360_to_bas_ca1"] = [25]
nsyns_inputs_neurons["ec3_360_to_cck_ca1"] = [25]

# syn
syn_inputs_neurons["sep_180_to_bas_ca1"] = [["somaGABA_sep180"]]
syn_inputs_neurons["sep_360_to_olm_ca1"] = [["somaGABA_sep360"]]
syn_inputs_neurons["sep_180_to_cck_ca1"] = [["somaGABA_sep180"]]
syn_inputs_neurons["ec3_180_to_pyr_ca1"] = [["Adend3AMPA_ec3180","Adend3NMDA_ec3180"]]
# syn_inputs_neurons["ec3_180_to_bas_ca1"] = [["somaAMPAf", "somaNMDA"]]
# syn_inputs_neurons["ec3_180_to_cck_ca1"] = [["somaAMPAf", "somaNMDA"]]
syn_inputs_neurons["ec3_360_to_pyr_ca1"] = [["Adend3AMPA_ec3360", "Adend3NMDA_ec3360"]]
syn_inputs_neurons["ec3_360_to_bas_ca1"] = [["somaAMPA_ec3360", "somaNMDA_ec3360"]]
syn_inputs_neurons["ec3_360_to_cck_ca1"] = [["somaAMPA_ec3360", "somaNMDA_ec3360"]]
syn_inputs_neurons["pyr_ca3_to_pyr_ca1"] = [["Adend1AMPA_pyrCA3", "Adend1NMDA_pyrCA3"]]
syn_inputs_neurons["pyr_ca3_to_bas_ca1"] = [["somaAMPA_pyrCA3", "somaNMDA_pyrCA3"]]
syn_inputs_neurons["pyr_ca3_to_cck_ca1"] = [["somaAMPA_pyrCA3", "somaNMDA_pyrCA3"]]

'''###########################################################################
                    External background noise
###########################################################################'''
weights_noise_neurons["bas_ca1"]["somaAMPA_noise"]   = 0.0125e-3
gain = np.linspace(1,10,21)
weights_noise_neurons["pyr_ca1"]["somaAMPA_noise"]   = 5*0.0125e-3
weights_noise_neurons["pyr_ca1"]["Adend3AMPA_noise"] = 0.0125e-3

##############################################################################
# remove the inhibition from the non-basket interneurons
#weights_neurons_neurons["olm_ca1_to_pyr_ca1"] = [[0]]
#weights_neurons_neurons["cck_ca1_to_pyr_ca1"] = [[0], [0]] 
#weights_neurons_neurons["cck_ca1_to_cck_ca1"] = [[0]]

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

'''###########################################################################
                           saving the currents
###########################################################################'''

synlist = [ "Adend3GABA_olm","Adend3AMPA_ec3360","Adend3NMDA_ec3360", "Adend3GABA_noise", "Adend3AMPA_noise",
             "Adend2GABA_cck",
             "Adend1AMPA_pyrCA3", "Adend1NMDA_pyrCA3",
             "somaGABA_bas", "somaGABA_cck", "somaAMPA_noise", "somaGABA_noise"]

'''############################################################################
                                The network
#############################################################################'''
record_mp = {"soma":False,"Bdend": False, "Adend1":False, "Adend2":False, "Adend3":False}
net = Network( weights_inputs_neurons  = weights_inputs_neurons,
               nsyns_inputs_neurons    = nsyns_inputs_neurons,
               syn_inputs_neurons      = syn_inputs_neurons,
               weights_neurons_neurons = weights_neurons_neurons,
               nsyns_neurons_neurons   = nsyns_neurons_neurons,
               syn_neurons_neurons     = syn_neurons_neurons,
               weights_noise_neurons   = weights_noise_neurons,
               delay_inputs_neurons    = delay_inputs_neurons,
               delay_neurons_neurons   = delay_neurons_neurons,
               iseed=3421*(2*iseed+1),  # semilla para los inputs
               bseed=27*(2*ibseed+1),   # semilla del ruido
               DoMakeNoise = DoMakeNoise,
               DoMakeExternalInputs = DoMakeExternalInputs,
               MakeCellStim = MakeCellStim,
               MakeNetStim  = MakeNetStim,
               inputs_folder = inputs_folder,
               external_inputs_iseed = iseed,
               external_inputs_ibseed = ibseed,
               n_pyr_ca3 = 800,
               n_bas_ca3 = 100,
               n_olm_ca3 = 30,
               n_pyr_ca1 = 800,
               n_bas_ca1 = 100,
               n_olm_ca1 = 30,
               n_cck_ca1 = 70,
               connections = True,
               nseg_pyr= 3, 
               record_mp = record_mp,
               resolution = time_resolution)

tau2seq = np.linspace(1,5,21)
tauNMDAseq = np.linspace(15,40,21)

for cell in net.cck_ca1.cell:
    cell.somaInj.amp = 32.5*1e-3 # fixed 15/11
for cell in net.pyr_ca1.cell:
    cell.Adend3GABA_olm.syn.tau2 = 20.0
# for cell in net.pyr_ca1.cell: # la putada es que aqui estoy cambiando tambien las sinapsis de las basket. 
#     # habria que usar el fichero de github
#     cell.Adend2GABA_cck.syn.tau2 = tau2seq[input1]
#     cell.somaGABA_cck.syn.tau2   = tau2seq[input1]
for cell in net.cck_ca1.cell: 
    cell.somaAMPA_ec3360.syn.tau2 = 2.5#tau2seq[input1]
for cell in net.pyr_ca1.cell: 
    cell.Adend1NMDA_pyrCA3.syn.tau1NMDA = 15#tauNMDAseq[input1]
    cell.Adend2GABA_cck.syn.tau2 = 4.0#tau2seq[input1]
    cell.somaGABA_cck.syn.tau2   = 4.0#tau2seq[input1]
    
# modification of the threshold of the cck 
for cell in net.cck_ca1.cell: 
    cell.spike_detector.threshold = -10.0

if net.DoMakeNoise:
    if net.MakeNetStim:
        net.set_noise_inputs(simulation_time) # set background noise and external input
        net.init_NetStims()  # init rngs of background

if net.DoMakeExternalInputs:
    if net.MakeCellStim:
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
        if key.endswith("ca1"):
            idv      = net.external_inputs_data[key]["idv"]
            spikes   = net.external_inputs_data[key]["spikes"]
            trg      = net.external_inputs_data[key]["population"]
            syn_list = net.external_inputs_data[key]["synapses"]
            delays   = net.external_inputs_data[key]["delays"]
            w_list   = net.external_inputs_data[key]["weights"]
            conn     = net.external_inputs_data[key]["connectivity"]
            if key.startswith("ec3_360"): # "quick solution to make cck spike before"
                spikes -= 10
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


if record_lfp:
    # electrode_y_coords = np.arange(0,850,50)
    zcoords = [-100,10,85,235, 385]
    # electrode_number = len(electrode_y_coords)
    electrode = {}
    electrode["ca1"] = [ ]

    for i,z in enumerate(zcoords):
        electrode["ca1"].append( LfpElectrode(x=25.0, y=25.0, z=z, sampling_period=1, neuron_type = "Pyramidal CA3"))

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
diff  = tock-tick
horas = int(diff/3600)
diff  = diff-horas*3600
mint  = int(diff/60)
diff  = diff-mint*60
seg   = int(diff)

print("Simulation done after", horas,'h', mint, 'min', seg, 's')
print("-------------------------------------------")
print("Preparing the data to be saved") 

# only saving the data if the simulation was run in the master node

if pc.id() == 0:

    if save_data_volt:
        #############################################
        print("Unifying the membrane potential data")
        #############################################
        data_volt = {}
        all_volt  = dict.fromkeys(["pyr_ca1","bas_ca1","olm_ca1","cck_ca1"])
        keys = [] 
        for key in record_mp.keys():
            if record_mp[key] == True:
                keys.append(key)

        if keys: 
            all_volt["pyr_ca1"] = dict.fromkeys(keys)
            for key in keys:        
                all_volt["pyr_ca1"][key] = unify_data(f"{key}_volt", cells=net.__dict__["pyr_ca1"].cell,nr=4)

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
            data_syn = {}
            all_syn = {}
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
            data_ica["ca1"] = {"ica":[],"electrode":[],"component":[]}
            for i, lfp_ in enumerate(electrode["ca1"]):
                for j,sublfp_ in enumerate(np.array(lfp_.values_per_section).T):
                    data_ica["ca1"]["ica"].extend( np.array(sublfp_) )
                    data_ica["ca1"]["electrode"].extend( [zcoords[i]]*len(sublfp_) )
                    data_ica["ca1"]["component"].extend( [zcoords[j]]*len(sublfp_) )

            data_ica["ca1"] = pd.DataFrame(data_ica["ca1"])
            for i in range(number_of_argvs):
                value = inputs_argvs[i]
                # data_ica["ca1"][f"input{i+1}"] = []
                data_ica["ca1"][column_labels[i]] =  [value]*len(data_ica["ca1"]) 

            title = f"ica_ca1_{argvs}.lzma"
            file_management.save_lzma(data_ica["ca1"], title, parent_dir=save_folder)
