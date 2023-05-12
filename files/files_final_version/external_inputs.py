"""
@author: Jaime

Generation of the ECII, ECII and DG inputs

# input generation with python

# date: 08/05/2023  
"""

import numpy as np
import pandas as pd
from inputs_parameters import *
import sys
import os
#sys.path.append("/home/jaime/Desktop/hippocampus/files")
import file_management

input1 = int(sys.argv[1]) # select the input category: baseline, ec2_faster_1, etc
input2 = int(sys.argv[2]) # select the label (population): sep, ec2, dg, etc

input_type = inputs[input1]
label_type = labels[input2]

# functions
def inputs_generation(sead, time_simulation, time_initial, sigma, period, ncells):
    idv, spikes = [],[]
    for i in range(ncells):
        np.random.seed(sead+i)# para que el tiempo no afecte # 16/06/2022
        tmean = time_initial
        tm, t0 = -1.0, 0.0
        while tmean <= time_simulation:
            while tm <= t0:
                tm = np.round( np.random.normal(loc=tmean, scale=sigma), 1 )
            spikes.append(tm)
            idv.append(i)
            tmean += period
            t0 = spikes[-1]
    idv = np.array(idv).astype(int)
    spikes = np.array(spikes)
    sead += ncells
    return sead, idv, spikes

def dg_burst_inputs(sead,time_simulation, time_initial, sigma, period, ncells, pburst=0.5, nspikes=6):
    std = sigma # sigma of gaussian
    std_ = 1.0   # decay of the exponential distributions for isi
    isimin = 3.5 # min isi according to literature
    idv, spikes = [],[]

    #tmean = time_initial+0.125*theta_rythm
    #time_ref = 1000.0 # tiempo a partir del cual cambio el pburst de 0.5 al que sea s
    pburst = [0.5,pburst]
    k=1
    for i in range(ncells):
        np.random.seed(sead+i)# para que el tiempo no afecte # 16/06/2022
        #tmean = time_initial+0.125*theta_rythm
        tmean = time_initial
        t0 = 0.0
        while tmean <= time_simulation:
            tm = -1.0
            while tm <= t0:
                tm = np.round(np.random.normal(loc=tmean, scale=std),1)
            x0 = np.random.rand(1)
            if x0<=pburst[k]:
                burst_sequence=nspikes
            else:
                burst_sequence=nspikes-1
            check_isi = np.array([1])
            kk=0
            while len(check_isi) > 0:
                isi = np.round(isimin+np.sort(np.random.exponential(scale=std_,size = burst_sequence-1)),1)
                check_isi = np.where(np.diff(isi)==0)[0]
            isi = np.append(0,isi)
            spikes.append(np.round(tm+isi,1))
            idv.append(np.ones(burst_sequence)*i)
            tmean += period
            t0 = spikes[-1][-1] # for avoding overlapping with the next spike of the next cycle

    idv = np.concatenate(idv).astype(int)
    spikes = np.concatenate(spikes)
    sead += ncells
    return sead, idv, spikes

# inputs path (do not forget to modify)
inputs_path = os.path.join( "/home/jaime/Desktop/hippocampus/external_inputs/", input_type)

ncells_per_label = [ inputs_params[input_type][lb]["ncells"] for lb in labels ]
ncells_total = np.sum(ncells_per_label)
ncells_per_label_cumsum = np.cumsum(ncells_per_label) 
ncells_per_label_cumsum = np.insert(ncells_per_label_cumsum,0,0)
ncells_per_label_cumsum = ncells_per_label_cumsum[:-1]

j = labels.index(label_type)

time_simulation = 30000.0 
columns = dict.fromkeys(["iseed","idvec","tvec"])
for key in columns.keys():
    columns[key] = []

ntrials = 20
for i in range(ntrials):
    sead = 1234 + i*ncells_total + ncells_per_label_cumsum[j]

    sigma      = inputs_params[input_type][label_type]["sigma"]
    ncells     = inputs_params[input_type][label_type]["ncells"]
    time_start = inputs_params[input_type][label_type]["time_start"]
    period     = inputs_params[input_type][label_type]["period"]

    if label_type == "dg_burst":
        sead,b,c = dg_burst_inputs(sead=sead,sigma=sigma,time_simulation=time_simulation,period=period,time_initial=time_start,ncells=ncells)
    else:
        sead,b,c = inputs_generation(sead=sead,sigma=sigma,time_simulation=time_simulation,period=period,time_initial=time_start,ncells=ncells)

    columns["iseed"].append(np.ones(len(b))*i)
    columns["idvec"].append(b)
    columns["tvec"].append(c)

columns["iseed"]  = np.concatenate(columns["iseed"]).astype(int)
columns["idvec"]  = np.concatenate(columns["idvec"])
columns["tvec"]   = np.concatenate(columns["tvec"])

df = pd.DataFrame(columns)
file_management.save_lzma(df, f"external_inputs_{label_type}", parent_dir=inputs_path)