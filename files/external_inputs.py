# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 12:43:43 2021

@author: Jaime

Generation of the ECII, ECII and DG inputs
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from neuron import h,gui
global current_folder
global data_folder
global data_file
global random_stream_offset_,idns

random_stream_offset_ =  1000
class RandomStream():
    def __init__(self, sead,generator):
        self.sead = sead
        self.generator = generator
        self.r = h.Random()

    def set_generator(self):
        if self.generator == 'gaussian':
            g = self.r.normal(0,1)
        if self.generator == 'exponential':
            g = self.r.negexp(1)

    def start(self):
        return self.r.MCellRan4(self.sead,self.sead)

    def repick(self):
        return self.r.repick()

def make_new_external_inputs(time_initial,sead):
    '''
    simdur: simulation time
    seed :  seed for the generation of random numbers
    relative_time: relative time between MS and ECII inputs
    time_initial: time of the first generation of a MS inputs. Altough
    MS inputs are not generated in this function, this parameters is used
    as a reference.
    '''
    theta_rythm = 125    # ms (8 Hz)

    ###############################################################################
    ''' medial septum '''
    ###############################################################################
    # MS 180 inhibition
    sep_180_num   = 1000           # number of SEP spikes
    sep_180_start = time_initial+theta_rythm/2   # time of first SEP 180 spike
    sep_180_int   = theta_rythm    # SEP spike ISI (during burst)
    sep_180_noise = 0.05		   # SEP ISI noise
    sep_180_bint  = theta_rythm/2  # SEP interburst interval
    sep_180_blen  = theta_rythm/2  # SEP burst length

    # MS 360 inhibition
    sep_360_num   = 1000		    # number of SEP spikes
    sep_360_start = time_initial	# time of first SEP spike
    sep_360_int   = theta_rythm 	# SEP spike ISI (during burst)
    sep_360_noise = 0.05			# SEP ISI noise
    sep_360_bint  = theta_rythm/2 	# SEP interburst interval
    sep_360_blen  = theta_rythm/2	# SEP burst length

    ############################################################################
    ''' ec2 '''
    ############################################################################
    # ec2-180 excitation
    ec2_180_num   = 1000			# number of ec2-180 spikes
    ec2_180_start = time_initial+theta_rythm/2    # time of first ec2-180 spike
    ec2_180_int   = theta_rythm	    # ec2-180 spike ISI
    ec2_180_noise = 0.15		    # ec2-180 ISI noise
    ec2_180_bint  = theta_rythm     # ec2-180 interburst interval
    ec2_180_blen  = theta_rythm/2   # ec2-180 burst length

    #// ec2-360 excitation
    ec2_360_num   = 1000			# number of ec2-360 spikes
    ec2_360_start = time_initial	#// time of first ec2-360 spike
    ec2_360_int   = theta_rythm 	# ec2-360 spike ISI
    ec2_360_noise = 0.15			# ec2-360 ISI noise
    ec2_360_bint = theta_rythm		# ec2-360 interburst interval
    ec2_360_blen = theta_rythm/2    # ec2-360 burst length

    ############################################################################
    ''' ec3  '''
    ############################################################################
    # ec3-180 excitation
    ec3_180_num   = 1000			# number of ec3-180 spikes
    ec3_180_start = time_initial+theta_rythm/2    # time of first ec3-180 spike
    ec3_180_int   = theta_rythm	    # ec3-180 spike ISI
    ec3_180_noise = 0.15		    # ec3-180 ISI noise
    ec3_180_bint  = theta_rythm     # ec3-180 interburst interval
    ec3_180_blen  = theta_rythm/2   # ec3-180 burst length

    #// ec3-360 excitation
    ec3_360_num   = 1000			# number of ec3-360 spikes
    ec3_360_start = time_initial	#// time of first ec3-360 spike
    ec3_360_int   = theta_rythm 	# ec3-360 spike ISI
    ec3_360_noise = 0.15			# ec3-360 ISI noise
    ec3_360_bint = theta_rythm		# ec3-360 interburst interval
    ec3_360_blen = theta_rythm/2    # ec3-360 burst length

    ###############################################################################
    ''' dg '''
    ###############################################################################
    # dg burst inputs >40 Hz (Buzsaki analysis 2009)
    dg_burst_num   = 1000			   # number of dg burst spikes
    dg_burst_start = time_initial+0.125*theta_rythm   # time of first dg burst spike
    dg_burst_int   = theta_rythm	   # dg burst spike ISI (during burst)
    dg_burst_noise = 0.12			   # dg burst ISI noise
    dg_burst_bint  = theta_rythm/2	   # dg burst interburst interval
    dg_burst_blen  = theta_rythm/2	   # dg burst burst length

    # dg regular inputs 8Hz (Buzsaki analysis 2009)
    dg_regular_num   = 1000			   # number of dg regular spikes
    dg_regular_start = time_initial+0.375*theta_rythm # time of first dg regular spike
    dg_regular_int   = theta_rythm	   # dg regular spike ISI (during burst)
    dg_regular_noise = 0.12			   # dg regular ISI noise

    sep_180_neurons = 10
    sep_360_neurons = 10
    ec2_180_neurons = 1000
    ec2_360_neurons = 1000
    ec3_180_neurons = 1000
    ec3_360_neurons = 1000
    dg_burst_neurons = 1000
    dg_regular_neurons = 1000

    nslist = [ [] for i in range(7)]
    nrlist = []
    nrseadlist = []
    nslabel = []
    ###########################################################################
    ''' medial septum '''
    ###########################################################################
    for i in range(sep_180_neurons):
        ns = h.BurstStim()

        ns.interval = sep_180_int
        ns.number   = sep_180_num
        ns.start    = sep_180_start
        ns.noise    = sep_180_noise
        ns.burstint = sep_180_bint
        ns.burstlen = sep_180_blen

        rs = RandomStream(sead,'gaussian')
        rs.set_generator()
        ns.noiseFromRandom(rs.r)
        rs.start()

        nslist[0].append(ns)
        nrlist.append(rs)
        nrseadlist.append(rs.sead)
        sead+=1

    nslabel.append(i+1)
    for i in range(sep_360_neurons):
        ns = h.BurstStim()

        ns.interval = sep_360_int
        ns.number   = sep_360_num
        ns.start    = sep_360_start
        ns.noise    = sep_360_noise
        ns.burstint = sep_360_bint
        ns.burstlen = sep_360_blen

        rs = RandomStream(sead,'gaussian')
        rs.set_generator()
        ns.noiseFromRandom(rs.r)
        rs.start()

        nslist[1].append(ns)
        nrlist.append(rs)
        nrseadlist.append(rs.sead)
        sead+=1
    nslabel.append(i+1)

    ###########################################################################
    ''' ec2 '''
    ###########################################################################
    for i in range(ec2_180_neurons):
        ns = h.RegnStim()

        ns.interval = ec2_180_int
        ns.number   = ec2_180_num
        ns.start    = ec2_180_start
        ns.noise    = ec2_180_noise

        rs = RandomStream(sead,'gaussian')
        rs.set_generator()
        ns.noiseFromRandom(rs.r)
        rs.start()

        nslist[2].append(ns)
        nrlist.append(rs)
        nrseadlist.append(rs.sead)
        sead+=1
    nslabel.append(i+1)

    for i in range(ec2_360_neurons):
        ns = h.RegnStim()

        ns.interval = ec2_360_int
        ns.number   = ec2_360_num
        ns.start    = ec2_360_start
        ns.noise    = ec2_360_noise

        rs = RandomStream(sead,'gaussian')
        rs.set_generator()
        ns.noiseFromRandom(rs.r)
        rs.start()

        nslist[3].append(ns)
        nrlist.append(rs)
        nrseadlist.append(rs.sead)
        sead+=1
    nslabel.append(i+1)

    ###########################################################################
    ''' ec3 '''
    ###########################################################################
    for i in range(ec3_180_neurons):
        ns = h.RegnStim()

        ns.interval = ec3_180_int
        ns.number   = ec3_180_num
        ns.start    = ec3_180_start
        ns.noise    = ec3_180_noise

        rs = RandomStream(sead,'gaussian')
        rs.set_generator()
        ns.noiseFromRandom(rs.r)
        rs.start()

        nslist[4].append(ns)
        nrlist.append(rs)
        nrseadlist.append(rs.sead)
        sead+=1
    nslabel.append(i+1)

    for i in range(ec3_360_neurons):
        ns = h.RegnStim()

        ns.interval = ec3_360_int
        ns.number   = ec3_360_num
        ns.start    = ec3_360_start
        ns.noise    = ec3_360_noise

        rs = RandomStream(sead,'gaussian')
        rs.set_generator()
        ns.noiseFromRandom(rs.r)
        rs.start()

        nslist[5].append(ns)
        nrlist.append(rs)
        nrseadlist.append(rs.sead)
        sead+=1
    nslabel.append(i+1)

    ###########################################################################
    ''' dg '''
    ###########################################################################
    #for i in range(dg_burst_neurons):
    #    ns = h.BurstStim()
    #
    #   ns.interval = dg_burst_int
    #   ns.number   = dg_burst_num
    #   ns.start    = dg_burst_start
    #   ns.noise    = dg_burst_noise
    #   ns.burstint = dg_burst_bint
    #   ns.burstlen = dg_burst_blen

    #   rs = RandomStream(sead,'gaussian')
    #   rs.set_generator()
    #   ns.noiseFromRandom(rs.r)
    #   rs.start()

    #    nslist.append(ns)
    #   nrlist.append(rs)
    #   nrseadlist.append(rs.sead)
    #   sead+=1
    #nslabel.append(i+1)

    for i in range(dg_regular_neurons):
        ns = h.RegnStim()
        ns.interval = dg_regular_int
        ns.number   = dg_regular_num
        ns.start    = dg_regular_start
        ns.noise    = dg_regular_noise

        rs = RandomStream(sead,'gaussian')
        rs.set_generator()
        ns.noiseFromRandom(rs.r)
        rs.start()

        nslist[6].append(ns)
        nrlist.append(rs)
        nrseadlist.append(rs.sead)
        sead+=1
    nslabel.append(i+1)

    #nslabel = np.cumsum(nslabel).astype(int)
    #nslabel = np.insert(nslabel, 0, 0, axis=0)

    labels = ['sep_180', 'sep_360', 'ec2_180', 'ec2_360', 'ec3_180','ec3_360','dg_regular', 'dg_burst']
    start_time = [sep_180_start, sep_360_start, ec2_180_start, ec2_360_start,
                  ec3_180_start, ec3_360_start, dg_regular_start, dg_burst_start]
    return sead, nslist, nrlist, nslabel, labels, start_time

##################################################################################
# actualizacion 15/07/2022 
##################################################################################
def inputs_generation(sead,time_simulation, time_initial, sigma, theta_rythm=125.0, ncells=1000):
    idv, spikes = [],[]
    for i in range(ncells):
        np.random.seed(sead+i)# para que el tiempo no afecte # 16/06/2022
        tmean = time_initial
        while tmean <= time_simulation:
            tm = np.round( np.random.normal(loc=tmean, scale=sigma), 1 )
            spikes.append(tm)
            idv.append(i)
            tmean += theta_rythm
    idv = np.array(idv).astype(int)
    spikes = np.array(spikes)

    return sead, idv, spikes

def dg_burst_inputs(sead,time_simulation, time_initial, sigma, pburst=0.5,theta_rythm=125.0, nspikes=6,ncells=1000):
    std = sigma # sigma of gaussian
    std_ = 1.0 # decay of the exponential distributions for isi
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
            while tm<=t0:
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
            tmean += theta_rythm
            t0 = spikes[-1][-1] # for avoding overlapping with the next spike of the next cycle

    idv = np.concatenate(idv).astype(int)
    spikes = np.concatenate(spikes)
    return sead, idv, spikes

def make_external_inputs_python(time_simulation, sead, pburst=0.5):
    
    labels = ['sep_180', 'sep_360', 'ec2_180', 'ec2_360', 'ec3_180','ec3_360','dg_regular', 'dg_burst']
    idvec  = dict.fromkeys(labels)
    tvec   = dict.fromkeys(labels)
    ncells = dict.fromkeys(labels)

    # constants 
    time_initial = 50
    theta_rythm = 125 
    time_start = dict.fromkeys(labels)
    ncells = dict.fromkeys(labels)
    sigma  = dict.fromkeys(labels)

    time_start["sep_180"] = time_initial+theta_rythm/2.0
    time_start["sep_360"] = time_initial
    time_start["ec2_180"] = time_initial+theta_rythm/2.0
    time_start["ec2_360"] = time_initial
    time_start["ec3_180"] = time_initial+theta_rythm/2.0
    time_start["ec3_360"] = time_initial
    time_start["dg_regular"] = time_initial+0.375*theta_rythm
    time_start["dg_burst"] = time_initial+0.125*theta_rythm
    
    ncells["sep_180"] = 10
    ncells["sep_360"] = 10
    ncells["ec2_180"] = 1000
    ncells["ec2_360"] = 1000
    ncells["ec3_180"] = 1000
    ncells["ec3_360"] = 1000
    ncells["dg_regular"] = 1000
    ncells["dg_burst"] = 1000 

    sigma["sep_180"] = 0.2*theta_rythm
    sigma["sep_360"] = 0.2*theta_rythm
    sigma["ec2_180"] = 0.2*theta_rythm
    sigma["ec2_360"] = 0.2*theta_rythm
    sigma["ec3_180"] = 0.2*theta_rythm
    sigma["ec3_360"] = 0.2*theta_rythm
    sigma["dg_regular"] = 0.2*theta_rythm
    sigma["dg_burst"] = 0.2*theta_rythm

    for lb in labels[:-1]: # except dg burst 
        sead, idvec[lb], tvec[lb] = inputs_generation(sead=sead,sigma=sigma[lb],time_simulation=time_simulation,
                                                    time_initial=time_start[lb],ncells=ncells[lb])
    lb = "dg_burst"
    sead, idvec[lb], tvec[lb] = dg_burst_inputs(sead=sead, time_simulation=time_simulation, sigma=sigma[lb],
                                                time_initial=time_start[lb], pburst=pburst, ncells=ncells[lb])
    
    start_time = [ time_start[lb] for lb in labels ]
    nslabel    = [ ncells[lb] for lb in labels ]
    
    return idvec, tvec, nslabel, labels, start_time

'''####################################################################################'''
# anadido 20/04/2022
def dg_burst_inputs_pre(sead,pburst,time_simulation,time_initial=50.0,theta_rythm=125.0,
                    nspikes=6,ncells=1000):

    std = 0.2*theta_rythm # sigma of gaussian
    std_ = 1.0 # decay of the exponential distributions for isi
    isimin = 3.5 # min isi according to literature
    spk,spikes = [],[]
    idv_, idv  = [],[]
    tmean = time_initial+0.125*theta_rythm
    time_ref = 1000.0 # tiempo a partir del cual cambio el pburst de 0.5 al que sea s
    pburst = [0.5,pburst]
    k=1
    for i in range(ncells):
        np.random.seed(sead+i)# para que el tiempo no afecte # 16/06/2022
        tmean = time_initial+0.125*theta_rythm
        t0 = 0.0
        while tmean <= time_simulation:
            tm = -1.0
            while tm<=t0:
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
            tmean += theta_rythm
            t0 = spikes[-1][-1] # for avoding overlapping with the next spike of the next cycle

    idv = np.concatenate(idv).astype(int)
    spikes = np.concatenate(spikes)
    return idv, spikes, ncells

# anadido 23/09/2021
def dg_burst_inputs_old(sead,pburst,time_simulation,time_initial=50.0,theta_rythm=125.0,
                    nspikes=6,ncells=1000):

    std = 0.5*0.2*theta_rythm # sigma of gaussian # ojo
    std_ = 1.0 # decay of the exponential distributions for isi
    isimin = 3.5 # min isi according to literature

    np.random.seed(sead)
    spk,spikes = [],[]
    idv_, idv  = [],[]
    tmean = time_initial+0.125*theta_rythm
    #burst_sequence = 5
    time_ref = 1000.0 # tiempo a partir del cual cambio el pburst de 0.5 al que sea s
    pburst = [0.5,pburst]
    for k,time_ in enumerate([time_ref, time_simulation]):
        while tmean <= time_:
            spk.append(np.random.normal(loc=tmean,scale=std, size=ncells))
            idv_.append(np.arange(ncells))
            for i,tm in zip(idv_[-1],spk[-1]):
                x0 = np.random.rand(1)
                if x0<=pburst[k]:
                    burst_sequence=nspikes
                elif x0>pburst[k]:
                    burst_sequence=nspikes-1
                isi = np.round(isimin+np.sort(np.random.exponential(scale=std_,size = burst_sequence-1)),1)
                kk=0
                while len( np.where(np.diff(isi)==0)[0] )>0:
                    isi = np.round(isimin+np.sort(np.random.exponential(scale=std_,size = burst_sequence-1)),1)
                    isi = np.round(isimin+np.sort(np.random.exponential(scale=std_,size = burst_sequence-1)),1)
                    #esto ultimo lo hago, por que ha dado la casualidad de que el input del intervalo anterior ha conincidio
                    # con el input del siguiente ciclo... para solo un ciclo y para una sola neurona artificial.. asi que lo
                    # mejor seria ir neurona por neurona para evitar todo este jaleo. Pero esto lo haremos mas adelante.
                    kk+=1
                isi = np.append(0,isi)
                spikes.append(np.sort(np.unique(np.round(tm+isi,1))))
                idv.append(np.ones(burst_sequence)*i)
                #np.random.normal(loc=t,scale=std_,size = 5)
                #if i==817:
                #    print(np.round(tm+isi,1))
            tmean += theta_rythm

    spk = np.concatenate(spk)
    idv_ = np.concatenate(idv_)

    spikes = np.concatenate(spikes)
    idv    = np.concatenate(idv)

    idv = idv[np.argsort(spikes)]
    spikes = np.sort(spikes)
    print(np.min(idv), np.max(idv))
    return idv, spikes, ncells

################################################################################
def make_external_inputs_ca3(simdur, seed,relative_time,time_initial):
    '''
    simdur: simulation time
    seed :  seed for the generation of random numbers
    relative_time: relative time between MS and ECII inputs
    time_initial: time of the first generation of a MS inputs. Altough
    MS inputs are not generated in this function, this parameters is used
    as a reference.
    '''
    theta_rythm = 125 #ms
    ISI = 2*theta_rythm
    time_limit = simdur
    time_offset = relative_time

    # data extracted form Mizuseki et al 2009 fig 8. Results from gaussian fits
    ec2_4_diff, ec2_4_sg =    0.0,  9.54 #ms
    ec3_4_diff, ec3_4_sg = -60.04, 13.20 #ms
    dg_4_diff,  dg_4_sg  =  12.19, 16.51 #ms
    dg_40_diff, dg_40_sg =  39.38, 14.69 #ms

    np.random.seed(seed)
    ###########################################################################
    ''' EC2 inputs '''
    ###########################################################################
    time_input = time_offset + time_initial
    time_sg = ec2_4_sg

    EC2_inputs_mean   = []
    EC2_input_label_a = []
    EC2_input_label_b = []
    EC2_input_a, EC2_input_b = [],[]
    while time_input < time_limit:
        EC2_input_a.append(np.random.normal(time_input,time_sg,size=100))
        EC2_input_b.append(np.random.normal(time_input+theta_rythm,time_sg,size=100))
        EC2_input_label_a+= list(range(100))
        EC2_input_label_b+= list(range(100,200))
        EC2_inputs_mean.append(time_input)
        EC2_inputs_mean.append(time_input+theta_rythm)
        time_input+=ISI

    EC2_input_a = np.concatenate(EC2_input_a)
    EC2_input_b = np.concatenate(EC2_input_b)
    EC2_input_a = np.round(EC2_input_a,2)
    EC2_input_b = np.round(EC2_input_b,2)

    ###########################################################################
    ''' DG inputs '''
    ###########################################################################
    time_input_4  = time_offset+time_initial+theta_rythm-dg_4_diff
    time_input_40 = time_offset+time_initial+theta_rythm-dg_40_diff
    time_sg_4     = dg_4_sg
    time_sg_40    = dg_40_sg

    DG_input_label   = np.arange(200)
    DG_input_label_a = DG_input_label[:83]
    DG_input_label_b = DG_input_label[83:166]
    DG_input_label_c = DG_input_label[166:]

    DG_input_label_a,  DG_input_a = [], []
    DG_input_label_b,  DG_input_b = [], []
    DG_input_label_c,  DG_input_c = [], []

    time_input = time_input_4
    while time_input < time_limit: #4Hz inputs
        DG_input_a.append(np.random.normal(time_input,time_sg_4,size=83))
        DG_input_b.append(np.random.normal(time_input+theta_rythm,time_sg_4,size=83))
        DG_input_label_a += list(range(83))
        DG_input_label_b += list(range(83,166))
        time_input+=ISI

    aux = np.random.rand(34)
    baux = np.zeros(34)
    baux[aux>=0.5] = 5
    baux[aux<0.5]  = 6
    baux = baux.astype(int)
    time_input = time_input_40
    while time_input < time_limit: #>40Hz inputs
        for i in range(34):
            DG_input_c.append(np.random.normal(time_input,time_sg_40,size=baux[i]))
            DG_input_label_c += list(np.repeat(166+i,baux[i]))
        time_input+=theta_rythm

    DG_input_a = np.round(np.concatenate(DG_input_a),2)
    DG_input_b = np.round(np.concatenate(DG_input_b),2)
    DG_input_c = np.round(np.concatenate(DG_input_c),2)

    EC2_list = [EC2_input_a, EC2_input_b, EC2_input_label_a, EC2_input_label_b]
    DG_list  = [DG_input_a, DG_input_b, DG_input_c, DG_input_label_a, \
                DG_input_label_b, DG_input_label_c]

    return EC2_list, DG_list

def make_external_inputs_ca1(simdur, seed,relative_time,time_initial,inputs):
    '''
    Here I generate the same inputs that I do in the previous function
    for not having the same peaks in EC3 than in EC2 (maybe this can be
    easily improved)
    simdur: simulation time
    seed :  seed for the generation of random numbers
    relative_time: relative time between MS and ECII inputs
    time_initial: time of the first generation of a MS inputs. Altough
    MS inputs are not generated in this function, this parameters is used
    as a reference.
    '''
    theta_rythm = 125 #ms
    ISI = 2*theta_rythm
    time_limit = simdur
    time_offset = relative_time

    # data extracted form Mizuseki et al 2009 fig 8. Results from gaussian fits
    ec2_4_diff, ec2_4_sg =    0.0,  9.54 #ms
    ec3_4_diff, ec3_4_sg = -60.04, 13.20 #ms
    dg_4_diff,  dg_4_sg  =  12.19, 16.51 #ms
    dg_40_diff, dg_40_sg =  39.38, 14.69 #ms

    np.random.seed(seed)
    ###########################################################################
    ''' EC2 inputs '''
    ###########################################################################
    time_input = time_offset + time_initial
    time_sg = ec2_4_sg

    EC2_inputs_mean   = []
    EC2_input_label_a = []
    EC2_input_label_b = []
    EC2_input_a, EC2_input_b = [],[]
    while time_input < time_limit:
        EC2_input_a.append(np.random.normal(time_input,time_sg,size=100))
        EC2_input_b.append(np.random.normal(time_input+theta_rythm,time_sg,size=100))
        EC2_input_label_a+= list(range(100))
        EC2_input_label_b+= list(range(100,200))
        EC2_inputs_mean.append(time_input)
        EC2_inputs_mean.append(time_input+theta_rythm)
        time_input+=ISI

    EC2_input_a = np.concatenate(EC2_input_a)
    EC2_input_b = np.concatenate(EC2_input_b)
    EC2_input_a = np.round(EC2_input_a,2)
    EC2_input_b = np.round(EC2_input_b,2)
    EC2_list = [EC2_input_a, EC2_input_b, EC2_input_label_a, EC2_input_label_b]
    ###########################################################################
    ''' DG inputs '''
    ###########################################################################
    time_input_4  = time_offset+time_initial+theta_rythm-dg_4_diff
    time_input_40 = time_offset+time_initial+theta_rythm-dg_40_diff
    time_sg_4     = dg_4_sg
    time_sg_40    = dg_40_sg

    DG_input_a = []
    DG_input_b = []
    DG_input_c = []

    time_input = time_input_4
    while time_input < time_limit: #4Hz inputs
        DG_input_a.append(np.random.normal(time_input,time_sg_4,size=83))
        DG_input_b.append(np.random.normal(time_input+theta_rythm,time_sg_4,size=83))
        time_input+=ISI

    aux = np.random.rand(34)
    baux = np.zeros(34)
    baux[aux>=0.5] = 5
    baux[aux<0.5]  = 6
    baux = baux.astype(int)
    time_input = time_input_40
    while time_input < time_limit: #>40Hz inputs
        for i in range(34):
            DG_input_c.append(np.random.normal(time_input,time_sg_40,size=baux[i]))
        time_input+=theta_rythm

    #######################################################################
    ''' EC3 inputs '''
    #######################################################################
    time_input = time_offset+time_initial-ec3_4_diff
    time_sg = ec3_4_sg

    EC3_inputs_mean   = []
    EC3_input_label_a = []
    EC3_input_label_b = []
    EC3_input_a, EC3_input_b = [],[]
    while time_input < time_limit:
        EC3_input_a.append(np.random.normal(time_input,time_sg,size=100))
        EC3_input_b.append(np.random.normal(time_input+theta_rythm,time_sg,size=100))
        EC3_input_label_a+= list(range(100))
        EC3_input_label_b+= list(range(100,200))
        EC3_inputs_mean.append(time_input)
        EC3_inputs_mean.append(time_input+theta_rythm)
        time_input+=ISI

    EC3_input_a = np.concatenate(EC3_input_a)
    EC3_input_b = np.concatenate(EC3_input_b)
    EC3_input_a = np.round(EC3_input_a,2)
    EC3_input_b = np.round(EC3_input_b,2)

    EC3_list = [EC3_input_a, EC3_input_b, EC3_input_label_a, EC3_input_label_b]

    ###########################################################################
    ''' pyr ca3 inputs '''
    ###########################################################################
    os.chdir(data_folder)
    file = data_file + str(inputs[0])+'_'+str(inputs[1])+'_'+str(inputs[2])+'.txt'

    x,y=[],[]
    for line in open(file):
        x.append( float(line.split(' ')[0]))
        y.append( int(line.split(' ')[1]))
    x,y = np.array(x),np.array(y)
    x = np.round(x,1)
    os.chdir(current_folder)
    pyrca3_list = [x,y]
    return EC3_list, pyrca3_list

def make_external_inputs(simdur, seed, relative_time, time_initial):
    '''
    This function creates all the inputs to the whole circuit
    of hippucampus. This will be used when a running the code
    where all the circuit will be

    simdur: simulation time
    seed :  seed for the generation of random numbers
    relative_time: relative time between MS and ECII inputs
    time_initial: time of the first generation of a MS inputs. Altough
    MS inputs are not generated in this function, this parameters is used
    as a reference.
    '''
    theta_rythm = 125 #ms
    ISI = 2*theta_rythm
    time_limit = simdur
    time_offset = relative_time

    # data extracted form Mizuseki et al 2009 fig 8. Results from gaussian fits
    ec2_4_diff, ec2_4_sg =    0.0,  9.54 #ms
    ec3_4_diff, ec3_4_sg = -60.04, 13.20 #ms
    dg_4_diff,  dg_4_sg  =  12.19, 16.51 #ms
    dg_40_diff, dg_40_sg =  39.38, 14.69 #ms

    np.random.seed(seed)
    ###########################################################################
    ''' EC2 inputs '''
    ###########################################################################
    time_input = time_offset + time_initial
    time_sg = ec2_4_sg

    EC2_inputs_mean   = []
    EC2_input_label_a = []
    EC2_input_label_b = []
    EC2_input_a, EC2_input_b = [],[]
    while time_input < time_limit:
        EC2_input_a.append(np.random.normal(time_input,time_sg,size=100))
        EC2_input_b.append(np.random.normal(time_input+theta_rythm,time_sg,size=100))
        EC2_input_label_a+= list(range(100))
        EC2_input_label_b+= list(range(100,200))
        EC2_inputs_mean.append(time_input)
        EC2_inputs_mean.append(time_input+theta_rythm)
        time_input+=ISI

    EC2_input_a = np.concatenate(EC2_input_a)
    EC2_input_b = np.concatenate(EC2_input_b)
    EC2_input_a = np.round(EC2_input_a,2)
    EC2_input_b = np.round(EC2_input_b,2)

    ###########################################################################
    ''' DG inputs '''
    ###########################################################################
    time_input_4  = time_offset+time_initial+theta_rythm-dg_4_diff
    time_input_40 = time_offset+time_initial+theta_rythm-dg_40_diff
    time_sg_4     = dg_4_sg
    time_sg_40    = dg_40_sg

    DG_input_label   = np.arange(200)
    DG_input_label_a = DG_input_label[:83]
    DG_input_label_b = DG_input_label[83:166]
    DG_input_label_c = DG_input_label[166:]

    DG_input_label_a,  DG_input_a = [], []
    DG_input_label_b,  DG_input_b = [], []
    DG_input_label_c,  DG_input_c = [], []

    time_input = time_input_4
    while time_input < time_limit: #4Hz inputs
        DG_input_a.append(np.random.normal(time_input,time_sg_4,size=83))
        DG_input_b.append(np.random.normal(time_input+theta_rythm,time_sg_4,size=83))
        DG_input_label_a += list(range(83))
        DG_input_label_b += list(range(83,166))
        time_input+=ISI

    aux = np.random.rand(34)
    baux = np.zeros(34)
    baux[aux>=0.5] = 5
    baux[aux<0.5]  = 6
    baux = baux.astype(int)
    time_input = time_input_40
    while time_input < time_limit: #>40Hz inputs
        for i in range(34):
            DG_input_c.append(np.random.normal(time_input,time_sg_40,size=baux[i]))
            DG_input_label_c += list(np.repeat(166+i,baux[i]))
        time_input+=theta_rythm

    DG_input_a = np.round(np.concatenate(DG_input_a),2)
    DG_input_b = np.round(np.concatenate(DG_input_b),2)
    DG_input_c = np.round(np.concatenate(DG_input_c),2)

    ############################################################################
    ''' EC3 inputs '''
    ############################################################################
    time_input = time_offset+time_initial-ec3_4_diff
    time_sg = ec3_4_sg

    EC3_inputs_mean   = []
    EC3_input_label_a = []
    EC3_input_label_b = []
    EC3_input_a, EC3_input_b = [],[]
    while time_input < time_limit:
        EC3_input_a.append(np.random.normal(time_input,time_sg,size=100))
        EC3_input_b.append(np.random.normal(time_input+theta_rythm,time_sg,size=100))
        EC3_input_label_a+= list(range(100))
        EC3_input_label_b+= list(range(100,200))
        EC3_inputs_mean.append(time_input)
        EC3_inputs_mean.append(time_input+theta_rythm)
        time_input+=ISI

    EC3_input_a = np.concatenate(EC3_input_a)
    EC3_input_b = np.concatenate(EC3_input_b)
    EC3_input_a = np.round(EC3_input_a,2)
    EC3_input_b = np.round(EC3_input_b,2)

    ############################################################################
    EC2_list = [[EC2_input_a, EC2_input_b], [EC2_input_label_a, EC2_input_label_b]]
    DG_list  = [[DG_input_a, DG_input_b, DG_input_c], [DG_input_label_a, \
                DG_input_label_b, DG_input_label_c]]
    EC3_list = [[EC3_input_a, EC3_input_b], [EC3_input_label_a, EC3_input_label_b]]

    return EC2_list, DG_list, EC3_list

################################################################################
''' El cacho que hay a continuacion corresponde con el metodo previo que implemente
para modificar el burst_noise durante la simulacion. Sin embargo, dada la incompataibilidad
del h.BurstStim con el gaussian generator me he visto en la obligacion de generar los inputs
con python, al menos solo los correspondientes con el burst '''
'''
nslabel = np.cumsum(nslabel).astype(int)
nslabel = np.insert(nslabel, 0, 0, axis=0)    #return nsl_sep_180, nsl_sep_360, nsl_ec2_180, nsl_ec2_360,nsl_ec3_180, nsl_ec3_360, nsl_dg_burst, nsl_dg_regular

################################################################################
# auxiliar dg burst source for modifying noise during simulation
# The same seed for the generation of the inputs will be used
################################################################################
aux_seadlist = nrseadlist[nslabel[6]:nslabel[7]]
for i in range(dg_burst_neurons):
    ns = h.BurstStim()

    ns.interval = dg_burst_int
    ns.number   = dg_burst_num
    ns.start    = dg_burst_start
    ns.noise    = dg_burst_noise_var
    ns.burstint = dg_burst_bint
    ns.burstlen = dg_burst_blen

    sead = aux_seadlist[i]
    rs = RandomStream(sead,'gaussian')
    rs.set_generator()
    ns.noiseFromRandom(rs.r)
    rs.start()

    nslist.append(ns)
    nrlist.append(rs)
    nrseadlist.append(rs.sead)

nslabel = np.append(nslabel, nslabel[-1]+dg_burst_neurons)
'''
