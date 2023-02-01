'''#############################################################################
Date 06/06/2022

Functions file
#############################################################################'''

import numpy as np
from neuron import h
import random
from collections import Counter

class Population:
    "Population of cells"
    # cell_type -- pyr, bas, olm
    # n -- number of cells in the population
    # x, y, z -- initial position for the first Cell
    # dx -- an increment of the x-position for the cell location
    # amp, dur, delay -- parameters for the IClamp in the soma
    def __init__(self, cell_type, n , x, y, z, dx, amp, dur, delay):
        self.cell = [] # List of cells in the population
        self.nc   = [] # NetCon list for recording spikes
        self.n    = n  # number of cells
        self.x    = x
        self.y    = y
        self.z    = z
        self.ltimevec = [] # list of Vectors for recording spikes, one per cell
        self.lidvec = []
        #self.tvec = h.Vector()
        #self.idvec = h.Vector()
        self.nssidx = {}
        self.nseidx = {}
        self.ncsidx = {}
        self.nceidx = {}
        for i in range(n):
            #self.cell.append(cell_type(x+i*dx,y,z,gGID))
            self.cell.append(cell_type(x+i*dx,y,z,i))
            self.cell[-1].somaInj.amp   = amp
            self.cell[-1].somaInj.dur   = dur
            self.cell[-1].somaInj.delay = delay
            
            self.cell[-1].spike_detector.record(self.cell[-1].spike_times)
            #self.nc.append(h.NetCon(self.cell[-1].soma(0.5)._ref_v, None, sec=self.cell[-1].soma))
            #self.ltimevec.append(h.Vector()) #NB: each NetCon gets own Vectors for recording. needed to avoid multithreading crash
            #self.lidvec.append(h.Vector())
            #self.nc[-1].record(self.ltimevec[-1],self.lidvec[-1],i)
            
        for section in self.cell[0].sec_list: # define a vector per each comparment/secion to later save the synaptic currents
            self.__dict__[section+"_isyn"] = {}

        for syn in self.cell[0].syn_list:
            self.__dict__[syn+"_mean"] = []
            self.__dict__[syn+"_std"]  = []

    def set_r(self, syn, r):
        for c in self.cell:
            c.__dict__[syn].syn.r = r

def init_NetStims(nrl,nrlsead):
    for i in range(len(nrl)):
        rds = nrl[i]
        sead = nrlsead[i]
        rds.MCellRan4(sead,sead)
        rds.negexp(1)

def make_NetStims(po, syn, delay,w,ISI, time_limit, seed):
    ''' noise for background '''
    #self.netstims_tvec[po.cell[0].label+"_"+syn] = h.Vector()
    #self.netstims_ivec[po.cell[0].label+"_"+syn] = h.Vector()
    nsl, ncl, nrl, nrlsead = [],[],[],[]
    for i in range(po.n):
        cell = po.cell[i]

        ns = h.NetStim()
        ns.interval = ISI
        ns.noise = 1
        ns.number = (1e3/ISI)*time_limit
        ns.start = 0

        nc = h.NetCon(ns,cell.__dict__[syn].syn)
        #print(po.pgidlist[i]-po.naux,po.naux)
        nc.delay = delay[i]
        nc.weight[0] = w
        #if pc.id() == 0:
        #    if self.record_noise:
        #        nc.record(self.netstims_tvec[po.cell[0].label+"_"+syn], self.netstims_ivec[po.cell[0].label+"_"+syn],i)

        rds = h.Random()
        rds.negexp(1)            # set random # generator using negexp(1) - avg interval in NetStim
        sead = seed+i #to avoid sead = 0 and having an element to control
        rds.MCellRan4(sead,sead) # seeds are in order, shouldn't matter
        ns.noiseFromRandom(rds)  # use random # generator for this NetStim

        nsl.append(ns) # netstims list
        ncl.append(nc) # netcons list
        nrl.append(rds) # generator list
        nrlsead.append(sead) #sead list

    return nsl, ncl, nrl, nrlsead, sead+1

def set_connectivity(nsrc, ntrg, conv):
    conn = make_conn(nsrc,ntrg,conv)
    return conn

def make_conn(preN, postN, conv)
    conn = np.zeros((postN,conv),dtype=np.int16)
    for i in range(postN):
        conn[i,:]=random.sample(range(preN),conv)
    return conn

def set_connections(src,trg,syn,delay,w,conn): # mirar esto con detalle
    #conn = self.make_conn(src.n,trg.ng,conv, trg)
    #conn = self.make_conn(src.n,trg.n,conv, trg)
    #conn = conn[ np.arange(pc.id(), trg.n, pc.nhost()),:]
    nc = []
    for post_id, all_pre in enumerate(conn):
        for j, pre_id in enumerate(all_pre):
            nc.append(h.NetCon(src.cell[pre_id].soma(0.5)._ref_v, trg.cell[post_id].__dict__[syn].syn, 0, delay[j,post_id], w, sec=src.cell[pre_id].soma))	
            
    return nc
    