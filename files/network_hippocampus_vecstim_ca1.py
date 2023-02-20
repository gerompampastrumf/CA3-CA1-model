"""
IMPORTANT!!!! 

In this version I use Pyr_CA3 for pyramidals in CA1. This 
is to understand better what can be happening. 

So, in this situation I will just focus on the connectios between the cck and pyramidal 
and leaving the rest of the paramteres related to the populations unchanged. 

For example, conductance from inputs from EC3 to the different cells will reamin as it was in EC2. 
"""

from neuron import h
import numpy as np
from spectrum import *
import random
#import functions
from neurons import *
import sys
import os
# import warnings
# from external_inputs import *
sys.path.append('/home/jaime/Desktop/hippocampus/files/')
import file_management

#warnings.simplefilter(action='ignore', category=FutureWarning)
#warnings.simplefilter(action='ignore', category=RuntimeWarning)

#h.nrnmpi_init()
pc = h.ParallelContext()

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

class RandomStream():
    def __init__(self, sead, generator):
        self.generator = generator
        self.r = h.Random()
        self.sead = sead
        self.start()

    def set_generator(self):
        if self.generator == 'gaussian':
            g = self.r.normal(0,1)
        if self.generator == 'exponential':
            g = self.r.negexp(1)
        return g

    def start(self):
        return self.r.MCellRan4(self.sead,self.sead)

    def repick(self):
        return self.r.repick()

class Population:
    "Population of cells"
    # cell_type -- pyr, bas, olm
    # n -- number of cells in the population
    # x, y, z -- initial position for the first Cell
    # dx -- an increment of the x-position for the cell location
    # amp, dur, delay -- parameters for the IClamp in the soma
    def __init__(self, cell_type, pgidlist, n , x, y, z, dx, amp, dur, delay, nseg, record_mp, pyr=False, naux=0,
    x_position = None, y_position = None):
        self.cell = [] # List of cells in the population
        self.nc   = [] # NetCon list for recording spikes
        self.pgidlist = pgidlist
        self.ng   = len(self.pgidlist) # number of cells in host
        self.n    = n  # number of cells
        self.x    = x
        self.y    = y
        self.z    = z
        #self.ltimevec = h.List() # list of Vectors for recording spikes, one per cell
        #self.lidvec = h.List()
        #self.tvec = h.Vector()
        #self.idvec = h.Vector()
        self.nssidx = {}
        self.nseidx = {}
        self.ncsidx = {}
        self.nceidx = {}
        self.ncd  = [] # only for pyramidal

        self.naux = int(naux) # to distinguish gid from i-th position of the neruon in the network
        self.nc_conn = []
        if pyr:
            for i in self.pgidlist: #range(n):
                #self.cell.append(cell_type(x+i*dx,y,z,gGID))
                self.cell.append(cell_type(x_position[i],y_position[i],z,i,nseg, record_mp))
                self.cell[-1].somaInj.amp   = amp
                self.cell[-1].somaInj.dur   = dur
                self.cell[-1].somaInj.delay = delay
                pc.cell(i, self.cell[-1].spike_detector)
        else:
            for i in self.pgidlist: #range(n):
                #self.cell.append(cell_type(x+i*dx,y,z,gGID))
                self.cell.append(cell_type(x+(i-self.naux)*dx,y,z,i, nseg, record_mp))
                self.cell[-1].somaInj.amp   = amp
                self.cell[-1].somaInj.dur   = dur
                self.cell[-1].somaInj.delay = delay
                pc.cell(i, self.cell[-1].spike_detector)
        for section in self.cell[0].sec_list: # define a vector per each comparment/secion to later save the synaptic currents
            self.__dict__[section+"_isyn"] = {}

        for syn in self.cell[0].syn_list:
            self.__dict__[syn+"_mean"] = []
            self.__dict__[syn+"_std"]  = []

    def set_r(self, syn, r):
        for c in self.cell:
            c.__dict__[syn].syn.r = r

''' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'''
class Network:

    def __init__(self, weights_inputs_neurons, nsyns_inputs_neurons, syn_inputs_neurons,
               weights_neurons_neurons, nsyns_neurons_neurons, syn_neurons_neurons,
               weights_noise_neurons, delay_neurons_neurons, delay_inputs_neurons,
               connections=True,
               DoMakeNoise=True,
               DoMakeExternalInputs=True,
               MakeNetStim=True,
               MakeCellStim=True,
               MakeCellStim_python = True,
               bseed=1234,
               wseed=4321,
               iseed=3421,
               jseed=2134,
               scale=1.0,
               SaveConn=False,
               noise_burst=0.5,
               create_ca1_network=False,
               record_noise=False,
               ISI_noise = 1.0,
               jitter=0.2,
               delays_in_noise=True,
               delays_in_inputs=True,
               delays_in_connections=True,
               delay_mean_noise=2.0,
               delay_mean_inputs=2.0,
               delay_mean_connections=2.0,
               save_input_connection=False,
               topology="rectangular",
               inputs_folder = '/home/jaime/Desktop/',
               burst_level_label = '10',
               external_inputs_label="0",
               background_noise_label="0",
               save_external_inputs = False,
               n_pyr_ca3 = 800,
               n_bas_ca3 = 100,
               n_olm_ca3 = 30,
               n_pyr_ca1 = 800,
               n_bas_ca1 = 100,
               n_olm_ca1 = 30,
               n_cck_ca1 = 70,
               nseg_pyr = 5, 
               record_mp = {"soma":True, "Bdend":False, "Adend1":False, "Adend2":False, "Adend3":False} ):

        self.nseg_pyr = seg_pyr
        print("Setting Cells")
        self.connections = connections
        self.weights_inputs_neurons  = weights_inputs_neurons
        self.nsyns_inputs_neurons    = nsyns_inputs_neurons
        self.syn_inputs_neurons      = syn_inputs_neurons
        self.weights_neurons_neurons = weights_neurons_neurons
        self.nsyns_neurons_neurons   = nsyns_neurons_neurons
        self.syn_neurons_neurons     = syn_neurons_neurons
        self.weights_noise_neurons   = weights_noise_neurons
        self.delay_neurons_neurons   = delay_neurons_neurons
        self.delay_inputs_neurons    = delay_inputs_neurons
        self.record_mp = record_mp ## new 17/01

        self.ISI_noise = ISI_noise
        self.create_ca1_network = create_ca1_network
        self.record_noise = record_noise
        self.save_input_connection = save_input_connection

        self.delays_in_noise  = delays_in_noise
        self.delays_in_inputs = delays_in_inputs
        self.delays_in_connections = delays_in_connections
        self.delay_mean_noise  = delay_mean_noise
        self.delay_mean_inputs = delay_mean_inputs
        self.delay_mean_connections = delay_mean_connections

        self.topology = topology
        self.n_pyr_ca3 = n_pyr_ca3
        self.n_bas_ca3 = n_bas_ca3
        self.n_olm_ca3 = n_olm_ca3
        self.n_pyr_ca1 = n_pyr_ca1
        self.n_bas_ca1 = n_bas_ca1
        self.n_olm_ca1 = n_olm_ca1
        self.n_cck_ca1 = n_cck_ca1
        self.set_pyramidal_topology(type=self.topology) #this will change the number of neurons

        # self.n_pop = [self.n_pyr_ca3, self.n_bas_ca3, self.n_olm_ca3]
        # self.cell_labels = ["pyr_ca3", "bas_ca3", "olm_ca3"]
        # if self.create_ca1_network:
        self.n_pop = [self.n_pyr_ca1, self.n_bas_ca1, self.n_olm_ca1, self.n_cck_ca1]
        self.cell_labels = ["pyr_ca1", "bas_ca1", "olm_ca1", "cck_ca1"]

        self.n_pop = np.cumsum(self.n_pop)
        self.n = self.n_pop[-1]

        self.set_gids()
        self.create_populations()

        self.bseed = bseed # seed for background noise should be ordened in order to keep same realization of when running ca3 and ca1 together
        self.wseed = wseed # seed for wiring
        self.iseed = iseed # seed for external inputs (EC2, EC3, MS and DG)
        self.jseed = jseed # wiring for external inputs
        self.dseed = 1993  # for distribution of delays on connections

        self.jitter = jitter # sigma of the distribution of delays
        self.DoMakeNoise = DoMakeNoise
        self.DoMakeExternalInputs = DoMakeExternalInputs
        self.MakeNetStim  = MakeNetStim
        self.MakeCellStim = MakeCellStim # added on 14/03/2022
        self.MakeCellStim_python = MakeCellStim_python # added on 16/07/2022
        self.SCGain   = 0.5
        self.SaveConn = SaveConn
        self.save_external_inputs = save_external_inputs

        self.time_initial = 50.0 # time in which external inputs are generated
        self.theta_rythm = 125.0 # 8Hz overall activit

        self.noise_burst = noise_burst

        self.inputs_folder = inputs_folder
        self.burst_level_label = burst_level_label
        self.external_inputs_label = external_inputs_label
        self.background_noise_label = background_noise_label # choose the correct pyr_ca3 file.

        self.set_delays_distribution()
        if self.connections:
            print("Setting Connections")
            self.f_conns()

    def set_pyramidal_topology(self, type):
        dx = 50.0
        if type == "circular":
            npyr = np.array([6*i for i in range(1,17)])[2:]
            radius = np.array([dx*i for i in range(1,17)])[2:]
            # self.pyr_ca3_x, self.pyr_ca3_y = [],[]
            # for i,n in enumerate(npyr):
            #     theta = np.linspace(0,2*np.pi, n+1)[:-1]
            #     self.pyr_ca3_x.append(radius[i]*np.cos(theta))
            #     self.pyr_ca3_y.append(radius[i]*np.sin(theta))
            # self.pyr_ca3_x = np.concatenate(self.pyr_ca3_x)
            # self.pyr_ca3_y = np.concatenate(self.pyr_ca3_y)
            # self.n_pyr_ca3 = np.sum(npyr)

            self.pyr_ca1_x, self.pyr_ca1_y = [],[]
            for i,n in enumerate(npyr):
                theta = np.linspace(0,2*np.pi, n+1)[:-1]
                self.pyr_ca1_x.append(radius[i]*np.cos(theta))
                self.pyr_ca1_y.append(radius[i]*np.sin(theta))
            self.pyr_ca1_x = np.concatenate(self.pyr_ca1_x)
            self.pyr_ca1_y = np.concatenate(self.pyr_ca1_y)
            self.n_pyr_ca1 = np.sum(npyr)

        elif type == "rectangular":
            # self.pyr_ca3_x = np.concatenate([np.arange(0,40,1)*dx for i in range(20)])-1000.0
            # self.pyr_ca3_y = np.concatenate([np.ones(40)*dx*i for i in range(20)])-500.0
            self.pyr_ca1_x = np.concatenate([np.arange(0,40,1)*dx for i in range(20)])-1000.0
            self.pyr_ca1_y = np.concatenate([np.ones(40)*dx*i for i in range(20)])-500.0
        else:
            pass
            print("All neurons placed in a row")

    def set_gids(self):
        """ set the gidlist on this host. """
        self.gidlist=[]
        for i in range(pc.id(), self.n ,pc.nhost()):
            self.gidlist.append(i)
            pc.set_gid2node(i, int(pc.id()) )

    def create_populations(self):
        # gid_pyr_ca3, gid_bas_ca3, gid_olm_ca3 = [],[],[]
        gid_pyr_ca1, gid_bas_ca1, gid_olm_ca1, gid_cck_ca1 = [],[],[],[] #,[],[],[]

        for i in self.gidlist:
            if np.logical_and(i>=0,i<self.n_pop[0]): gid_pyr_ca1.append(i)
            if np.logical_and(i>=self.n_pop[0],i<self.n_pop[1]): gid_bas_ca1.append(i)
            if np.logical_and(i>=self.n_pop[1],i<self.n_pop[2]): gid_olm_ca1.append(i)
            if np.logical_and(i>=self.n_pop[2],i<self.n_pop[3]): gid_cck_ca1.append(i)

        # esto hay que retocarlo un poco.
        # self.pyr_ca3 = Population(cell_type=PyrAdr_CA3, pgidlist=gid_pyr_ca3, n=self.n_pyr_ca3, x= 0, y=0, z=0, dx=50.0, amp= 40e-3, dur=1e9, delay=2*h.dt, nseg = self.nseg_pyr, naux = 0, pyr=True, x_position=self.pyr_ca3_x,y_position=self.pyr_ca3_y) # Estas corrientes se han anadido
        # self.bas_ca3 = Population(cell_type=Bwb,        pgidlist=gid_bas_ca3, n=self.n_bas_ca3, x=10, y=0, z=0, dx=50.0, amp= 0,     dur=1e9, delay=2*h.dt, nseg = 1, naux=self.n_pop[0]) # para obtener el baseline state
        # self.olm_ca3 = Population(cell_type=Ow,         pgidlist=gid_olm_ca3, n=self.n_olm_ca3, x=20, y=0, z=0, dx=50.0, amp=-25e-3, dur=1e9, delay=2*h.dt, nseg = 1, naux=self.n_pop[1]) # dur tiene que se mayor o igual que el tiempo de simulacion
        # self.cell_ca3 = self.pyr_ca3.cell+self.bas_ca3.cell+self.olm_ca3.cell
        # self.pop_ca3  = [self.pyr_ca3, self.bas_ca3, self.olm_ca3]

        self.pyr_ca1 = Population(cell_type=PyrAdr_CA3, pgidlist=gid_pyr_ca1, n=self.n_pyr_ca1, x=0,  y=0, z=0, dx=50, amp= 40e-3, dur=1e9, delay=2*h.dt, nseg=self.nseg_pyr, naux=0, record_mp = self.record_mp) # Estas corrientes se han anadido
        self.bas_ca1 = Population(cell_type=Bwb,        pgidlist=gid_bas_ca1, n=self.n_bas_ca1, x=10, y=0, z=0, dx=50, amp= 0,     dur=1e9, delay=2*h.dt, nseg=1, naux=self.n_pop[0], record_mp = self.record_mp) # para obtener el baseline state
        self.olm_ca1 = Population(cell_type=Ow,         pgidlist=gid_olm_ca1, n=self.n_olm_ca1, x=20, y=0, z=0, dx=50, amp=-25e-3, dur=1e9, delay=2*h.dt, nseg=1, naux=self.n_pop[1], record_mp = self.record_mp) # dur tiene que se mayor o igual que el tiempo de simulacion
        self.cck_ca1 = Population(cell_type=Cck_cell,   pgidlist=gid_cck_ca1, n=self.n_cck_ca1, x=30, y=0, z=0, dx=50, amp= 10e-3, dur=1e9, delay=2*h.dt, nseg=1, naux=self.n_pop[2], record_mp = self.record_mp) # dur tiene que se mayor o igual que el tiempo de simulacion
        self.cell_ca1 = self.pyr_ca1.cell+self.bas_ca1.cell+self.olm_ca1.cell+self.cck_ca1.cell
        self.pop_ca1  = [self.pyr_ca1, self.bas_ca1, self.olm_ca1, self.cck_ca1]

    def set_noise_inputs(self,simdur):
        if self.DoMakeNoise:
            if self.MakeNetStim:
                self.make_all_NetStims(simdur)
            else:
                print("Function to load external noise not implemented yet")
        else:
            print("Simulation without external noise")

        if self.DoMakeExternalInputs:
            print('1')
            if self.MakeCellStim:
                print('2')
                if self.MakeCellStim_python:
                    print('3')
                    self.make_all_CellStims_python(simdur)
                else:
                    self.make_all_CellStims(simdur)
                    print('4')
            else:
                self.load_all_CellStims(simdur)
                print('5')
        else:
            print('6')
            print("Simulation without external inputs")

        print("Done!")

    #both should be called @ beginning of each sim
    def init_NetStims(self):
        # h.mcell_ran4_init(self.iseed)
        for i in range(len(self.nrl)):
            rds = self.nrl[i]
            sead = self.nrlsead[i]
            rds.MCellRan4(sead,sead)
            rds.negexp(1)

    def init_CellStims(self):
        print('---------------------')
        print('INIT CELLSTIMS')
        # print(len(self.nrl_))
        print('---------------------')
        for i in range(len(self.nrl_)):
            rds = self.nrl_[i]
            rds.start()
            rds.set_generator()
            # print(i,rds.sead)
        # print('------------')

    def make_NetStims(self, po, syn, delay, w, ISI, time_limit, add_to_sead): #added delay in 24th
        ''' noise for background '''
        po.nssidx[syn] = len(self.nsl) #index into net.nsl
        po.ncsidx[syn] = len(self.ncl) #index into net.
        self.netstims_tvec[po.cell[0].label+"_"+syn] = h.Vector()
        self.netstims_ivec[po.cell[0].label+"_"+syn] = h.Vector()
        for i in range(po.ng):
            cell = po.cell[i]

            ns = h.NetStim()
            ns.interval = ISI
            ns.noise = 1
            ns.number = (1e3/ISI)*time_limit
            ns.start = 0

            nc = h.NetCon(ns,cell.__dict__[syn].syn)
            #print(po.pgidlist[i]-po.naux,po.naux)
            nc.delay = delay[po.pgidlist[i]-po.naux]
            nc.weight[0] = w
            if pc.id() == 0:
                if self.record_noise:
                    nc.record(self.netstims_tvec[po.cell[0].label+"_"+syn], self.netstims_ivec[po.cell[0].label+"_"+syn],i)

            rds = h.Random()
            rds.negexp(1)            # set random # generator using negexp(1) - avg interval in NetStim
            sead = po.pgidlist[i] + self.bseed + add_to_sead #to avoid sead = 0 and having an element to control
            rds.MCellRan4(sead,sead) # seeds are in order, shouldn't matter
            ns.noiseFromRandom(rds)  # use random # generator for this NetStim

            #ns.start = rds.discunif(0,1e3) # start inputs random time btwn 0-1e3 ms to avoid artificial sync
            #rds.MCellRan4(sead,sead) # reinit rand # generator
            #nc_list.append(nc)
            self.nsl.append(ns)
            self.ncl.append(nc)
            cell.ncl.append(nc)
            self.nrl.append(rds) #esta lista de random es la que no entiendo muy bien
            self.nrlsead.append(sead)

        add_to_sead = add_to_sead + po.n

        po.nseidx[syn] = len(self.nsl)-1
        po.nceidx[syn] = len(self.ncl)-1
        return add_to_sead

    def set_delays_distribution(self):
        # I generate delays for every input and noise source
        self.sigma_noise  = 0.0
        self.sigma_inputs = 0.0
        self.sigma_connections = 0.0

        if self.delays_in_noise:
            self.sigma_noise = 1.0
        if self.delays_in_inputs:
            self.sigma_inputs = 1.0
        if self.delays_in_connections:
            self.sigma_connections = 1.0

        rdtmp = self.dseed
        rdtmp = self.set_delays_to_synaptic_noise(delay_mean=self.delay_mean_noise, sigma=self.sigma_noise,seed=rdtmp)
        rdtmp = self.set_delays_to_synaptic_inputs( sigma=self.sigma_inputs, seed=rdtmp)
        rdtmp = self.set_delays_to_synpatic_connections( sigma=self.sigma_connections,seed=rdtmp)

    def set_delays_to_synpatic_connections(self,sigma, seed):
        ''' delays neurons -> neurons connections '''
        print("Setting the delays: neurons to neurons")
        self.delays_neurons_neurons = dict.fromkeys(self.weights_neurons_neurons.keys() )
        for key in self.delays_neurons_neurons.keys():
            if key.endswith("ca1"):
                np.random.seed(seed)
                nsyn = self.nsyns_neurons_neurons[key]
                po = self.__dict__[key.split("_to_")[1]]
                self.delays_neurons_neurons[key] = []
                delay_mean = self.delay_neurons_neurons[key][0]
                for n in nsyn:
                    self.delays_neurons_neurons[key].append( np.abs(delay_mean+sigma*np.random.normal(0,self.jitter,(n,po.n))) )
                    seed += 1
        return seed

    def set_delays_to_synaptic_inputs(self, sigma, seed):
        ''' delays inputs -> neurons connections '''
        print("Setting the delays: inputs to neurons")
        self.delays_inputs_neurons = dict.fromkeys(self.weights_inputs_neurons.keys())
        for key in self.delays_inputs_neurons.keys():
            if key.endswith("ca1"):
                np.random.seed(seed)
                nsyn = self.nsyns_inputs_neurons[key]
                po = self.__dict__[key.split("_to_")[1]]
                self.delays_inputs_neurons[key] = []
                delay_mean = delay_inputs_neurons[key][0]
                for n in nsyn:
                    self.delays_inputs_neurons[key].append( np.abs(delay_mean+sigma*np.random.normal(0,self.jitter,(n,po.n))) )
                    seed += 1
        return seed

    def set_delays_to_synaptic_noise(self,delay_mean,sigma,seed):
        ''' delays inputs noise background -> neurons '''
        print("Setting the delays: noise to neurons")
        self.delays_noise_neurons = dict.fromkeys(self.cell_labels)
        for key in self.cell_labels:
            if key.endswith("ca1"):
                np.random.seed(seed)
                po = self.__dict__[key]
                if key.split("_")[0] == "pyr":
                    keys = ["somaAMPA_noise","somaGABA_noise","Adend3AMPA_noise","Adend3GABA_noise","Adend3NMDA_noise"]
                else:
                    keys = ["somaAMPA_noise","somaGABA_noise"]
                self.delays_noise_neurons[key] = dict.fromkeys(keys)
                for k in keys:
                    self.delays_noise_neurons[key][k] = np.abs(delay_mean+sigma*np.random.normal(0.0,self.jitter,po.n))
                seed+=1
        return seed

    def make_all_NetStims(self,simdur): #,rdmseed):
        print("Making NetStims")
        self.netstims_tvec = {} # saving all
        self.netstims_ivec = {}
        self.nsl = [] # NetStim List
        self.ncl = [] # NetCon List
        self.nrl = [] # Random List for NetStims
        self.nrlsead = [] #List of seeds for NetStim random

        print("Making Background Noise")
        rdtmp = 0 # add_to_sead value - incremented in make_NetStims
        for cell in self.cell_labels:
            print("to "+cell.split("_")[0]+" "+cell.split("_")[1])
            po = self.__dict__[cell]
            for syn in self.weights_noise_neurons[cell].keys():
                delays = self.delays_noise_neurons[cell][syn]
                weight = self.weights_noise_neurons[cell][syn]
                rdtmp  = self.make_NetStims(po=po, syn=syn, delay=delays, w=weight,  ISI=self.ISI_noise,  time_limit=simdur, add_to_sead=rdtmp)#, sead=rdtmp)

    def set_CellStims_connections(self, ncellstims):
        keys = ["sep_180","sep_360","ec3_180","ec3_360","pyr_ca3"]
        inputs = dict.fromkeys(keys)
        for i, key in enumerate(keys):
            inputs[key] = ncellstims[i]

        random.seed(self.jseed) # for wiring
        keys = self.weights_inputs_neurons.keys()
        self.connectivity_inputs_neurons = dict.fromkeys(keys)
        for key in keys:
            if key.endswith("ca1"):
                ninputs = inputs[key.split("_to_")[0]]
                po = self.__dict__[key.split("_to_")[1]]
                # n = self.nsyns_inputs_neurons[key]
                # self.connectivity_inputs_neurons[key] = self.set_connectivity(ninputs, po.n, n)
                self.connectivity_inputs_neurons[key] = []
                for n in self.nsyns_inputs_neurons[key]:
                    self.connectivity_inputs_neurons[key].append( self.set_connectivity(ninputs, po.n, n) )

    def make_all_CellStims_python(self,simdur):
        '''
        This function load the external inputs created by  make_external_inputs_python
        '''
        self.ncl_ =[]
        print("Making external inputs")
        self.idvec, self.tvec, self.nsl_n, labels, self.input_time = make_external_inputs_python(simdur, self.iseed, pburst=0.5)

        print("Making cellstims connectionss")
        self.set_CellStims_connections(self.nsl_n)

        print("Creating dictionary with all parameters")
        keys1 = self.weights_inputs_neurons.keys()
        keys2 = ["idv", "spikes", "population", "synapses", "delays", "weights", "connectivity"]
        self.external_inputs_data = dict.fromkeys(keys1)
        for key in keys1:
            if key.endswith("ca1"):
                self.external_inputs_data[key] = dict.fromkeys(keys2)
                inputs = key.split("_to_")[0]
                po = self.__dict__[key.split("_to_")[1]]
                self.external_inputs_data[key]["idv"]          = self.idvec[inputs]
                self.external_inputs_data[key]["spikes"]       = self.tvec[inputs]
                self.external_inputs_data[key]["population"]   = po
                self.external_inputs_data[key]["synapses"]     = self.syn_inputs_neurons[key]
                self.external_inputs_data[key]["delays"]       = self.delays_inputs_neurons[key]
                self.external_inputs_data[key]["weights"]      = self.weights_inputs_neurons[key]
                self.external_inputs_data[key]["connectivity"] = self.connectivity_inputs_neurons[key]

    def make_all_CellStims(self,simdur):
        '''
        This function load the external inputs created by make_new_external_inputs
        and dg_burst_inputs. The ones generated by the first function are then
        connected to their correspondent populations while burst inputs are just
        loaded to be further sent to the main file to create the connections.
        This is necessary for the effective VecStim implementation.
        '''
        print("Making CellStims")
        self.nsl_  = [] #CellStims List
        self.nsl_n = []
        self.ncl_  = [] #NetCon List
        self.input_labels = []
        self.input_time = []  #to further generate the mean of each spike inputs
        print("Making external inputs")
        sead, self.nsl_, self.nrl_, self.nsl_n, self.input_labels, self.input_time =  make_new_external_inputs(self.time_initial,self.iseed)
        # idvec_dg_burst, tvec_dg_burst, ncells_burst = dg_burst_inputs(sead,self.noise_burst,simdur)
        # self.nsl_.append(ncells_burst)
        # self.nsl_n.append(ncells_burst)

        keys = ["sep_180","sep_360","ec3_180", "ec3_360","pyr_ca3"]
        self.tvec  = dict.fromkeys(keys)
        self.idvec = dict.fromkeys(keys)
        if pc.id() == 0: # only record in host 0
            if self.save_external_inputs:
                for key in keys[:-1]:
                    self.tvec[key] = h.Vector()
                    self.idvec[key] = h.Vector()

        print("Making cellstims connectionss")
        print("-----------------------------")
        self.set_CellStims_connections(self.nsl_n)

        inputs = dict.fromkeys(keys)
        for i, key in enumerate(keys):
            inputs[key] = self.nsl_[i]
        print("Making external inputs connections")
        print("----------------------------------")

        for key in self.weights_inputs_neurons.keys():
            if key.endswith("ca1"):
                if key.split("_to_")[0] != "dg_burst":
                    # print(key)
                    inp     = inputs[key.split("_to_")[0]]
                    po      = self.__dict__[key.split("_to_")[1]]
                    delays  = self.delays_inputs_neurons[key]
                    syn     = self.syn_inputs_neurons[key]
                    weights = self.weights_inputs_neurons[key]
                    conn    = self.connectivity_inputs_neurons[key]
                    self.set_input_connection(inp, po, syn, delays, weights, conn)

        # Save external inputs
        if pc.id() == 0:
            if self.save_external_inputs:
                for key in list(inputs.keys())[:-1]:
                    self.record_external_inputs(inputs[key], self.tvec[key], self.idvec[key])
                    #self.record_external_inputs( self.tvec[key], self.idvec[key])

    def record_external_inputs(self, CellStim, tvec, idvec):
        nil = h.Section(name='nil')
        nil.insert('hh')
        syn_ = h.ExpSyn(nil(0.5))
        for i in range(len(CellStim)):
            nc = h.NetCon(CellStim[i],syn_)
            nc.record(tvec,idvec,i)

    def set_input_connection(self,ns_src,trg,syn_list,delay,w_list,conn, save=False, name='noname'):
        '''Important: syn_list and w_list must be lists.'''
        for k,conn_ in enumerate(conn):
            for post_id, all_pre in enumerate(conn_):
                for j, pre_id in enumerate(all_pre):
                    for syn,w in zip(syn_list[k],w_list[k]):
                        self.ncl_.append(h.NetCon(ns_src[pre_id], trg.cell[post_id].__dict__[syn].syn, 0, delay[j,post_id], w))

    def make_conn(self, preN, postN, conv):
        conn = np.zeros((postN,conv),dtype=np.int16)
        for i in range(postN):
            conn[i,:]=random.sample(range(preN),conv)
        return conn

    def set_connectivity(self, nsrc, ntrg, conv):
        conn = self.make_conn(nsrc,ntrg,conv)
        conn_ = conn[ np.arange(pc.id(), ntrg, pc.nhost()),:]
        return conn_

    def set_cells_connectivity(self):
        random.seed(self.wseed)
        ''' Set the connectivity '''
        keys = self.weights_neurons_neurons.keys()
        self.connectivity_neurons_neurons = dict.fromkeys(keys)
        for key in keys:
            if key.endswith("ca1"):
                from_po = self.__dict__[key.split("_to_")[0]]
                to_po   = self.__dict__[key.split("_to_")[1]]
                # n = self.nsyns_neurons_neurons[key]
                # self.connectivity_neurons_neurons[key] = self.set_connectivity(from_po.n, to_po.n, n)
                self.connectivity_neurons_neurons[key] = []
                for n in self.nsyns_neurons_neurons[key]:
                    print( from_po.n, to_po.n, n)
                    self.connectivity_neurons_neurons[key].append( self.set_connectivity(from_po.n, to_po.n, n) )

    def f_conns(self):
        self.set_cells_connectivity()
        print("Making connections among populations")
        print("------------------------------------")
        for key in self.weights_neurons_neurons.keys():
            if key.endswith("ca1"):
                from_po = self.__dict__[key.split("_to_")[0]]
                to_po   = self.__dict__[key.split("_to_")[1]]
                delays  = self.delays_neurons_neurons[key]
                syn     = self.syn_neurons_neurons[key]
                weights = self.weights_neurons_neurons[key]
                conn    = self.connectivity_neurons_neurons[key]

                print( key, syn, weights )
                self.set_connections(from_po, to_po, syn, delays, weights, conn)

                # print(key.split("_to_")[0], key.split("_to_")[1], syn, weights)
                # print(conn)
                # print(' ')
                # print(delays)
                # print('-----------------------------------')

    def set_conn_weight(self, conn, weight):
        for nc in conn:
            nc.weight[0] = weight

    def set_connections(self,src,trg,syn_list,delay,w_list,conn, print_=False): # mirar esto con detalle
        #conn = self.make_conn(src.n,trg.ng,conv, trg)
        #conn = self.make_conn(src.n,trg.n,conv, trg)
        #conn = conn[ np.arange(pc.id(), trg.n, pc.nhost()),:]
        for k, conn_ in enumerate(conn):
            for post_id, all_pre in enumerate(conn_):
                #trg_id = trg.pgidlist[post_id]
                #delays = delay[:,post_id]
                #print(post_id,delays)
                for j, pre_id in enumerate(all_pre):
                    src_id = pre_id+src.naux
                    for syn,w in zip(syn_list[k],w_list[k]):
                        #trg.nc.append(pc.gid_connect(src_id, pc.gid2cell(trg_id).__dict__[syn].syn))
                        trg.nc.append(pc.gid_connect(pre_id+src.naux, trg.cell[post_id].__dict__[syn].syn))
                        trg.nc[-1].weight[0] = w
                        trg.nc[-1].delay = delay[k][j,post_id]
                        trg.nc[-1].threshold = 0.0

        if self.SaveConn:
            try:
                print(self.nqcon.size())
            except:
                self.nqcon = h.NQS("id1","id2","w","syn")
                self.nqcon.strdec("syn")
            for post_id, all_pre in enumerate(conn):
                for j, pre_id in enumerate(all_pre):
                    self.nqcon.append(src.cell[pre_id].id,trg.cell[post_id].id,w,syn)

    def calc_lfp(self): # lfp is modeled as a difference between voltages in distal apical and basal compartemnts
        self.vlfp_ca1 = h.Vector(self.pyr_ca1.cell[0].Adend3_volt.size()) #lfp in neuron Vector
        for cell in self.pyr_ca1.cell:
            self.vlfp_ca1.add(cell.Adend3_volt)
            self.vlfp_ca1.sub(cell.Bdend_volt)
        #self.vlfp_ca3.div(len(self.pyr_ca3.cell)) # normalize lfp by amount of pyr cells
        self.vlfp_ca1=np.array(self.vlfp_ca1.to_python())

    def record_synaptic_currents(self): #optional
        ''' save synaptic currents that every neuron receives '''
        for cell in self.cell_ca1: #the three population unified
            cell.record_synapses()

    def get_synaptic_currents(self):
        '''
        It must be called after simulation and record must have been called previously
        '''
        polist = [self.pyr_ca1, self.bas_ca1, self.olm_ca1, self.cck_ca1] # added 11/03/2022
        for po in polist:
            cell = po.cell[0]
            nsecs = len(cell.sec_list)
            seclist = cell.sec_list
            synlist = cell.syn_list

            syn_per_sec = {sec: [] for sec in seclist}
            for sec in seclist:
                for syn in synlist:
                    if sec in syn:
                        syn_per_sec[sec].append(syn)
                        po.__dict__[sec+"_isyn"][syn.replace(sec,"")] = []

            for cell in po.cell:
                for sec, syn_list in syn_per_sec.items():
                    for syn in syn_list:
                        syn_ = syn.replace(sec,"")
                        po.__dict__[sec+"_isyn"][syn_].append( np.array(cell.__dict__[syn+"_i"].to_python()) )
                        if "NMDA" in syn:
                            po.__dict__[sec+"_isyn"][syn_][-1] += np.array( cell.__dict__[syn+"_iNMDA"].to_python() )
                        #if "NMDA" not in syn:
                        #   po.__dict__[sec+"_isyn"][syn_].append( np.array(cell.__dict__[syn+"_i"].to_python()) )
                        #else:
                        #    po.__dict__[sec+"_isyn"][syn_].append( np.array(cell.__dict__[syn+"_i"].to_python() )+ cell.__dict__[syn+"_iNMDA"].to_python())
            for sec, syn_list in syn_per_sec.items():
                #po.__dict__[sec+"_isyn"]["all"] = 0
                for syn in syn_list:
                    syn_ = syn.replace(sec,"")
                    po.__dict__[sec+"_isyn"][syn_] = np.sum( np.array(po.__dict__[sec+"_isyn"][syn_]), axis=0 )
                    # the mean will be computed in host 0 when all elements where placed
                    #po.__dict__[sec+"_isyn"]["all"] += po.__dict__[sec+"_isyn"][syn_]

    def get_synaptic_currents_v2(self):

        #polist = [self.pyr_ca3, self.bas_ca3, self.olm_ca3] # added 16/03/2022
        polist = [self.pyr_ca1, self.bas_ca1, self.olm_ca1]

        for po in polist:
            cell = po.cell[0]
            nsecs = len(cell.sec_list)
            seclist = cell.sec_list
            synlist = cell.syn_list

            for cell in po.cell:
                for syn in synlist:
                    if "NMDA" in syn:
                        po.__dict__[syn+"_mean"].append( np.array( cell.__dict__[syn+"_iNMDA"].to_python()) )
                    else:
                        po.__dict__[syn+"_mean"].append( np.array( cell.__dict__[syn+"_i"].to_python()) )

            for syn in synlist:
                po.__dict__[syn+"_mean"] = np.sum( np.array(po.__dict__[syn+"_mean"]), axis=0 )

    def load_all_CellStims(self,simdur):

        self.ncl_ = []
        print("Loading CellsStims")
        # file_folder = os.path.join(os.getcwd(),"external_inputs")
        file_folder = '/home/jaime/Desktop/hippocampus/external_inputs/'#+self.external_inputs_subfolder
        keys = ['sep_180','sep_360','ec3_180','ec3_360','pyr_ca3']

        self.tvec = dict.fromkeys(keys)
        self.idvec = dict.fromkeys(keys)
        t0 = 0 # 4050 # 50 + thermalization time
        tf = simdur+t0
        extra_label = "_"+self.external_inputs_label
        ncellstims = []
        for key in keys[:-1]:
            file = "external_inputs_"+key+extra_label+".lzma"
            data = file_management.load_lzma(os.path.join(self.inputs_folder,file))
            x,y = data["tvec"], data["idvec"]
            window = np.logical_and(x>=t0,x<=tf)
            self.tvec[key]   = x[window]-t0+125
            self.idvec[key]  = y[window]
            ncellstims.append(np.max(self.idvec[key])+1)

        key=keys[-1]
        extra_label = "_"+self.external_inputs_label +"_"+self.background_noise_label
        file = "external_inputs_"+key+extra_label+".lzma"
        burst_folder = os.path.join(self.inputs_folder, self.burst_level_label)
        data = file_management.load_lzma(os.path.join(burst_folder,file))
        x,y = data["tvec"], data["idvec"]

        window = np.logical_and(x>=t0,x<=tf)
        self.tvec[key]   = x[window]-t0+125
        self.idvec[key]  = y[window]
        ncellstims.append(np.max(self.idvec[key])+1)

        print("Setting the connections")
        self.set_CellStims_connections(ncellstims)

        print("Creating dictionary with all parameters")
        keys1 = self.weights_inputs_neurons.keys()
        keys2 = ["idv", "spikes", "population", "synapses", "delays", "weights", "connectivity"]
        self.external_inputs_data = dict.fromkeys(keys1)
        for key in keys1:
            if key.endswith("ca1"):
                self.external_inputs_data[key] = dict.fromkeys(keys2)
                inputs = key.split("_to_")[0]
                po = self.__dict__[key.split("_to_")[1]]
                self.external_inputs_data[key]["idv"]          = self.idvec[inputs]
                self.external_inputs_data[key]["spikes"]       = self.tvec[inputs]
                self.external_inputs_data[key]["population"]   = po
                self.external_inputs_data[key]["synapses"]     = self.syn_inputs_neurons[key]
                self.external_inputs_data[key]["delays"]       = self.delays_inputs_neurons[key]
                self.external_inputs_data[key]["weights"]      = self.weights_inputs_neurons[key]
                self.external_inputs_data[key]["connectivity"] = self.connectivity_inputs_neurons[key]
