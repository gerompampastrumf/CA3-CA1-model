# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 11:11:21 2021

@author: Asus

All the measurements related with the circuit of hippocampus CA3
(probably will need a further update to consider CA1 too)

Actualizacion de 17/09/2021
voy a anadir una funcion donde obtengo varios plots relaciandos con la nueva version
del modelo de CA3 el cual esta optimizado (hasta que se demuestre lo contrario):
hippocampus_network3.txt
"""

import numpy as np
import matplotlib.pyplot as plt
from neuron import h,gui
import matplotlib.gridspec as gridspec
from spectrum import *
import scipy
from scipy.signal import hilbert, chirp
from collections import Counter
import copy

################################################################################
''' Actualizacion de 17/09/2021 '''
################################################################################
def multitapper(x,fs,NW = 2.5, k=4):
    '''
    x: signal 
    fs: sampling frequency (Hz)
    '''
    nsamples = len(x)
    ft = np.fft.fft( (x-x.mean())/x.std())
    ft = abs(ft[:nsamples//2])**2 / nsamples/fs
    xf = np.fft.fftfreq(nsamples, 1/fs)[:nsamples//2]
    df = np.diff(xf)[0]

    [tapers, eigen] = dpss(nsamples, NW, k)
    Sk_complex, weights, eigenvalues=pmtm( (x-x.mean())/x.std(), e=eigen, v=tapers, NFFT=nsamples, show=False)
    Sk = abs(Sk_complex)**2
    Sk = np.mean(Sk*np.transpose(weights), axis=0)/fs
    Sk = Sk[:nsamples//2]
    return xf, ft, Sk

def get_plots_for_preliminar_analysis(net,t,all_spikes, all_volt, lfp, tmin,tmax,tp,
    figtitle1,figtitle2,figtitle3,figtitle4):
    for po in net.pop_ca3:
        seclist = po.cell[0].sec_list
        for sec in seclist:
            all = np.ones(len(t))*0
            for key in po.__dict__[sec+"_isyn"].keys():
                po.__dict__[sec+"_isyn"][key] = np.sum(po.__dict__[sec+"_isyn"][key], axis=0)/po.n
                all += po.__dict__[sec+"_isyn"][key]
            po.__dict__[sec+"_isyn"]["all"] = all

    lfp = np.sum(lfp, axis=0)/net.pyr_ca3.n

    data_spikes = {}
    for process_data in all_spikes:
        data_spikes.update(process_data)

    ############################################################################
    ''' Spikes and inputs'''
    ############################################################################
    #tmin, tmax = 0.0, simulation_time
    window_plot = np.logical_and(t>=tmin, t<=tmax)

    color, labels = np.ones(net.n).astype(str), np.ones(net.n).astype(str)
    color[:net.n_pop[0]], labels[:net.n_pop[0]] = 'tab:blue', 'pyr'
    color[net.n_pop[0]:net.n_pop[1]], labels[net.n_pop[0]:net.n_pop[1]] = 'tab:orange', 'basket'
    color[net.n_pop[1]:net.n_pop[2]], labels[net.n_pop[1]:net.n_pop[2]] = 'tab:green', 'olm'
    loc = "center left"

    fig = plt.figure(1,figsize=(16,8))
    gs  = gridspec.GridSpec(7,7,height_ratios=[1,0.5,0.5,0.5,0.5,0.1,0.5])
    ax  = [ plt.subplot(gs[i,:6]) for i in range(5)]
    ax.append( plt.subplot(gs[6,:6]))
    for i, spike_times in data_spikes.items():
        ax[0].plot(spike_times, i*np.ones(len(spike_times)), 'o', markersize=.7, color=color[i])
        #ax[0].set_title('Number of hosts: ' + str(int(pc.nhost())))
    ax[0].set_ylabel('# Neuron')
    ax[0].set_xticklabels([ ])
    ax[0].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)

    # Raster external inputs: septum and enthorrinal
    if net.two_septal_inputs:
        if net.ms180_gain > 0:
            ax[1].plot(net.tvec_ms_180,net.idvec_ms_180,'o',color='tab:orange',label='to basket',markersize=1.5)
        if net.ms360_gain > 0:
            ax[1].plot(net.tvec_ms_360,net.idvec_ms_360,'o',color='tab:green', label='to olm',markersize=1.5)
    else:
        if np.logical_and(net.MSGain_basket != 0, net.MSGain_olm !=0):
            label='to basket and olm'
        elif np.logical_and(net.MSGain_basket !=0, net.MSGain_olm == 0):
            label='to basket'
        elif np.logical_and(net.MSGain_basket == 0 , net.MSGain_olm != 0):
            label='to olm'
        ax[1].plot(net.tvec_ms,np.ones(len(net.tvec_ms)),'o',color='black',label=label)

    ax[1].set_ylabel('MS')
    ax[1].set_xticklabels([ ])
    ax[1].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)

    #if net.tvec_dg_burst: #net.tvec_dg_burst.to_python():
    ax[2].plot(net.tvec_dg_burst,net.idvec_dg_burst,'o',color='tab:red',label='>40 Hz',markersize=1.5)
    ax[2].plot(net.tvec_dg_regular,net.idvec_dg_regular,'o',color='tab:purple',label='4 Hz',markersize=1.5)
    ax[2].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)
    ax[2].set_ylabel('DG')
    ax[2].set_xticklabels([ ])

    if net.tvec_ec2_180.to_python():
        ax[3].plot(net.tvec_ec2_180,net.idvec_ec2_180,'o',color='tab:red',label='180',markersize=1.5)
    if net.tvec_ec2_360.to_python():
        ax[3].plot(net.tvec_ec2_360,net.idvec_ec2_360,'o',color='tab:purple',label='360',markersize=1.5)

    ax[3].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)
    ax[3].set_ylabel('EC2')
    ax[3].set_xticklabels([ ])

    ax[4].plot(t[window_plot],lfp[window_plot],color='black')
    ax[4].set_ylabel('LFP [mV]')
    ax[4].set_xlabel('Time [ms]')

    x,y = tappers(t[window_plot],lfp[window_plot])
    ax[5].plot(x,y, color='black')
    ax[5].set_xlabel("Frequency [Hz]")
    ax[5].set_ylabel("Power [a.u.]")
    ax[5].set_xlim([0,80])
    #ax[5].set_yscale('log')
    ax[5].set_ylim([1e-4, np.max(y)*1.05])

    fig.align_labels()
    for i in range(5):
        ax[i].set_xlim([tmin,tmax])
        ax[i].grid(True)
    #plt.show()
    plt.savefig(figtitle1)
    plt.close()
    ############################################################################
    '''synaptic currents plot '''
    ############################################################################
    fig = plt.figure(2,figsize=(16,8))
    gs = gridspec.GridSpec(7,7,height_ratios=[0.5,0.5,0.5,0.5,0.5,0.5,0.5])#, width_ratios=[1,1])
    ax = [ plt.subplot(gs[i,:6]) for i in range(7)]
    index = 0
    color = {"AMPA": "tab:blue", "NMDA": "tab:green", "GABA": "red", "all": "black"}
    po = net.pyr_ca3
    seclist = ["Adend3", "Adend2", "Adend1", "soma", "Bdend"]
    labels = seclist + ["Basket", "OLM"]
    for sec in seclist:
        for key, isyn  in po.__dict__[sec+"_isyn"].items():
            for keyc in color.keys():
                if keyc in key:
                    color_=color[keyc]
            ax[index].plot(t[window_plot], isyn[window_plot], label=key, color=color_)
        ax[index].set_ylabel(labels[index])
        ax[index].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1), fontsize=8)

        index+=1
    for po in net.pop_ca3[1:]:
        seclist = po.cell[0].sec_list
        for sec in seclist:
            for key, isyn  in po.__dict__[sec+"_isyn"].items():
                for keyc in color.keys():
                    if keyc in key:
                        color_=color[keyc]
                        if key == "GABAss":
                            color_="tab:orange"
                ax[index].plot(t[window_plot], isyn[window_plot], label=key, color=color_)
            ax[index].set_ylabel(labels[index])
            ax[index].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)
            index+=1

    for i in range(6):
        ax[i].set_xticklabels([ ])
    for i in range(7):
        ax[i].set_xlim([tmin,tmax])
        ax[i].grid(True)
    ax[-1].set_xlabel('Time [ms]')
    fig.align_labels()
    plt.savefig(figtitle2)
    plt.close()

    ############################################################################
    '''time_phase_spike_distribution'''
    ###########################################################################
    t0=tmin
    binsize=18 #18 da dt = 6.25 ms sacado del paper

    pd = dict.fromkeys(net.cell_labels)
    pd_inputs = dict.fromkeys(net.input_labels)
    pd_sum = dict.fromkeys(net.cell_labels)

    t_ref = net.input_time[-3] #ec3_360
    dt = binsize/360.*net.theta_rythm
    tf = tmax#simulation_time
    tf = tf-np.mod(tf-t0,dt)
    ntime = int((tf-t0)/dt)+1
    time_phase = np.linspace(t0,tf,ntime)

    spiking_array = np.zeros((net.n, len(time_phase)))

    pyr, bas, olm = [],[],[]
    for i, spike_times in data_spikes.items():
        t = np.array(spike_times)
        for tt in t[t>t0]:
            if i < 800:
                pyr.append(tt)
            if np.logical_and(i>=800,i<1000):
                bas.append(tt)
            if np.logical_and(i>=1000,i<1200):
                olm.append(tt)
    # external inputs:
    '''
    def spikes_hist(x,y,ncells,ntime):
        spiking_array = np.zeros((ncells, len(time_phase)))
        x_ = np.round(np.array(x.to_python()),1)
        y_ = np.array(y.to_python()).astype(int)
        for i, spike_time in zip(y_,x_):
            j = int((spike_time-t0)/dt)
            spiking_array[i,j] = 1
        return np.sum(spiking_array,axis=0)
    '''
    bins = time_phase

    pd["pyr_ca3"],x  = np.histogram(pyr, bins=bins)
    pd["bas_ca3"],x  = np.histogram(bas, bins=bins)
    pd["olm_ca3"],x  = np.histogram(olm, bins=bins)

    pd_inputs["sep_180"],x    = np.histogram(net.tvec_ms_180,     bins=bins)#spikes_hist(net.tvec_ms_180, net.idvec_ms_180, 10, ntime)
    pd_inputs["sep_360"],x    = np.histogram(net.tvec_ms_360,     bins=bins)#spikes_hist(net.tvec_ms_360, net.idvec_ms_360, 10, ntime)
    pd_inputs["ec2_180"],x    = np.histogram(net.tvec_ec2_180,    bins=bins)#spikes_hist(net.tvec_ec2_180, net.idvec_ms_180, 100, ntime)
    pd_inputs["ec2_360"],x    = np.histogram(net.tvec_ec2_360,    bins=bins)#spikes_hist(net.tvec_ec2_360, net.idvec_ms_360, 100, ntime)
    pd_inputs["dg_regular"],x = np.histogram(net.tvec_dg_regular, bins=bins)#spikes_hist(net.tvec_dg_regular, net.idvec_dg_regular, 100, ntime)
    pd_inputs["dg_burst"],x   = np.histogram(net.tvec_dg_burst,   bins=bins)#spikes_hist(net.tvec_dg_burst, net.idvec_dg_burst, 100, ntime)
    #pd_inputs["dg_burst"], x  = np.histogram(net.tvec_dg_burst_var, bins=bins)
    x = x[:-1]
    #pd_inputs["dg_burst"] = np.append(pd_inputs["dg_burst"][x<tp],pd_inputs["dg_burst_var"][x>=tp])
    time_cycle = []
    ncycles = []
    cycles = 1
    tx = t_ref
    while tx<=tf:
        if tx>=t0:
            time_cycle.append(tx-t0)
            ncycles.append(cycles)
        tx+=net.theta_rythm
        cycles+=1

    nspikes_cycle_cells  = dict.fromkeys(net.cell_labels)
    nspikes_cycle_inputs = dict.fromkeys(net.input_labels)
    for key in pd.keys():
        nspikes_cycle_cells[key] = []
    for key in pd_inputs.keys():
        nspikes_cycle_inputs[key] = []

    for t1,t2 in zip(time_cycle[:-1],time_cycle[1:]):
        window=np.logical_and(x>=t1,x<t2)
        for key in pd.keys():
            nspikes_cycle_cells[key].append(np.sum(pd[key][window]))
        for key in pd_inputs.keys():
            try:
            #print(np.array(pd_inputs[key])[window])
                nspikes_cycle_inputs[key].append(np.sum(pd_inputs[key][window]))
            except:
                nspikes_cycle_inputs[key].append(0)

    plt.figure(3,figsize=(6,4))
    gs = gridspec.GridSpec(3,7,height_ratios=[0.5,0.5,0.5])#, width_ratios=[1,1])
    ax = [ plt.subplot(gs[i,:6]) for i in range(3)]
    ax[0].step(x,pd["pyr_ca3"],           where='mid', color='tab:blue',   label='pyr')
    ax[1].step(x,pd_inputs["dg_regular"], where='mid', color='tab:red',    label='> 40 Hz')
    ax[1].step(x,pd_inputs["dg_burst"],   where='mid', color='tab:purple', label='4 Hz')
    ax[2].step(x,pd_inputs["ec2_180"],    where='mid', color='tab:red',    label='180')
    ax[2].step(x,pd_inputs["ec2_360"],    where='mid', color='tab:purple', label='360')
    ax[0].fill_between(x,pd["pyr_ca3"],           step="mid", color='tab:blue',   alpha=0.4)
    ax[1].fill_between(x,pd_inputs["dg_regular"], step="mid", color='tab:red',    alpha=0.4)
    ax[1].fill_between(x,pd_inputs["dg_burst"],   step="mid", color='tab:purple', alpha=0.4)
    ax[2].fill_between(x,pd_inputs["ec2_180"],    step="mid", color='tab:red',    alpha=0.4)
    ax[2].fill_between(x,pd_inputs["ec2_360"],    step="mid", color='tab:purple', alpha=0.4)

    ax[1].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1))
    ax[2].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1))

    labels = ['Pyr', 'DG', 'EC2']
    for i in range(3):
        ax[i].set_ylabel(labels[i])
        ax[i].grid(True)
    for i in range(2):
        ax[i].set_xticklabels([ ])
    ax[-1].set_xlabel('Time [ms]')
    plt.savefig(figtitle3)
    plt.close()

    n0 = int((tp-t_ref)/net.theta_rythm)+1
    plt.figure(4,figsize=(6,4))
    plt.plot(ncycles[:-1],nspikes_cycle_cells["pyr_ca3"],'o-',color='tab:blue', label='pyr')
    plt.plot(ncycles[:-1],nspikes_cycle_cells["bas_ca3"],'o-',color='tab:orange', label='bas')
    plt.plot(ncycles[:-1],nspikes_cycle_inputs["dg_burst"],'o-',color='tab:red', label='dg burst')
    plt.plot(ncycles[:-1],nspikes_cycle_inputs["dg_regular"],'o-',color='tab:purple', label='dg regular')
    plt.legend(loc=loc)
    plt.grid(True)
    #plt.ylim([0,400])
    plt.axvline(n0)
    plt.savefig(figtitle4)
    plt.close()

################################################################################
''' Previous fucntions '''
################################################################################

def density(x,y,t0):
    x = np.array(x)
    y = np.array(y)
    timebin = 1
    time_binned = np.linspace(t0+timebin, h.tstop, int((h.tstop-t0)/timebin))
    window = int(timebin/h.dt)
    spk  = np.zeros(len(time_binned))
    for ii in range(len(time_binned)):
        time_binned[ii] = round(time_binned[ii],2)
        spk[ii] = len(np.where(np.logical_and( x>=t0+timebin*ii, x<t0+timebin*(ii+1)))[0])
    time_r = time_binned[window//2:-window//2]
    spk_r = spk[window//2:-window//2]
    return time_r, spk_r

def get_spikes(net):
    size = np.zeros(len(net.cell_labels))
    for i in range(len(net.cell_labels)):
        size[i]=net.__dict__[net.cell_labels[i]].ng
    size = np.cumsum(size)
    size = np.insert(size,0,0)

    spikes = dict.fromkeys(net.cell_labels)
    for i in range(len(net.cell_labels)):
        spk_t, spk_i = [],[]
        population = net.__dict__[net.cell_labels[i]]
        for j in range(population.ng): #range(population.n):
            x = np.array(population.ltimevec[j].to_python())
            y = np.full(len(x),j+size[i])
            spk_t.append(x)
            spk_i.append(y)
        #spk_t, spk_i = np.array(spk_t), np.array(spk_i)
        x,y = np.concatenate(spk_t), np.concatenate(spk_i)
        spikes[net.cell_labels[i]] = [np.sort(x),y[np.argsort(x)]]
    return size,spikes

def get_phase_distribution(net,ncicles=2):
    pd = dict.fromkeys(net.cell_labels)
    size,spikes = get_spikes(net)
    for i,label in enumerate(net.cell_labels):
        try:
            yaux = np.array(spikes[label][0])
            yaux = yaux[yaux>=400]
            yp = np.mod( (yaux-net.time_initial-net.relative_time)/net.theta_rythm*360,360) #two theta cicles
            a,b = np.histogram(yp,bins=20,density=True,range=(0,360))
            b = b[1:]-np.diff(b)[0]
            a,b = np.append(a,a[0]),np.append(b,360)
            b1,a1 = b[b>=270]-360, a[b>=270]
            b2,a2 = b[b<=90]+360,  a[b<=90]
            b = b1.tolist()+b.tolist()+b2.tolist()
            a = a1.tolist()+a.tolist()+a2.tolist()
            pd[label] = [np.array(b),np.array(a)]
        except:
            pass

    if net.EC3_inputs:
        yaux = np.sort(net.EC3_inputs)
        yaux = yaux[yaux>=200]
        yp = np.mod( (yaux-net.time_initial-net.relative_time)/net.theta_rythm*360,360) #two theta cicles
        a,b = np.histogram(yp,bins=20,density=True,range=(0,360))
        b = b[1:]-np.diff(b)[0]
        a,b = np.append(a,a[0]),np.append(b,360)

        b1,a1 = b[b>=270]-360, a[b>=270]
        b2,a2 = b[b<=90]+360,  a[b<=90]
        b = b1.tolist()+b.tolist()+b2.tolist()
        a = a1.tolist()+a.tolist()+a2.tolist()
        pd["EC3_inputs_ca1"] = [np.array(b),np.array(a)]

    if net.EC2_inputs:
        yaux = np.sort(net.EC2_inputs)
        yaux = yaux[yaux>=200]
        yp = np.mod( (yaux-net.time_initial-net.relative_time)/net.theta_rythm*360,360) #two theta cicles
        a,b = np.histogram(yp,bins=20,density=True,range=(0,360))
        b = b[1:]-np.diff(b)[0]
        a,b = np.append(a,a[0]),np.append(b,360)
        b1,a1 = b[b>=270]-360, a[b>=270]
        b2,a2 = b[b<=90]+360,  a[b<=90]
        b = b1.tolist()+b.tolist()+b2.tolist()
        a = a1.tolist()+a.tolist()+a2.tolist()
        pd["EC2_inputs_ca3"] = [np.array(b),np.array(a)]

    if net.DG_inputs:
        yaux = np.sort(net.DG_inputs)
        yaux = yaux[yaux>=200]
        yp = np.mod( (yaux-net.time_initial-net.relative_time)/net.theta_rythm*360,360) #two theta cicles
        a,b = np.histogram(yp,bins=20,density=True,range=(0,360))
        b = b[1:]-np.diff(b)[0]
        a,b = np.append(a,a[0]),np.append(b,360)
        b1,a1 = b[b>=270]-360, a[b>=270]
        b2,a2 = b[b<=90]+360,  a[b<=90]
        b = b1.tolist()+b.tolist()+b2.tolist()
        a = a1.tolist()+a.tolist()+a2.tolist()
        pd["DG_inputs_ca3"] = [np.array(b),np.array(a)]

    net.phase_distributions = pd

def plots_new(net,ti,tf):
    ca3_cell_labels = []
    ca1_cell_labels = []
    for element in net.cell_labels:
        if element.endswith("ca3"):
            ca3_cell_labels.append(element)
        if element.endswith("ca1"):
            ca1_cell_labels.append(element)

    color = ['#1f77b4', '#ff7f0e', '#2ca02c','#9467bd']
    net.calc_lfp()
    freq_spectrum,power_spectrum = tappering_analysis(net)[:2]
    ## a, b, c not going to be used here
    size,spikes = get_spikes(net)
    ymin,ymax = size[0],size[3]

    xt = np.linspace(0,h.tstop,int(h.tstop/h.dt)+1)
    plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(5,1,height_ratios=[1,0.5,0.5,0.1,0.5])#, width_ratios=[1,1])

    ax = []
    ax.append(plt.subplot(gs[0,0]))
    ax.append(plt.subplot(gs[1,0]))
    ax.append(plt.subplot(gs[2,0]))
    ax.append(plt.subplot(gs[4,0]))

    #ax.append(plt.subplot(gs[0,1]))
    #ax.append(plt.subplot(gs[1,1]))
    #ax.append(plt.subplot(gs[2,1]))
    #ax.append(plt.subplot(gs[4,1]))

    labels = [ca3_cell_labels, ca1_cell_labels]
    inputs = [ [net.DG_inputs,net.EC2_inputs], net.EC3_inputs]
    inputs_label =[ [ net.DG_inputs_label, net.EC2_inputs_label], net.EC3_inputs_label]

    iax = 0
    for k in range(1): #CA3 and CA1 solo ca3
        # Raster plot
        for i,cell_type in enumerate(labels[k]):
            x,y = spikes[cell_type]
            ax[iax].plot(x,y,'o',color=color[i],markersize=.7)

        ax[iax].set_ylabel('# Neuron',fontsize=15)
        ax[iax].grid(True)
        ax[iax].set_xticklabels([ ])
        ax[iax].set_xlim([ti,tf])
        iax+=1

        # Raster external inputs: septum and enthorrinal
        n_inputs = len(net.MS_inputs)
        xms = np.array(net.MS_inputs.to_python())
        yms = np.repeat(50,n_inputs)
        ax[iax].plot(xms,yms,'o',markersize=10,color='black')
        ax[iax].set_xticklabels([ ])
        ax[iax].grid(True)
        ax[iax].set_ylabel('Inputs',fontsize=15)

        if inputs[k]:
            if len(inputs[k]) == 1:
                n_inputs = len(inputs[k])
                xms = np.array(inputs[k])
                yms = np.array(inputs_label[k])
                ax[iax].plot(xms[xms>=ti],yms[xms>=ti],'o',markersize=.7,color='#d62728')
                ax[iax].set_ylim([0,np.max(yms)])
            else:
                y0 = 0
                for kk in range(len(inputs[k])):
                    if inputs[k][kk]:
                        n_inputs = len(inputs[k][kk])
                        xms = np.array(inputs[k][kk])
                        yms = y0+np.array(inputs_label[k][kk])
                        ax[iax].plot(xms[xms>=ti],yms[xms>=ti],'o',markersize=.7,color=['#d62728','#9467bd'][kk])
                        y0 = y0+np.max(yms)
                        ax[iax].set_ylim([0,np.max(yms)])
        ax[iax].set_xlim(ti,tf)
        iax+=1

        y = net.vlfp_ca3#,net.vlfp_ca1][k]
        ax[iax].plot(xt,y,color='black')
        ax[iax].set_ylabel('LFP [mV]',fontsize=15)
        ax[iax].set_xlabel('time [ms]',fontsize=15)
        ax[iax].grid(True)
        ax[iax].set_xlim([ti,tf])
        iax+=1

        ax[iax].plot(freq_spectrum[k],power_spectrum[k])
        ax[iax].set_xlim([0,120])
        ax[iax].set_xlabel('Frequency [Hz]',fontsize=15)
        ax[iax].set_ylabel('Power [a.u.]', fontsize=15)
        ax[iax].grid(True)
        iax+=1

    #ax[4].set_xlim([xmin,xmax])
    #ax[5].set_xlim([xmin,xmax])
    #ax[6].set_xlim([xmin,xmax])

def firing_rate_old(net,m=8,sg=20,t0=1000,smooth=False,inputs_as_sin=False):
    spikes = get_spikes(net)[1]
    timebin, runtime = 0.5, h.tstop
    time_binned = np.linspace(t0+timebin, runtime, int((runtime-t0)/timebin))
    n = len(time_binned)
    window = m*int(sg/timebin)

    rate = dict.fromkeys(net.cell_labels)
    spk  = dict.fromkeys(net.cell_labels)

    time_r = time_binned[window//2:-window//2]
    ttb = np.arange(-window/2,window/2+1)
    for label in net.cell_labels:
        rate[label] = np.zeros(len(time_binned))
        spk[label]  = np.zeros(len(time_binned))
        x = spikes[label][0]
        nexc = net.__dict__[label].n
        for ii in range(len(time_binned)):
            spk[label][ii] = len(np.where(np.logical_and(x >= t0+timebin*ii, x < t0+timebin*(ii+1)))[0])
        gaussian = np.exp(-ttb**2/(2.0*(sg/timebin)**2))
        gaussian /= np.sum(gaussian)
        for ii in range(window//2, n-window//2):
            ss  = spk[label][ii-window//2:ii+window//2+1]
            rate[label][ii] = np.sum(gaussian*ss)

        rate[label] = rate[label][window//2:-window//2]
        rate[label] = rate[label]/timebin/nexc*1e3
        rate[label] -= rate[label].mean()
        if rate[label].std()>1e-8: #in case one of the population does not fire because of the parameter scan
            rate[label] /= rate[label].std()

    arg = np.arcsin(1)-2.0*np.pi*net.time_initial/net.theta_rythm
    rate["ms"] = np.sin(arg+2.0*np.pi/net.theta_rythm*time_binned)
    rate["ms"] = rate["ms"][window//2:-window//2]
    rate["ms"] /= rate["ms"].std()

    if inputs_as_sin:
        arg = np.arcsin(1)-2.0*np.pi*(net.time_initial+net.relative_time)/net.theta_rythm
        rate["ec2"] = np.sin(arg+2.0*np.pi/net.theta_rythm*time_binned)
        rate["ec2"] = rate["ec2"][window//2:-window//2]
        rate["ec2"] /= rate["ec2"].std()

        arg = np.arcsin(1)-2.0*np.pi*(net.time_initial+net.relative_time-15.63+net.theta_rythm)/net.theta_rythm
        rate["dg"] = np.sin(arg+2.0*np.pi/net.theta_rythm*time_binned)
        rate["dg"] = rate["dg"][window//2:-window//2]
        rate["dg"] /= rate["dg"].std()
    else:
        rate["ec2"] = np.zeros(len(time_binned))
        spk["ec2"]  = np.zeros(len(time_binned))
        x = np.array(net.EC2_inputs)
        nexc = 800
        for ii in range(len(time_binned)):
            spk["ec2"][ii] = len(np.where(np.logical_and(x >= t0+timebin*ii, x < t0+timebin*(ii+1)))[0])
        gaussian = np.exp(-ttb**2/(2.0*(sg/timebin)**2))
        gaussian /= np.sum(gaussian)
        for ii in range(window//2, n-window//2):
            ss  = spk["ec2"][ii-window//2:ii+window//2+1]
            rate["ec2"][ii] = np.sum(gaussian*ss)

        rate["ec2"] = rate["ec2"][window//2:-window//2]
        rate["ec2"] = rate["ec2"]/timebin/nexc*1e3
        rate["ec2"] -= rate["ec2"].mean()
        rate["ec2"] /= rate["ec2"].std()

        rate["dg"] = np.zeros(len(time_binned))
        spk["dg"]  = np.zeros(len(time_binned))
        x = np.array(net.DG_inputs)
        nexc = 1600
        for ii in range(len(time_binned)):
            spk["dg"][ii] = len(np.where(np.logical_and(x >= t0+timebin*ii, x < t0+timebin*(ii+1)))[0])
        gaussian = np.exp(-ttb**2/(2.0*(sg/timebin)**2))
        gaussian /= np.sum(gaussian)
        for ii in range(window//2, n-window//2):
            ss  = spk["dg"][ii-window//2:ii+window//2+1]
            rate["dg"][ii] = np.sum(gaussian*ss)

        rate["dg"]  = rate["dg"][window//2:-window//2]
        rate["dg"]  = rate["dg"]/timebin/nexc*1e3
        rate["dg"] -= rate["dg"].mean()
        rate["dg"] /= rate["dg"].std()

    time = time_r
    rates = rate
    if smooth:
        time, rates = smoothing(time_r,rate)

    return time, rates

def smoothing(x,y,m=150):
    rates = y
    for label in y.keys():
        n = len(rates[label])
        xaux = np.zeros(n-m)
        for j in range(m//2,n-m//2):
            xaux[j-m//2] = np.mean( rates[label][j-m//2:j+m//2] )
        rates[label] = xaux
    time = x[m//2:n-m//2]
    return time,rates

def maxima(x,y):
    from scipy.signal import find_peaks
    n = len(y)
    # smoothing
    m = 30
    xaux = np.zeros(n-m)
    for j in range(m//2,n-m//2):
        xaux[j-m//2] = np.mean( y[j-m//2:j+m//2] )
    #peaks, _ = find_peaks(xaux, height = np.mean(xaux),distance = 100)# int(0.5*(14.0/0.5)))
    peaks, _ = find_peaks(xaux,distance = 100)# int(0.5*(14.0/0.5)))

    return m//2 + peaks.astype(int)

def phi(net):
    peaks = dict.fromkeys(net.cell_labels)
    phase = dict.fromkeys(net.cell_labels)
    imin  = dict.fromkeys(net.cell_labels)
    imax  = dict.fromkeys(net.cell_labels)

    time, rates = net.firing_rate(net,smooth=True)
    n = len(time)
    for j,label in enumerate(net.cell_labels):
        x,y = time, rates[label]
        peaks[label] = net.maxima(x,y)
        phase[label] = np.zeros(n)
        for i in range(len(peaks[label])-1):
            a, b = peaks[label][i], peaks[label][i+1]
            bins = b-a+1
            theta = np.linspace(0,2.0*np.pi,bins)
            phase[label][ peaks[label][i] : peaks[label][i+1]+1 ] = theta

        imin[label] = peaks[label][0]
        imax[label] = peaks[label][-1]

    peaks_ms = net.maxima(time,rates["ms"])
    phase_ms = np.zeros(n)
    for i in range(len(peaks_ms)-1):
        a,b = peaks_ms[i], peaks_ms[i+1]
        bins = b-a+1
        theta = np.linspace(0,2.0*np.pi,bins)
        phase_ms[peaks_ms[i]:peaks_ms[i+1]+1] = theta

    phase["ms"] = phase_ms
    imin["ms"] = peaks_ms[0]
    imax["ms"] = peaks_ms[-1]

    peaks_ec2 = net.maxima(time,rates["ec2"])
    phase_ec2 = np.zeros(n)
    for i in range(len(peaks_ec2)-1):
        a,b = peaks_ec2[i],peaks_ec2[i+1]
        bins = b-a+1
        theta = np.linspace(0,2.0*np.pi,bins)
        phase_ec2[peaks_ec2[i]:peaks_ec2[i+1]+1] = theta

    phase["ec2"] = phase_ec2
    imin["ec2"]  = peaks_ec2[0]
    imax["ec2"]  = peaks_ec2[-1]
    return imin, imax, phase

def phase_diff(net):
    imin, imax, phase = net.phi(net)
    dphase = dict.fromkeys(net.cell_labels)
    ms_min, ms_max = imin["ec2"], imax["ec2"] #changed from ms to ec2
    phase_ms = phase["ec2"]
    for label in net.cell_labels:
        imin_ = np.max([imin[label],ms_min])
        imax_ = np.min([imax[label],ms_max])
        dphase[label] = np.mod((phase_ms-phase[label]),2.0*np.pi)[imin_:imax_]
    return dphase

def tappering_analysis(net,t0=1000,NW = 2.5, k=4):
    lowcut  = [30.0,12.0, 3.0]
    highcut = [80.0,30.0, 12.0]

    net.calc_lfp()
    Power_theta = []
    Power_beta  = []
    Power_gamma = []
    theta_max, theta_max_amp = [],[]
    beta_max,  beta_max_amp  = [],[]
    gamma_max, gamma_max_amp = [],[]

    freq_spectrum  = []
    Power_spectrum = []

    time = np.linspace(0,h.tstop, int(h.tstop/h.dt)+1)
    for y in [net.vlfp_ca3[time>=t0]]:#, net.vlfp_ca1][:1]: #solo ca3
        nsamples = len(y)
        # classical FFT
        yf = np.fft.fft(y-y.mean())
        yf = abs(yf)**2 / nsamples * h.dt
        xf = np.linspace(0.0, 1.0/(2.0*h.dt*1e-3), nsamples//2)
        df = np.diff(xf)[0]

        #NW=2.5 #k=4
        [tapers, eigen] = dpss(nsamples, NW, k)
        Sk_complex, weights, eigenvalues=pmtm(y-y.mean(), e=eigen, v=tapers, NFFT=nsamples, show=False)
        Sk = abs(Sk_complex)**2
        Sk = np.mean(Sk*np.transpose(weights), axis=0)*h.dt

        # comparing total power
        #norm1,norm2 = np.sum(abs(yf[0:nsamples//2])), np.sum(Sk[0:nsamples//2])
        theta_band = np.logical_and(xf>=lowcut[2],xf<=highcut[2])
        beta_band  = np.logical_and(xf>=lowcut[1],xf<=highcut[1])
        gamma_band = np.logical_and(xf>=lowcut[0],xf<=highcut[0])
        Sk = Sk[0:nsamples//2]

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

    return freq_spectrum, Power_spectrum, power, fmax, amax

def synchro_index(net): #synchronization index
    si = dict.fromkeys(net.cell_labels)
    for i,label in enumerate(net.cell_labels):
        n = net.__dict__[label].n
        sgv, sgvi = 0, np.zeros(n)
        for j in range(n):
            sgv+=np.array(net.__dict__[label].cell[j].soma_volt.to_python())
            sgvi[j] = np.array(net.__dict__[label].cell[j].soma_volt.to_python()).std()
        sgv /= n
        sgv = sgv.std()
        si[label] = n*sgv/np.sum(sgvi)
    return si

def activity_index(net,t0=1000):
    spikes = get_spikes(net)[1]
    timebin, runtime = 5, h.tstop
    time_binned = np.linspace(t0+timebin, runtime, int((runtime-t0)/timebin))
    n = len(time_binned)

    tc = time_binned[0]
    spk  = dict.fromkeys(net.cell_labels)
    act  = dict.fromkeys(net.cell_labels)
    act_std = dict.fromkeys(net.cell_labels)
    par  = dict.fromkeys(net.cell_labels)
    for label in net.cell_labels:
        spk[label]  = np.zeros(len(time_binned))
        x = spikes[label][0]
        nexc = net.__dict__[label].n
        for ii in range(len(time_binned)):
            spk[label][ii] = len(np.where(np.logical_and(x >= t0+timebin*ii, x < t0+timebin*(ii+1)))[0])

        act[label]=[]
        act_std[label]=[]
        tc = time_binned[0]
        while tc < runtime:
            tw = np.logical_and(time_binned>=tc, time_binned<tc+net.theta_rythm)
            act[label].append(np.max(spk[label][tw])/net.__dict__[label].n)
            tc+=net.theta_rythm
        act[label] = np.array(act[label])
        act_std[label]= np.std(act[label])
        #par[label] = np.log(np.prod(act[label]/np.max(act[label])))
        par[label] = np.log(np.prod(act[label]))
        if par[label]==float('-inf'):
            par[label]=1
        else:
            par[label]-= len(act[label])*np.log(np.max(act[label]))
        act[label] = np.mean(act[label])

    return act, act_std, par

################################################################################
''' eddited on 12th may '''
################################################################################
def phasediff_corr(net):
    time, rates, input_rate = firing_rate(net, t0=1000,smooth=True)
    n = len(time)
    tlagmax = (n-1)*0.5
    tlag = np.linspace(-tlagmax, tlagmax, 2*n-1)
    phase = dict.fromkeys(net.cell_labels)
    window = np.logical_and(tlag>=-net.theta_rythm/2., tlag<net.theta_rythm/2.)
    for label in net.cell_labels:
        cc = np.correlate(rates[label],input_rate["ec3_360"],'full')/net.__dict__[label].n
        phase[label] = tlag[window][np.argmax(cc[window])] #anadido este simbolo, porque de alguna manera estamos sacando el antisimetrico
        phase[label] = np.mod(phase[label]/net.theta_rythm,1)

    return phase

################################################################################
''' added on 12th may '''
################################################################################
def phasediff_hilbert(net):
    #esta parte igual es excesiva pero puede ser util tenerla ya para el futuro
    time,rates, input_rate = firing_rate(net)
    phase         = dict.fromkeys(net.cell_labels)
    dphase        = dict.fromkeys(net.cell_labels)
    dphase_median = dict.fromkeys(net.cell_labels)
    dphase_mean   = dict.fromkeys(net.cell_labels)
    dphase_std    = dict.fromkeys(net.cell_labels)

    phase0 = hilbert_phase(input_rate["ec3_360"])/(2*np.pi)
    #plt.plot(phase0)
    #plt.plot(rates["pyr_ca3"])
    for label in net.cell_labels:
        phase[label] = hilbert_phase(rates[label])/(2*np.pi)
        dphase[label] = np.mod(phase0-phase[label],1)
    #plt.plot(phase["bas_ca3"])
    #plt.plot(rates["bas_ca3"])

    for label in net.cell_labels:
        phi0 = np.median(dphase[label]) #important to center the dphases ditributions
        w1 = np.where(dphase[label]<phi0-0.5)
        w2 = np.where(dphase[label]>=phi0+0.5)
        dphase[label][w1]+=1
        dphase[label][w2]-=1

        #plt.plot(dphase[label])
        dphase_median[label] = np.mod(np.median(dphase[label]),1)
        dphase_mean[label]   = np.mod(np.mean(dphase[label]),1)
        dphase_std[label]    = np.std(dphase[label])
        #plt.axhline(dphase_median[label])

    return dphase_mean, dphase_median, dphase_std

def hilbert_phase(signal):
    signal-=signal.mean()
    #signal/=signal.std()
    analytic_signal = hilbert(signal)
    instantaneous_phase = np.angle(analytic_signal)
    return np.pi+instantaneous_phase

def firing_rate(net,m=8,sg=20,t0=1000,smooth=False,inputs_as_sin=False):
    spikes = get_spikes(net)[1]
    timebin, runtime = 0.5, h.tstop
    time_binned = np.linspace(t0+timebin, runtime, int((runtime-t0)/timebin))
    n = len(time_binned)
    window = m*int(sg/timebin)

    rate = dict.fromkeys(net.cell_labels)
    spk  = dict.fromkeys(net.cell_labels)

    time_r = time_binned[window//2:-window//2]
    ttb = np.arange(-window/2,window/2+1)
    for label in net.cell_labels:
        rate[label] = np.zeros(len(time_binned))
        spk[label]  = np.zeros(len(time_binned))
        x = spikes[label][0]
        nexc = net.__dict__[label].n
        for ii in range(len(time_binned)):
            spk[label][ii] = len(np.where(np.logical_and(x >= t0+timebin*ii, x < t0+timebin*(ii+1)))[0])
        gaussian = np.exp(-ttb**2/(2.0*(sg/timebin)**2))
        gaussian /= np.sum(gaussian)
        for ii in range(window//2, n-window//2):
            ss  = spk[label][ii-window//2:ii+window//2+1]
            rate[label][ii] = np.sum(gaussian*ss)

        rate[label] = rate[label][window//2:-window//2]
        rate[label] = rate[label]/timebin/nexc*1e3
        rate[label] -= rate[label].mean()
        if rate[label].std()>1e-8: #in case one of the population does not fire because of the parameter scan
            rate[label] /= rate[label].std()

        input_rate = dict.fromkeys(net.input_labels)
        t1,t2 = time_r[0],time_r[-1]
        for i in range(len(net.input_time)):
            t0,tf = net.input_time[i][0],h.tstop
            dt = timebin
            t = np.linspace(t0,tf, int((tf-t0)/dt)+1)
            w = 2.0*np.pi/net.theta_rythm
            input_rate[net.input_labels[i]] = np.cos(w*(t-t0))[np.logical_and(t>=t1,t<=t2)]

    return time_r, rate, input_rate

################################################################################
''' added on 13th may '''
################################################################################
def time_phase_spike_distribution(net,t0=1000,binsize=18): #18 da dt = 6.25 ms
    ''' phase spike distribution of cells and artificial cells.
        I obatain here the mean, median and deviation in each cycle '''
    pd = dict.fromkeys(net.cell_labels)
    pd_inputs = dict.fromkeys(net.input_labels)
    pd_sum = dict.fromkeys(net.cell_labels)

    t_ref = net.input_time[-3] #ec3_360
    tf = h.tstop
    dt = binsize/360.*net.theta_rythm
    time_phase = np.linspace(t0,tf,int((tf-t0)/dt)+1)

    x,y = np.array(net.tvec.to_python()), np.array(net.idvec.to_python())
    for i,label in enumerate(net.input_labels):
        selectioned = np.logical_and(y>=net.nsl_n[i],y<net.nsl_n[i+1])
        x_aux, y_aux = x[selectioned], y[selectioned]

        counts = []
        for j in range(len(time_phase)-1):
            window = np.logical_and(x_aux>=time_phase[j], x_aux<time_phase[j+1])
            counts.append(len(y_aux[window]))
        pd_inputs[label] = np.array(counts).astype(int)

    size, spikes = get_spikes(net)
    for i, label in enumerate(net.cell_labels):
        x_aux, y_aux = spikes[label][0],spikes[label][1]
        counts = []
        for j in range(len(time_phase)-1):
            window = np.logical_and(x_aux>=time_phase[j], x_aux<time_phase[j+1])
            counts.append(len(y_aux[window]))
        pd[label] = np.array(counts).astype(int)

    time_phase -= t0
    time_phase -= t_ref #fijamos el 0 en el 0 grados del ciclo de EC3-360
    time_phase = time_phase/125. #normalizado, una vuelta es 0(0)-1(2pi)
    time_phase = time_phase[:-1] + binsize/360/2
    ############################################################################
    ''' Overlapp of the densities '''
    ############################################################################
    ncicles = int(time_phase[-1])
    window = np.logical_and(time_phase>=0, time_phase<1)
    cycle_phase = time_phase[window]
    two_cycle_phase = np.array(list(cycle_phase) + list(1+cycle_phase))

    phase_sum_mean = dict.fromkeys(net.cell_labels)
    phase_sum_median = dict.fromkeys(net.cell_labels)
    phase_sum_std    = dict.fromkeys(net.cell_labels)

    #plt.figure()
    for label in net.cell_labels:
        counts = 0
        phases = []
        for i in range(ncicles):
            window = np.logical_and(time_phase>=i, time_phase<i+1)
            counts += pd[label][window]
            for j, count in enumerate(pd[label][window]):
                phases.append(np.ones(count)*cycle_phase[j])
        phases=np.concatenate(phases)
        mdn = np.median(phases)
        print(mdn)
        w1 = np.where(phases>=mdn+0.5)[0]
        w2 = np.where(phases<mdn-0.5)[0]
        phases[w1] -= 1
        phases[w2] += 1
        pd_sum[label] = counts

        phase_sum_mean[label] = np.mod(np.mean(phases),1)
        phase_sum_median[label] = np.mod(np.median(phases),1)
        phase_sum_std[label] = np.std(phases)

        pd_sum[label] = np.array(list(pd_sum[label])+list(pd_sum[label]))
        #plt.step(two_cycle_phase,pd_sum[label], where='mid')
        #plt.fill_between(two_cycle_phase,pd_sum[label],step="mid", alpha=0.4)
        #plt.axvline(phase_sum_median[label],color='black')

    ############################################################################
    ''' Dynamic measurements '''
    ############################################################################
    phase_m_mean = dict.fromkeys(net.cell_labels)
    phase_m_median = dict.fromkeys(net.cell_labels)
    phase_m_std    = dict.fromkeys(net.cell_labels)
    phase_s_mean = dict.fromkeys(net.cell_labels)
    phase_s_median = dict.fromkeys(net.cell_labels)
    phase_s_std    = dict.fromkeys(net.cell_labels)

    for label in net.cell_labels:
        mean_, median_, std_ = [],[],[]
        for i in range(ncicles):
            window = np.logical_and(time_phase>= phase_sum_median[label]-0.5+i, time_phase<phase_sum_median[label]+i+0.5)
            window_phases = []
            for j,count in enumerate(pd[label][window]):
                window_phases.append( np.ones(count)*time_phase[window][j])

            window_phases = np.concatenate(window_phases)
            #median_.append( np.mod(np.median(window_phases),1 ))
            #mean_.append( np.mod(np.mean(window_phases),1) )
            #std_.append( np.std(window_phases))
            median_.append( np.median(window_phases))
            mean_.append (np.mean(window_phases) )
            std_.append( np.std(window_phases))

        phase_m_mean[label]   = np.mean(np.mod(mean_,1))
        phase_m_median[label] = np.mean(np.mod(median_,1))
        phase_m_std[label]    = np.mean(std_)
        phase_s_mean[label]   = np.std(mean_)
        phase_s_median[label] = np.std(median_)
        phase_s_std[label]    = np.std(std_)

        plt.figure(figsize=(6,3))
        gs = gridspec.GridSpec(4,1,height_ratios=[0.5,0.5,0.5,0.5])#, width_ratios=[1,1])
        ax = [ plt.subplot(gs[i]) for i in range(4)]
        ax[0].step(time_phase,pd_inputs["sep_180"], where='mid')
        ax[0].step(time_phase,pd_inputs["sep_360"], where='mid')
        ax[1].step(time_phase,pd_inputs["ec2_180"], where='mid')
        ax[1].step(time_phase,pd_inputs["ec2_360"], where='mid')
        ax[2].step(time_phase,pd_inputs["ec3_180"], where='mid')
        ax[2].step(time_phase,pd_inputs["ec3_360"], where='mid')
        ax[3].step(time_phase,pd_inputs["dg_regular"], where='mid')
        ax[3].step(time_phase,pd_inputs["dg_burst"], where='mid')#return time_phase, pd, pd_inputs
        ax[0].fill_between(time_phase,pd_inputs["sep_180"],step="mid", alpha=0.4)
        ax[0].fill_between(time_phase,pd_inputs["sep_360"],step="mid", alpha=0.4)
        ax[1].fill_between(time_phase,pd_inputs["ec2_180"],step="mid", alpha=0.4)
        ax[1].fill_between(time_phase,pd_inputs["ec2_360"],step="mid", alpha=0.4)
        ax[2].fill_between(time_phase,pd_inputs["ec3_180"],step="mid", alpha=0.4)
        ax[2].fill_between(time_phase,pd_inputs["ec3_360"],step="mid", alpha=0.4)
        ax[3].fill_between(time_phase,pd_inputs["dg_regular"],step="mid", alpha=0.4)
        ax[3].fill_between(time_phase,pd_inputs["dg_burst"],step="mid", alpha=0.4)#return time_phase, pd, pd_inputs

        plt.figure(figsize=(8,6))
        gs = gridspec.GridSpec(3,1,height_ratios=[0.5,0.5,0.5])#, width_ratios=[1,1])
        ax = [ plt.subplot(gs[i]) for i in range(3)]
        ax[0].step(time_phase,pd["pyr_ca3"], where='mid')
        ax[0].fill_between(time_phase,pd["pyr_ca3"], step="pre", alpha=0.4)

        ax[1].step(time_phase,pd["bas_ca3"], where='mid')
        ax[1].fill_between(time_phase,pd["bas_ca3"], step="pre", alpha=0.4)
        ax[2].step(time_phase,pd["olm_ca3"], where='mid')
        ax[2].fill_between(time_phase,pd["olm_ca3"], step="pre", alpha=0.4)
        plt.show()
    #return phase_sum_mean, phase_sum_median, phase_sum_std, phase_m_mean, phase_m_median, phase_m_std, phase_s_mean, phase_s_median, phase_s_std

def get_phases_all_methods(net):
    a1 = phasediff_corr(net)
    b1,b2,b3 = phasediff_hilbert(net)
    c1,c2,c3,d1,d2,d3,d4,d5,d6 = time_phase_spike_distribution(net)

    data = np.zeros((len(net.cell_labels),13))
    for i,l in enumerate(net.cell_labels):
        data[i] = np.array([a1[l],b1[l],b2[l],b3[l],c1[l],c2[l],c3[l],d1[l],d2[l],d3[l],d4[l],d5[l],d6[l]])

    return np.round(data,4)


'''
# Added on 02/03/2022 (provisional)
# last update 11/03/2022
'''

def raster_plot(net,spikes,tmin,tmax,figtitle):

    ax_list = [['(a)'], ['(a)'], ['(b)'], ['(b)'], ['(c)'], ['(d)'], ['(e)'], ['(f)']]
    fig, ax = plt.subplot_mosaic(ax_list, constrained_layout=True, figsize=(16,8))

    '''
    color, labels = np.ones(net.n).astype(str), np.ones(net.n).astype(str)
    color[:net.n_pop[0]], labels[:net.n_pop[0]] = 'tab:blue', 'pyr'
    color[net.n_pop[0]:net.n_pop[1]], labels[net.n_pop[0]:net.n_pop[1]] = 'tab:orange', 'basket'
    color[net.n_pop[1]:net.n_pop[2]], labels[net.n_pop[1]:net.n_pop[2]] = 'tab:green', 'olm'
    color[net.n_pop[2]:net.n_pop[3]], labels[net.n_pop[2]:net.n_pop[3]] = 'tab:blue', 'pyr'
    color[net.n_pop[3]:net.n_pop[4]], labels[net.n_pop[3]:net.n_pop[4]] = 'tab:orange', 'basket'
    color[net.n_pop[4]:net.n_pop[5]], labels[net.n_pop[4]:net.n_pop[5]] = 'tab:green', 'olm'
    '''

    color = ["tab:blue","tab:orange","tab:green"]
    label = ["pyr", "basket", "olm"]
    loc = "center left"
    data_spikes = dict.fromkeys(spikes.keys())
    for region in data_spikes.keys():
        data_spikes[region] = {}
        for process_data in spikes[region]:
            data_spikes[region].update(process_data)

    if net.create_ca1_network:
        index, spikes = [],[]
        for i, spike_times in data_spikes["ca1"].items():
            index.append( i*np.ones(len(spike_times)) )
            spikes.append( spike_times)
        index=np.concatenate(index).astype(int)
        spikes=np.concatenate(spikes)
        int_pyr = np.logical_and(index>=net.n_pop[2], index<net.n_pop[3])
        int_bas = np.logical_and(index>=net.n_pop[3], index<net.n_pop[4])
        int_olm = np.logical_and(index>=net.n_pop[4], index<net.n_pop[5])
        for i, int_ in enumerate([int_pyr,int_bas,int_olm]):
            ax['(a)'].plot(spikes[int_], index[int_], 'o', markersize=.7, color=color[i],label=label[i])
        ax['(a)'].set_ylabel('# Neurons CA1')
        ax['(a)'].set_xticklabels([ ])
        ax['(a)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)

    index, spikes = [],[]
    for i, spike_times in data_spikes["ca3"].items():
        index.append( i*np.ones(len(spike_times)) )
        spikes.append( spike_times)
    index=np.concatenate(index).astype(int)
    spikes=np.concatenate(spikes)
    int_pyr = np.logical_and(index>=0, index<net.n_pop[0])
    int_bas = np.logical_and(index>=net.n_pop[0], index<net.n_pop[1])
    int_olm = np.logical_and(index>=net.n_pop[1], index<net.n_pop[2])
    for i, int_ in enumerate([int_pyr,int_bas,int_olm]):
        ax['(b)'].plot(spikes[int_], index[int_], 'o' , markersize=.7, color=color[i],label=label[i])
    ax['(b)'].set_ylabel('# Neurons CA3')
    ax['(b)'].set_xticklabels([ ])
    ax['(b)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)

    # Septal inputs
    if net.two_septal_inputs:
        ax['(c)'].plot(net.tvec_ms_180,net.idvec_ms_180,'o',color='tab:orange',label='to basket',markersize=1.5)
        ax['(c)'].plot(net.tvec_ms_360,net.idvec_ms_360,'o',color='tab:green', label='to olm',markersize=1.5)
    else:
        if np.logical_and(net.MSGain_basket != 0, net.MSGain_olm !=0):
            label='to basket and olm'
        elif np.logical_and(net.MSGain_basket !=0, net.MSGain_olm == 0):
            label='to basket'
        elif np.logical_and(net.MSGain_basket == 0 , net.MSGain_olm != 0):
            label='to olm'
        ax['(c)'].plot(net.tvec_ms,np.ones(len(net.tvec_ms)),'o',color='black',label=label)
    ax['(c)'].set_ylabel('MS inputs')
    ax['(c)'].set_xticklabels([ ])
    ax['(c)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)

    # DG inpunts
    ax['(d)'].plot(net.tvec_dg_burst,  net.idvec_dg_burst,'o',color='tab:red',label='>40 Hz',markersize=1.5)
    ax['(d)'].plot(net.tvec_dg_regular,net.idvec_dg_regular,'o',color='tab:purple',label='4 Hz',markersize=1.5)
    ax['(d)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)
    ax['(d)'].set_ylabel('DG inputs')
    ax['(d)'].set_xticklabels([ ])

    # EC2 inputs
    ax['(e)'].plot(net.tvec_ec2_180,net.idvec_ec2_180,'o',color='tab:red',label='180',markersize=1.5)
    ax['(e)'].plot(net.tvec_ec2_360,net.idvec_ec2_360,'o',color='tab:purple',label='360',markersize=1.5)
    ax['(e)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)
    ax['(e)'].set_ylabel('EC2 inputs')
    ax['(e)'].set_xticklabels([ ])

    # EC3 inputs
    if net.create_ca1_network:
        ax['(f)'].plot(net.tvec_ec3_180,net.idvec_ec3_180,'o',color='tab:red',label='180',markersize=1.5)
        ax['(f)'].plot(net.tvec_ec3_360,net.idvec_ec3_360,'o',color='tab:purple',label='360',markersize=1.5)
        ax['(f)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)
        ax['(f)'].set_ylabel('EC3 inputs')
        #ax['(f)'].set_xticklabels([ ])
        ax['(f)'].set_xlabel('Time [ms]')

    for label in ax.keys():
        ax[label].set_xlim([tmin,tmax])
        ax[label].grid(True)
    fig.align_labels()
    plt.savefig(figtitle+".png", bbox_inches='tight')
    plt.close()

def plot_lfp(net,x,y,xmin,xmax,figtitle):
    fig, ax = plt.subplot_mosaic([['(a)','(a)','(b)','(b)'], ['(c)','(c)','(d)','(d)']],
                                    constrained_layout=True, figsize=(16,5))

    window_plot = np.logical_and(x>=xmin, x<=xmax)
    if net.create_ca1_network:
        ax['(a)'].plot(x[window_plot],y["ca1"][window_plot],color='black')
        ax['(a)'].set_ylabel('LFP [mV]')
        ax['(a)'].set_xlabel('Time [ms]')
        ax['(a)'].grid(True)
        ax['(a)'].set_xlim([xmin,xmax])

        xf, yf = tappers(x[window_plot],y["ca1"][window_plot])
        ax['(c)'].plot(xf,yf, color='black')
        ax['(c)'].set_xlabel("Frequency [Hz]")
        ax['(c)'].set_ylabel("Power [a.u.]")
        ax['(c)'].set_xlim([0,80])
        #ax['(c)'].set_yscale('log')
        ax['(c)'].grid(True)

    ax['(b)'].plot(x[window_plot],y["ca3"][window_plot],color='black')
    ax['(b)'].set_xlabel('Time [ms]')
    ax['(b)'].grid(True)
    ax['(b)'].set_xlim([xmin,xmax])

    xf,yf = tappers(x[window_plot],y["ca3"][window_plot])
    ax['(d)'].plot(xf,yf, color='black')
    ax['(d)'].set_xlabel("Frequency [Hz]")
    ax['(d)'].set_xlim([0,80])
    #ax['(d)'].set_yscale('log')
    ax['(d)'].grid(True)
    #ax[5].set_ylim([1e-4, np.max(y)*1.05])
    fig.align_labels()
    plt.savefig(figtitle+".png",bbox_inches='tight')
    plt.close()

def incoming_currents_plots(net,x,xmin,xmax,figtitle):
    window_plot = np.logical_and(x>=xmin, x<=xmax)
    loc = "center left"
    ax_list = [['(a)'], ['(b)'], ['(c)'], ['(d)'], ['(e)'], ['(f)'], ['(g)']]

    regions, region_labels = [net.pop_ca3],["ca3"]
    if net.create_ca1_network:
        regions, region_labels = [net.pop_ca3, net.pop_ca1], ["ca3, ca1"]
    for region, region_label in zip(regions, region_labels):
        fig, ax = plt.subplot_mosaic(ax_list, constrained_layout=True, figsize=(16,8))
        #labels: Adend3,Adend2,Adend1,soma,Bdend, Basket, Olm (#7)
        #labels: Adend3,Adend2,Adend1,soma,Bdend, Basket, Olm (#7, and eventually 8 because of CCK-inhibitory neurons)
        color = {"AMPA": "tab:blue", "NMDA": "tab:green", "GABA": "red", "all": "black"}
        po = region[0]
        index = 0
        seclist = ["Adend3", "Adend2", "Adend1", "soma", "Bdend"]
        labels = seclist + ["Basket", "OLM"]
        for sec,label in zip(seclist,ax_list[:5]):
            for key, isyn in po.__dict__[sec+"_isyn"].items():
                for keyc in color.keys():
                    if keyc in key:
                        color_=color[keyc]
                ax[label[0]].plot(x[window_plot], isyn[window_plot], label=key, color=color_)
            ax[label[0]].set_ylabel(labels[index])
            ax[label[0]].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1), fontsize=8)
            index+=1

        for po,label in zip(region[1:], ax_list[5:]):
            seclist = po.cell[0].sec_list
            for sec in seclist:
                for key, isyn  in po.__dict__[sec+"_isyn"].items():
                    for keyc in color.keys():
                        if keyc in key:
                            color_=color[keyc]
                            if key == "GABAss":
                                color_="tab:orange"
                    ax[label[0]].plot(x[window_plot], isyn[window_plot], label=key, color=color_)
                ax[label[0]].set_ylabel(labels[index])
                ax[label[0]].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1),fontsize=8)
                index+=1

        for label in ax_list[:6]:
            ax[label[0]].set_xticklabels([ ])
        for label in ax_list:
            ax[label[0]].set_xlim([xmin,xmax])
            ax[label[0]].grid(True)

        ax[label[0]].set_xlabel('Time [ms]')
        fig.align_labels()
        plt.savefig(figtitle+'_'+region_label+".png", bbox_inches='tight')
        plt.close()

def spiking_distributions(net, spikes, tmin, tf, figtitle="spiking_distribution.png"):
    data_spikes = dict.fromkeys(spikes.keys())
    for region in data_spikes.keys():
        data_spikes[region] = {}
        for process_data in spikes[region]:
            data_spikes[region].update(process_data)

    t0=tmin
    binsize=50.0#18 #18 da dt = 6.25 ms sacado del paper

    pd = dict.fromkeys(net.cell_labels)
    pd_inputs = dict.fromkeys(net.input_labels)
    pd_sum = dict.fromkeys(net.cell_labels)

    t_ref = net.input_time[-3] #ec3_360
    dt = binsize/360.*net.theta_rythm
    tf = tf-np.mod(tf-t0,dt)
    ntime = int((tf-t0)/dt)+1
    time_phase = np.linspace(t0,tf,ntime)

    spiking_array = np.zeros((net.n, len(time_phase)))
    regions = ["ca3"]
    if net.create_ca1_network:
        regions = ["ca3, ca1"]
    for ii, region in enumerate(regions):
        pyr, bas, olm = [],[],[]
        for i, spike_times in data_spikes[region].items():
            t = np.array(spike_times)
            for tt in t[t>t0]:
                if np.logical_and(i>=1200*ii, i<1200*ii+800):
                    pyr.append(tt)
                if np.logical_and(i>=1200*ii+800, i<1200*ii+1000):
                    bas.append(tt)
                if np.logical_and(i>=1200*ii+1000,i<1200*ii+1200):
                    olm.append(tt)

        bins = time_phase
        pd["pyr_"+region],x  = np.histogram(pyr, bins=bins)
        pd["bas_"+region],x  = np.histogram(bas, bins=bins)
        pd["olm_"+region],x  = np.histogram(olm, bins=bins)

    pd_inputs["sep_180"],x    = np.histogram(net.tvec_ms_180,     bins=bins)#spikes_hist(net.tvec_ms_180, net.idvec_ms_180, 10, ntime)
    pd_inputs["sep_360"],x    = np.histogram(net.tvec_ms_360,     bins=bins)#spikes_hist(net.tvec_ms_360, net.idvec_ms_360, 10, ntime)
    pd_inputs["ec2_180"],x    = np.histogram(net.tvec_ec2_180,    bins=bins)#spikes_hist(net.tvec_ec2_180, net.idvec_ms_180, 100, ntime)
    pd_inputs["ec2_360"],x    = np.histogram(net.tvec_ec2_360,    bins=bins)#spikes_hist(net.tvec_ec2_360, net.idvec_ms_360, 100, ntime)
    pd_inputs["ec3_180"],x    = np.histogram(net.tvec_ec3_180,    bins=bins)#spikes_hist(net.tvec_ec2_180, net.idvec_ms_180, 100, ntime)
    pd_inputs["ec3_360"],x    = np.histogram(net.tvec_ec3_360,    bins=bins)#spikes_hist(net.tvec_ec2_360, net.idvec_ms_360, 100, ntime)
    pd_inputs["dg_regular"],x = np.histogram(net.tvec_dg_regular, bins=bins)#spikes_hist(net.tvec_dg_regular, net.idvec_dg_regular, 100, ntime)
    pd_inputs["dg_burst"],x   = np.histogram(net.tvec_dg_burst,   bins=bins)#spikes_hist(net.tvec_dg_burst, net.idvec_dg_burst, 100, ntime)
    x = x[:-1]
    time_cycle = []
    ncycles = []
    cycles = 1
    tx = t_ref
    while tx<=tf:
        if tx>=t0:
            time_cycle.append(tx-t0)
            ncycles.append(cycles)
        tx+=net.theta_rythm
        cycles+=1

    nspikes_cycle_cells  = dict.fromkeys(net.cell_labels)
    nspikes_cycle_inputs = dict.fromkeys(net.input_labels)
    for key in pd.keys():
        nspikes_cycle_cells[key] = []
    for key in pd_inputs.keys():
        nspikes_cycle_inputs[key] = []

    for t1,t2 in zip(time_cycle[:-1],time_cycle[1:]):
        window=np.logical_and(x>=t1,x<t2)
        for key in pd.keys():
            nspikes_cycle_cells[key].append(np.sum(pd[key][window]))
        for key in pd_inputs.keys():
            try:
                nspikes_cycle_inputs[key].append(np.sum(pd_inputs[key][window]))
            except:
                nspikes_cycle_inputs[key].append(0)

    loc = "center left"
    ax_list = [['(a)'], ['(b)'], ['(c)'], ['(d)'], ['(e)']]
    fig, ax = plt.subplot_mosaic(ax_list, constrained_layout=True, figsize=(16,5))
    if net.create_ca1_network:
        ax['(e)'].step(x,pd["pyr_ca1"],           where='mid', color='tab:blue',   label='pyr')
        ax['(b)'].step(x,pd_inputs["ec3_180"],    where='mid', color='tab:red',    label='180')
        ax['(b)'].step(x,pd_inputs["ec3_360"],    where='mid', color='tab:purple', label='360')
        ax['(e)'].fill_between(x,pd["pyr_ca1"],           step='mid', color='tab:blue',   alpha=0.4)
        ax['(b)'].fill_between(x,pd_inputs["ec3_180"],    step='mid', color='tab:red',    alpha=0.4)
        ax['(b)'].fill_between(x,pd_inputs["ec3_360"],    step='mid', color='tab:purple', alpha=0.4)
        ax['(b)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1))

    ax['(d)'].step(x,pd["pyr_ca3"],           where='mid', color='tab:blue',   label='pyr')
    ax['(c)'].step(x,pd_inputs["dg_regular"], where='mid', color='tab:red',    label='4 Hz')
    ax['(c)'].step(x,pd_inputs["dg_burst"],   where='mid', color='tab:purple', label='> 40 Hz')
    ax['(a)'].step(x,pd_inputs["ec2_180"],    where='mid', color='tab:red',    label='180')
    ax['(a)'].step(x,pd_inputs["ec2_360"],    where='mid', color='tab:purple', label='360')
    ax['(d)'].fill_between(x,pd["pyr_ca3"],           step='mid', color='tab:blue',   alpha=0.4)
    ax['(c)'].fill_between(x,pd_inputs["dg_regular"], step='mid', color='tab:red',    alpha=0.4)
    ax['(c)'].fill_between(x,pd_inputs["dg_burst"],   step='mid', color='tab:purple', alpha=0.4)
    ax['(a)'].fill_between(x,pd_inputs["ec2_180"],    step='mid', color='tab:red',    alpha=0.4)
    ax['(a)'].fill_between(x,pd_inputs["ec2_360"],    step='mid', color='tab:purple', alpha=0.4)

    ax['(a)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1))
    ax['(c)'].legend(loc=loc, bbox_to_anchor=(1.02, 0.5,0.7,0.1))

    labels = ['EC2', 'EC3', 'DG', 'CA3', 'CA1']
    for i,key in enumerate(ax_list):
        ax[key[0]].set_ylabel(labels[i])
        ax[key[0]].grid(True)
    for key in ax_list[:4]:
        ax[key[0]].set_xticklabels([ ])
    ax['(e)'].set_xlabel('Time [ms]')
    fig.align_labels()
    plt.savefig(figtitle+'.png',bbox_inches='tight')
    plt.close()

################################################################################
# Implemented 12/06/2022 (updated 18/08)
################################################################################
def process_spike_data( data_unified, net, only_ca1=False): # Version actualizada 13/06
    npyr_ca3 = net.n_pyr_ca3 
    nbas_ca3 = net.n_bas_ca3
    nolm_ca3 = net.n_olm_ca3
    npyr_ca1 = net.n_pyr_ca1
    nbas_ca1 = net.n_bas_ca1
    nolm_ca1 = net.n_olm_ca1
    ncck_ca1 = net.n_cck_ca1
    nca3 = npyr_ca3+nbas_ca3+nolm_ca3
    nca1 = npyr_ca1+nbas_ca1+nolm_ca1+ncck_ca1

    if only_ca1:
        keys = ["pyr","bas","olm","cck"]
        neurons = dict.fromkeys(keys)
        neurons["pyr"] = np.arange(npyr_ca1)
        neurons["bas"] = np.arange(npyr_ca1,npyr_ca1+nbas_ca1)
        neurons["olm"] = np.arange(npyr_ca1+nbas_ca1,npyr_ca1+nbas_ca1+nolm_ca1)
        neurons["cck"] = np.arange(npyr_ca1+nbas_ca1+nolm_ca1,nca1)
        data={}
        for process_data in data_unified["ca1"]:
            data.update(process_data)
        spikes = {}
        spikes["ca1"] = dict.fromkeys(keys)
        for key in keys:
            spikes["ca1"][key] = dict.fromkeys(["tvec","idvec"])
            spikes["ca1"][key]["idvec"] = []
            spikes["ca1"][key]["tvec"]  = []
            for i in neurons[key]:
                if data[i]:
                    n = len(data[i])
                    spikes["ca1"][key]["idvec"].append( i*np.ones(n) )
                    spikes["ca1"][key]["tvec"].append( np.sort(data[i]) )

            if spikes["ca1"][key]["idvec"]:
                spikes["ca1"][key]["idvec"] = np.concatenate(spikes["ca1"][key]["idvec"]).astype(int)
                spikes["ca1"][key]["tvec"]  = np.round(np.concatenate(spikes["ca1"][key]["tvec"]),1)
                #spikes["ca1"][key]    = np.column_stack(spikes["ca1"][key])
    else:
        keys = ["pyr","bas","olm"]
        neurons = dict.fromkeys(keys)
        neurons["pyr"] = np.arange(npyr_ca3)
        neurons["bas"] = np.arange(npyr_ca3,npyr_ca3+nbas_ca3)
        neurons["olm"] = np.arange(npyr_ca3+nbas_ca3,nca3)

        data={}
        for process_data in data_unified["ca3"]:
            data.update(process_data)
        spikes = dict.fromkeys(["ca3","ca1"])
        spikes["ca3"] = dict.fromkeys(keys)
        for key in keys:
            spikes["ca3"][key] = dict.fromkeys(["tvec","idvec"])
            spikes["ca3"][key]["idvec"] = []
            spikes["ca3"][key]["tvec"]  = []
            for i in neurons[key]:
                if data[i]:
                    n = len(data[i])
                    spikes["ca3"][key]["idvec"].append(i*np.ones(n))
                    spikes["ca3"][key]["tvec"].append( np.sort(data[i]) )

            if spikes["ca3"][key]["tvec"]:
                spikes["ca3"][key]["idvec"] = np.concatenate(spikes["ca3"][key]["idvec"]).astype(int)
                spikes["ca3"][key]["tvec"]  = np.round(np.concatenate(spikes["ca3"][key]["tvec"]),1)
                #spikes["ca3"][key]    = np.column_stack(spikes["ca3"][key])

        if net.create_ca1_network:
            keys = ["pyr","bas","olm","cck"]
            neurons = dict.fromkeys(keys)
            neurons["pyr"] = np.arange(nca3,nca3+npyr_ca1)
            neurons["bas"] = np.arange(nca3+npyr_ca1,nca3+npyr_ca1+nbas_ca1)
            neurons["olm"] = np.arange(nca3+npyr_ca1+nbas_ca1,nca3+npyr_ca1+nbas_ca1+nolm_ca1)
            neurons["cck"] = np.arange(nca3+npyr_ca1+nbas_ca1+nolm_ca1,nca3+nca1)
            
            data={}
            for process_data in data_unified["ca1"]:
                data.update(process_data)
            spikes["ca1"] = dict.fromkeys(keys)
            for key in keys:
                spikes["ca1"][key] = dict.fromkeys(["tvec","idvec"])
                spikes["ca1"][key]["idvec"] = []
                spikes["ca1"][key]["tvec"]  = []
                for i in neurons[key]:
                    if data[i]:
                        n = len(data[i])
                        spikes["ca1"][key]["idvec"].append( i*np.ones(n) )
                        spikes["ca1"][key]["tvec"].append( np.sort(data[i]) )

                if spikes["ca1"][key]["idvec"]:
                    spikes["ca1"][key]["idvec"] = np.concatenate(spikes["ca1"][key]["idvec"]).astype(int)
                    spikes["ca1"][key]["tvec"]  = np.round(np.concatenate(spikes["ca1"][key]["tvec"]),1)
                    #spikes["ca1"][key]    = np.column_stack(spikes["ca1"][key])
    return spikes

def process_volt_data( data_unified, net, only_ca1 = False):
    npyr_ca3 = net.n_pyr_ca3 
    nbas_ca3 = net.n_bas_ca3
    nolm_ca3 = net.n_olm_ca3
    npyr_ca1 = net.n_pyr_ca1
    nbas_ca1 = net.n_bas_ca1
    nolm_ca1 = net.n_olm_ca1
    ncck_ca1 = net.n_cck_ca1
    nca3 = npyr_ca3+nbas_ca3+nolm_ca3
    nca1 = npyr_ca1+nbas_ca1+nolm_ca1+ncck_ca1

    if only_ca1:
        keys = ["pyr","bas","olm","cck"]
        neurons = dict.fromkeys(keys)
        neurons["pyr"] = np.arange(npyr_ca1)
        neurons["bas"] = np.arange(npyr_ca1,npyr_ca1+nbas_ca1)
        neurons["olm"] = np.arange(npyr_ca1+nbas_ca1,npyr_ca1+nbas_ca1+nolm_ca1)
        neurons["cck"] = np.arange(npyr_ca1+nbas_ca1+nolm_ca1,nca1)

        data={}
        for process_data in data_unified["ca1"]:
            data.update(process_data)
        volt = {}
        volt["ca1"] = dict.fromkeys(keys)
        for key in keys:
            volt["ca1"][key] = []
            volt_aux = []
            for i in neurons[key]:
                volt_aux.append(data[i])
            volt["ca1"][key].append( np.mean(volt_aux,axis=0) )
            volt["ca1"][key].append( np.std(volt_aux,axis=0) )

    else:

        keys = ["pyr","bas","olm"]
        neurons = dict.fromkeys(keys)
        ns = np.arange(1200)
        neurons["pyr"] = np.arange(npyr_ca3)
        neurons["bas"] = np.arange(npyr_ca3,npyr_ca3+nbas_ca3)
        neurons["olm"] = np.arange(npyr_ca3+nbas_ca3,nca3)

        data={}
        for process_data in data_unified["ca3"]:
            data.update(process_data)

        volt=dict.fromkeys(["ca3","ca1"])
        volt["ca3"] = dict.fromkeys(keys)
        for key in keys:
            volt["ca3"][key] = []
            volt_aux = []
            for i in neurons[key]:
                volt_aux.append(data[i])
            volt["ca3"][key].append( np.mean(volt_aux,axis=0) )
            volt["ca3"][key].append( np.std(volt_aux,axis=0) )
        print(volt.keys())

        if net.create_ca1_network:
            keys = ["pyr","bas","olm","cck"]
            neurons = dict.fromkeys(keys)
            neurons["pyr"] = np.arange(nca3,nca3+npyr_ca1)
            neurons["bas"] = np.arange(nca3+npyr_ca1,nca3+npyr_ca1+nbas_ca1)
            neurons["olm"] = np.arange(nca3+npyr_ca1+nbas_ca1,nca3+npyr_ca1+nbas_ca1+nolm_ca1)
            neurons["cck"] = np.arange(nca3+npyr_ca1+nbas_ca1+nolm_ca1,nca3+nca1)

            data={}
            for process_data in data_unified["ca1"]:
                data.update(process_data)

            volt["ca1"] = dict.fromkeys(keys)
            for key in keys:
                volt["ca1"][key] = []
                volt_aux = []
                for i in neurons[key]:
                    volt_aux.append(data[i])
                volt["ca1"][key].append( np.mean(volt_aux,axis=0) )
                volt["ca1"][key].append( np.std(volt_aux,axis=0) )
    return volt

#####################################################################################
## November 2022
#####################################################################################
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
    kurtosis = []
    for x, y in zip(xlist, ylist):
        t1, t2 = initial_time, initial_time+theta_cycle
        
        while t2<=final_time:
            w = np.logical_and(x>=t1, x<t2)
            nspikes.append( len(x[w]) ) 
            active_neurons.append( len(np.unique(y[w])) )
            
            kurtosis.append( scipy.stats.kurtosis( x[w] ))
            yaux = Counter(y[w]) 
            yaux = np.array(list(yaux.values()))
            active_same_neurons.append( len(yaux[np.where(yaux>1)]) )
            
            t1+=theta_cycle
            t2+=theta_cycle
    
    n = len(nspikes)
    nspikes_std = np.round( np.std(nspikes)/np.sqrt(n),2)
    nspikes_mean = np.round( np.mean(nspikes), 2)
    active_neurons_std = np.round( np.std(active_neurons)/np.sqrt(n), 2)
    active_neurons_mean = np.round( np.mean(active_neurons), 2)
    active_same_neurons_std = np.round( np.std(active_same_neurons)/np.sqrt(n), 2)
    active_same_neurons_mean = np.round( np.mean(active_same_neurons), 2)  
    kurtosis_std  = np.round( np.std(kurtosis)/np.sqrt(n), 2)
    kurtosis_mean = np.round( np.mean(kurtosis), 2)
    
    return nspikes_mean, nspikes_std, active_neurons_mean, active_neurons_std, active_same_neurons_mean, active_same_neurons_std, kurtosis_mean, kurtosis_std

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

    vdiff_mean = np.mean( vdiff )
    vdiff_sem  = np.std( vdiff )/np.sqrt(len(vdiff))
    vmin_mean  = np.mean(vmin)                                     
    vmin_sem   = np.std(vmin)/np.sqrt(len(vdiff))
    vmax_mean  = np.mean(vmax)                                     
    vmax_sem   = np.std(vmax)/np.sqrt(len(vdiff))
    return vmin_mean, vmin_sem, vmax_mean, vmax_sem, vdiff_mean, vdiff_sem