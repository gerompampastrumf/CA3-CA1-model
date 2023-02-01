import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy 
import copy 
'''
from neurodsp.filt import filter_signal
from neurodsp.plts import plot_time_series
from neurodsp.sim import sim_combined
from neurodsp.utils import set_random_seed, create_times

from bycycle.features import compute_features
from bycycle.cyclepoints import find_extrema, find_zerox
from bycycle.cyclepoints.zerox import find_flank_zerox
from bycycle.plts import plot_burst_detect_summary, plot_cyclepoints_array
from bycycle.utils.download import load_bycycle_data
'''
import os 
import sys
# import file_management
from spectrum import *
from scipy.signal import hilbert, chirp 

def hilbert_phase(x):
    sig=x-x.mean()
    std = sig.std()
    sig/=std
    analytic_sig = hilbert(sig)
    instantaneous_phase = np.angle(analytic_sig)
    amplitude_envelope  = np.abs(analytic_sig)
    return np.mod(instantaneous_phase,2*np.pi), (amplitude_envelope+x.mean())*std

def get_theta_phase_old(xs, fs, f0, df=1.0, n_seconds=1, plot_filter_properties=False): 
    # signal: original signal
    # fs: sampling frequency 
    # f0: theta frequency 
    # df: [f0-df,f0+df] filter bandwith
    
    f_range = (f0-df, f0+df)
    sig_band = filter_signal(xs, fs, "bandpass", f_range, n_seconds=n_seconds)
    no_nan_window = ~np.isnan(sig_band)
    phase = np.ones(len(xs))*np.nan
    phase_, envelope_ = hilbert_phase(sig_band[no_nan_window])
    phase[no_nan_window] = phase_
    
    #if plot_filter_properties:
    #    print("Theta filter:")
    #    sig_filt = filter_signal(sig_band, fs, "bandpass", f_range,n_seconds=n_seconds, plot_properties=plot_filter_properties)
    
    return sig_band+xs.mean(), phase
    
def get_gamma_amplitude_old(xs, fs, f0, df=10, n_seconds=1, plot_filter_properties=False, return_nans=False):
    f_range = (f0-df, f0+df)
    sig_band = filter_signal(xs, fs, "bandpass", f_range, n_seconds=n_seconds)
    no_nan_window = ~np.isnan(sig_band)
    envelope = np.ones(len(xs))*np.nan
    phase_, envelope_ = hilbert_phase(sig_band[no_nan_window])
    envelope[no_nan_window] = envelope_
    
    if plot_filter_properties:
        print("Gamma filter:")
        sig_filt = filter_signal(sig_band, fs, "bandpass", f_range,n_seconds=0.1, plot_properties=plot_filter_properties)
        
    if return_nans:
        return sig_band, envelope,  no_nan_window
    else:   
        return sig_band, envelope

def modulation_index_old(theta_phase, gamma_amplitude, fs=1e4, nbins=20, nsurrogates=100,surrogate_test=False): 

    if np.mod(nbins,4) == 0:
        nbins = nbins + (4-np.mod(nbins,4)) # Increase nbins to the next multiple of four
        
    n = len(theta_phase)
    time = np.linspace(0,(n-1)*fs*1e-3,n)
    
    starting_cycles = np.where( np.abs(np.diff( theta_phase ) )> 1.9*np.pi)[0]
    ncycles = len(starting_cycles)-1

    amplitude = np.zeros((ncycles, nbins))
    for i in range(ncycles): 
        n1,n2 = starting_cycles[i]+1, starting_cycles[i+1]
        t1,t2 = time[n1],time[n2]
                               
        phi1, phi2 = theta_phase[n1], theta_phase[n2]
                               
        time_aux = np.linspace(t1,t2,nbins+1)
        phi_aux  = np.linspace(phi1,phi2, nbins+1)
                               
        dphi = np.diff(phi_aux[:2])
        phi_seq = phi_aux[:-1]+dphi/2.0
                               
        for j in range(len(time_aux[:-1])): 
            tj1,tj2 = time_aux[j],time_aux[j+1]
                               
            w= np.logical_and(time>=tj1, time<tj2)
            amplitude[i,j] = np.mean(gamma_amplitude[w])
    
    amplitude = np.mean(amplitude,axis=0)
    pi = amplitude/np.sum(amplitude)
    entropy_max = np.log(nbins)
    entropy = -np.sum(pi*np.log(pi))
    modulation_index = (entropy_max-entropy)/entropy_max
    
    ### surrogate_text
    if surrogate_test:
        np.random.seed(1003)
        surrogates = np.full(( nsurrogates, len(theta_phase)), np.nan)
        no_nan_window = ~np.isnan(gamma_amplitude)
        gamma_amplitude_aux = gamma_amplitude[no_nan_window]
        ng = len(gamma_amplitude_aux)
        for s in range(nsurrogates):
            ic = np.random.randint(ng)
            surrogates[s][no_nan_window] = np.array( gamma_amplitude_aux[ic:].tolist()+ gamma_amplitude_aux[:ic].tolist())
        
        modulation_index_surrogates = np.zeros(nsurrogates)
        pi_surrogates = np.zeros((nsurrogates,nbins))
        for s in range(nsurrogates):
                            
            amplitude = np.zeros((ncycles, nbins))
            for i in range(ncycles): 
                n1,n2 = starting_cycles[i]+1, starting_cycles[i+1]
                t1,t2 = time[n1],time[n2]

                phi1, phi2 = theta_phase[n1], theta_phase[n2]

                time_aux = np.linspace(t1,t2,nbins+1)
                phi_aux  = np.linspace(phi1,phi2, nbins+1)

                dphi = np.diff(phi_aux[:2])
                phi_seq = phi_aux[:-1]+dphi/2.0

                for j in range(len(time_aux[:-1])): 
                    tj1,tj2 = time_aux[j],time_aux[j+1]

                    w= np.logical_and(time>=tj1, time<tj2)
                    amplitude[i,j] = np.mean(surrogates[s][w])
                            
            amplitude = np.mean(amplitude,axis=0)
            pi_surrogates[s] = amplitude/np.sum(amplitude)
            entropy = -np.sum(pi_surrogates[s]*np.log(pi_surrogates[s]))
            modulation_index_surrogates[s] = (entropy_max-entropy)/entropy_max
        
        pvalue = 1-len(modulation_index_surrogates[modulation_index_surrogates<= modulation_index])/nsurrogates
        plt.figure()
        plt.hist(modulation_index_surrogates)
        if pvalue <= 0.05: 
            modulation_index_sig = modulation_index
            plt.axvline(modulation_index,color='red')
            plt.title("pvlue less than 0.05")
        else: 
            modulation_index_sig = 0
            plt.title("statistically insignificant") 
            plt.axvline(modulation_index,color='red')

        return modulation_index , modulation_index_sig
    else: 
        return modulation_index

def compute_mi(timeseries, ftheta, fgamma, fs=1e3, nbins=20, surrogate_test=False, nsurrogates=100, surrogate_method = "block boostrapping", surrogate_nsplits=99): 
    
    ntheta,  ngamma = len(ftheta), len(fgamma)
    
    mi = np.zeros((ntheta,ngamma))
    mi_surrogates = np.zeros((nsurrogates, ntheta, ngamma))

    sig_band1 = [ [] for i in range(ntheta)]
    sig_band2 = [ [] for i in range(ngamma)]
    phase     = [ [] for i in range(ntheta)] 
    envelope  = [ [] for i in range(ngamma)]
    
    norder = 4
    for i,fth in enumerate(ftheta):
        sig_band1[i], _, phase[i] = bandpass_filter_and_hilbert_transform(timeseries, fs=fs, f0=fth, df=1.0, norder=norder)
    for i,fg in enumerate(fgamma):
        sig_band2[i], envelope[i], _ = bandpass_filter_and_hilbert_transform(timeseries, fs=fs, f0=fg, df=10.0, norder=norder)

    entropy_max = np.log(nbins)
    
    np.random.seed(19)
    for i in range(ntheta):
        theta_phase = phase[i]
        for j in range(ngamma):
            gamma_amplitude = envelope[j]
            
            pbins = np.linspace(0,2*np.pi,nbins+1)
            amplitude = np.zeros(len(pbins)-1)       
            for k in range(len(pbins)-1):         
                pl = pbins[k]
                pr = pbins[k+1]                   
                indices=(theta_phase>=pl) & (theta_phase<pr)  
                amplitude[k] = np.mean(gamma_amplitude[indices]) 

            pi = amplitude/np.sum(amplitude)
            entropy = -np.sum(pi*np.log(pi))
            modulation_index = (entropy_max-entropy)/entropy_max
            mi[i,j] = modulation_index
            
            if surrogate_test: 
                surrogates = surrogate_generation( gamma_amplitude, nsurrogates=nsurrogates, method = surrogate_method, nsplits=surrogate_nsplits)
                for l, surrogate in enumerate(surrogates): 
                    amplitude = np.zeros(len(pbins)-1)       
                    for k in range(len(pbins)-1):  
                        pl = pbins[k]
                        pr = pbins[k+1]                   
                        indices=(theta_phase>=pl) & (theta_phase<pr) 

                        amplitude[k] = np.mean(surrogate[indices]) 
                        
                    pi = amplitude/np.sum(amplitude)
                    entropy = -np.sum(pi*np.log(pi))
                    modulation_index = (entropy_max-entropy)/entropy_max
                    mi_surrogates[l,i,j] = modulation_index
                    
    return ftheta, fgamma, mi, mi_surrogates

# new 15/01/2022
def compute_mi_parallel(timeseries, ftheta, fgamma, fs=1e3, nbins=20, surrogate_test=False, nsurrogates=100, surrogate_method = "block boostrapping", surrogate_nsplits=99): 
    
    ngamma = len(fgamma)
    mi = np.zeros(ngamma)
    mi_surrogates = np.zeros((nsurrogates, ngamma))

    sig_band2 = [ [] for i in range(ngamma)]
    envelope  = [ [] for i in range(ngamma)]
    
    norder = 4
    sig_band1, _, phase = bandpass_filter_and_hilbert_transform(timeseries, fs=fs, f0=ftheta, df=1.0, norder=norder)
    for i,fg in enumerate(fgamma):
        sig_band2[i], envelope[i], _ = bandpass_filter_and_hilbert_transform(timeseries, fs=fs, f0=fg, df=10.0, norder=norder)
        
    entropy_max = np.log(nbins)
    np.random.seed(19)
    theta_phase = phase
    
    for j in range(ngamma):
        gamma_amplitude = envelope[j]
        pbins = np.linspace(0,2*np.pi,nbins+1)
        amplitude = np.zeros(len(pbins)-1)       
        for k in range(len(pbins)-1):         
            pl = pbins[k]
            pr = pbins[k+1]                   
            indices=(theta_phase>=pl) & (theta_phase<pr)  
            amplitude[k] = np.mean(gamma_amplitude[indices]) 

        pi = amplitude/np.sum(amplitude)
        entropy = -np.sum(pi*np.log(pi))
        modulation_index = (entropy_max-entropy)/entropy_max
        mi[j] = modulation_index

        if surrogate_test: 
            surrogates = surrogate_generation( gamma_amplitude, nsurrogates=nsurrogates, method = surrogate_method, nsplits = surrogate_nsplits)
            for l, surrogate in enumerate(surrogates): 
                amplitude = np.zeros(len(pbins)-1)       
                for k in range(len(pbins)-1):  
                    pl = pbins[k]
                    pr = pbins[k+1]                   
                    indices=(theta_phase>=pl) & (theta_phase<pr) 

                    amplitude[k] = np.mean(surrogate[indices]) 

                pi = amplitude/np.sum(amplitude)
                entropy = -np.sum(pi*np.log(pi))
                modulation_index = (entropy_max-entropy)/entropy_max
                mi_surrogates[l,j] = modulation_index
                
    return mi, mi_surrogates 

def compute_mi_lists(timeseries, ftheta, fgamma, fs=1e3, nbins=20, surrogate_test=False, nsurrogates=100, surrogate_method = "block boostrapping", surrogate_nsplits=99): 
    '''
    timeseries is now a list
    '''
    ntrials = len(timeseries)
    ntheta,  ngamma = len(ftheta), len(fgamma)
    mi = np.zeros((ntheta,ngamma))
    mi_surrogates = np.zeros((nsurrogates, ntheta, ngamma))

    sig_band1 = [ [ [] for j in range(ntrials) ] for i in range(ntheta)]
    sig_band2 = [ [ [] for j in range(ntrials) ] for i in range(ngamma)]
    phase     = [ [ [] for j in range(ntrials) ] for i in range(ntheta)] 
    envelope  = [ [ [] for j in range(ntrials) ] for i in range(ngamma)]
    
    norder = 4
    for i,fth in enumerate(ftheta):
        for j in range(ntrials):
            sig_band1[i][j], _, phase[i][j] = bandpass_filter_and_hilbert_transform(timeseries[j], fs=fs, f0=fth, df=1.0, norder=norder)
    for i,fg in enumerate(fgamma):
        for j in range(ntrials):
            sig_band2[i][j], envelope[i][j], _ = bandpass_filter_and_hilbert_transform(timeseries[j], fs=fs, f0=fg, df=10.0, norder=norder)

    entropy_max = np.log(nbins)
    
    np.random.seed(19)
    pbins = np.linspace(0,2*np.pi,nbins+1)
    for i in range(ntheta):
        #theta_phase = phase[i]
        for j in range(ngamma):
            #gamma_amplitude = envelope[j]
            amplitude = np.zeros((ntrials,nbins))
            for k in range(ntrials):
                theta_phase = phase[i][k]
                gamma_amplitude = envelope[j][k]
                for l in range(len(pbins)-1):         
                    pl = pbins[l]
                    pr = pbins[l+1]                   
                    indices=(theta_phase>=pl) & (theta_phase<pr)  
                    amplitude[k,l] = np.mean(gamma_amplitude[indices]) 
            amplitude = np.concatenate(amplitude, axis=0)
            pi = amplitude/np.sum(amplitude)
            entropy = -np.sum(pi*np.log(pi))
            modulation_index = (entropy_max-entropy)/entropy_max
            mi[i,j] = modulation_index
            
    if surrogate_test: 
        for i in range(ntheta):
            for j in range(ngamma):
                surrogates = surrogate_generation_list( envelope[j], nsurrogates=nsurrogates, method = surrogate_method, nsplits=surrogate_nsplits)
                for s, surrogate in enumerate(surrogates):
                    amplitude = np.zeros((ntrials, nbins))
                    for k in range(ntrials):
                        theta_phase = phase[i][k]
                        for l in range(len(pbins)-1):
                            pl = pbins[l]
                            pr = pbins[l+1]                   
                            indices=(theta_phase>=pl) & (theta_phase<pr)  
                            amplitude[k,l] = np.mean(surrogate[k][indices]) 

                    amplitude = np.concatenate(amplitude, axis=0)
                    pi = amplitude/np.sum(amplitude)
                    entropy = -np.sum(pi*np.log(pi))
                    modulation_index = (entropy_max-entropy)/entropy_max
                    mi[s,i,j] = modulation_index

    return ftheta, fgamma, mi, mi_surrogates

def compute_mi_lists_parallel(timeseries, ftheta, fgamma, fs=1e3, nbins=20, surrogate_test=False, nsurrogates=100, surrogate_method = "block boostrapping", surrogate_nsplits=99): 
    '''
    timeseries is now a list
    ftheta is scalar not a list/array
    '''
    ntrials = len(timeseries)
    ntheta,  ngamma = len(ftheta), len(fgamma)
    mi = np.zeros((ntheta,ngamma))
    mi_surrogates = np.zeros((nsurrogates, ntheta, ngamma))

    sig_band1 = [ [] for j in range(ntrials) ] 
    sig_band2 = [ [ [] for j in range(ntrials) ] for i in range(ngamma)]
    phase     = [ [] for j in range(ntrials) ] 
    envelope  = [ [ [] for j in range(ntrials) ] for i in range(ngamma)]
    
    norder = 4
    for j in range(ntrials):
        sig_band1[j], _, phase[j] = bandpass_filter_and_hilbert_transform(timeseries[j], fs=fs, f0=ftheta, df=1.0, norder=norder)
    for i,fg in enumerate(fgamma):
        for j in range(ntrials):
            sig_band2[i][j], envelope[i][j], _ = bandpass_filter_and_hilbert_transform(timeseries[j], fs=fs, f0=fg, df=10.0, norder=norder)

    entropy_max = np.log(nbins)
    np.random.seed(19)
    pbins = np.linspace(0,2*np.pi,nbins+1)
    # theta_phase = phase[i]
    for j in range(ngamma):
        # gamma_amplitude = envelope[j]
        amplitude = np.zeros((ntrials,nbins))
        for k in range(ntrials):
            theta_phase = phase[k]
            gamma_amplitude = envelope[j][k]
            for l in range(len(pbins)-1):         
                pl = pbins[l]
                pr = pbins[l+1]                   
                indices=(theta_phase>=pl) & (theta_phase<pr)  
                amplitude[k,l] = np.mean(gamma_amplitude[indices]) 
        amplitude = np.concatenate(amplitude, axis=0)
        pi = amplitude/np.sum(amplitude)
        entropy = -np.sum(pi*np.log(pi))
        modulation_index = (entropy_max-entropy)/entropy_max
        mi[j] = modulation_index
        
    if surrogate_test: 
        for j in range(ngamma):
            surrogates = surrogate_generation_list( envelope[j], nsurrogates=nsurrogates, method = surrogate_method, nsplits=surrogate_nsplits)
            for s, surrogate in enumerate(surrogates):
                amplitude = np.zeros((ntrials, nbins))
                for k in range(ntrials):
                    theta_phase = phase[k]
                    for l in range(len(pbins)-1):
                        pl = pbins[l]
                        pr = pbins[l+1]                   
                        indices=(theta_phase>=pl) & (theta_phase<pr)  
                        amplitude[k,l] = np.mean(surrogate[k][indices]) 

                amplitude = np.concatenate(amplitude, axis=0)
                pi = amplitude/np.sum(amplitude)
                entropy = -np.sum(pi*np.log(pi))
                modulation_index = (entropy_max-entropy)/entropy_max
                mi[s,j] = modulation_index

    return ftheta, fgamma, mi, mi_surrogates

#####################################################################################################
def compute_c_function(xseq, yseq): 
    a = np.sum( xseq*np.conj(yseq), axis=0 )
    b = np.sum( np.abs(xseq)**2, axis=0)
    c = np.sum( np.abs(yseq)**2, axis=0)
    return a/np.sqrt(b*c)

def compute_nu_function(xseq, yseq): 
    a = np.abs( np.sum( xseq*np.conj(yseq), axis=0 ) )
    b = np.sum( np.abs(xseq)**2, axis=0 )
    c = np.sum( np.abs(yseq)**2, axis=0 )
    return a/np.sqrt(b*c)

def compute_psi_old(xseq,yseq,iftheta, ibeta):
    psi = np.zeros(len(iftheta))
    c_function = compute_c_function(xseq,yseq)
    for ii,i in enumerate(iftheta):
        c1 = np.conj( c_function[i-ibeta//2:i+1+ibeta//2] )
        c2 = c_function[i+1-ibeta//2:i+2+ibeta//2]
        psi[ii] = np.imag( np.sum(c1*c2) )
    return psi

def compute_psi_function(xseq,yseq,n,fr):
    ibeta = int(2/np.diff(fr)[0])
    psi, fr_new = [],[]
    c_function = compute_c_function(xseq,yseq)
    # esto igual se puede mejorar haciendo algun producto matricial 
    for i in range(ibeta//2,n-ibeta//2-1): 
        c1 = np.conj( c_function[i-ibeta//2:i+1+ibeta//2] )
        c2 = c_function[i+1-ibeta//2:i+2+ibeta//2]
        psi.append(  np.imag( np.sum(c1*c2) ) )
        fr_new.append(fr[i])

    return np.array(fr_new), np.array(psi)

def plot_spectrum(s):
    f = np.fft.rfftfreq(len(s))
    plt.loglog(f, np.abs(np.fft.rfft(s)))

def noise_psd(N, psd = lambda f: 1):
        X_white = np.fft.rfft(np.random.randn(N))
        S = psd(np.fft.rfftfreq(N))
        # Normalize S
        S = S / np.sqrt(np.mean(S**2))
        X_shaped = X_white * S
        return np.fft.irfft(X_shaped)

def PSDGenerator(f):
    return lambda N: noise_psd(N, f)

@PSDGenerator
def white_noise(f):
    return 1

@PSDGenerator
def blue_noise(f):
    return np.sqrt(f)

@PSDGenerator
def violet_noise(f):
    return f

@PSDGenerator
def brownian_noise(f):
    return 1/np.where(f == 0, float('inf'), f)

@PSDGenerator
def pink_noise(f):
    return 1/np.where(f == 0, float('inf'), np.sqrt(f))

################################################################
## Con nuevos filtros ##
################################################################

def bandpass_filter(xs, norder, f_range, fs=1e4): 
    sos = scipy.signal.butter(N=norder, Wn=f_range, btype="bandpass", fs=fs, output="sos")
    return scipy.signal.sosfiltfilt(sos, xs)

def lowpass_filter(xs, norder, f_range, fs=1e4):
    sos = scipy.signal.butter(N=norder, Wn=f_range, btype="lowpass", fs=fs, output="sos")
    return scipy.signal.sosfiltfilt(sos, xs)

def get_theta_phase(sig, fs, f0, norder=4, df=1.0, n_seconds=1, plot_filter_properties=False): 
    # signal: original signal
    # fs: sampling frequency 
    # f0: theta frequency 
    # df: [f0-df,f0+df] filter bandwith
    
    f_range = (f0-df, f0+df)
    sig_band = bandpass_filter(sig,norder,f_range,fs)
    no_nan_window = ~np.isnan(sig_band)
    phase = np.ones(len(sig))*np.nan
    phase_, envelope_ = hilbert_phase(sig_band[no_nan_window])
    phase[no_nan_window] = phase_
    
    if plot_filter_properties:
        print("Theta filter:")
        sig_filt = filter_signal(sig_band, fs, "bandpass", f_range,n_seconds=n_seconds, plot_properties=plot_filter_properties)
    
    return sig_band+sig.mean(), phase
    
def get_gamma_amplitude(sig, fs, f0, norder=4,df=10, n_seconds=1, plot_filter_properties=False, return_nans=False):
    f_range = (f0-df, f0+df)
    sig_band = bandpass_filter(sig, norder, f_range,fs)
    no_nan_window = ~np.isnan(sig_band)
    envelope = np.ones(len(sig))*np.nan
    phase_, envelope_ = hilbert_phase(sig_band[no_nan_window])
    envelope[no_nan_window] = envelope_
    
    if plot_filter_properties:
        print("Gamma filter:")
        sig_filt = filter_signal(sig_band, fs, "bandpass", f_range,n_seconds=0.1, plot_properties=plot_filter_properties)
        
    if return_nans:
        return sig_band+sig.mean(), envelope,  no_nan_window
    else:   
        return sig_band+sig.mean(), envelope  
    
def bandpass_filter_and_hilbert_transform(sig, fs, f0, df, norder):
    scale = sig.mean() 
    f_range = (f0-df, f0+df)
    sig_band = bandpass_filter(sig, norder, f_range, fs)
    phase, envelope= hilbert_phase(sig_band)
    return sig_band+scale, envelope+scale,  phase

def get_nfft(m):
    nfft = 2
    k = 0 
    while m > nfft:
        k+=1 
        nfft = 2**k
    return nfft

def surrogate_generation(timeseries, nsurrogates, method ="block boostrapping" , replace = False, nsplits = 99): 
    surrogate_list = [] 
    n = len(timeseries)

    if method == "2block boostraping":
        for _ in range(nsurrogates):
            ir = np.random.randint(low=0,high=len(timeseries))
            surrogate = np.array( list(timeseries[ir:])+list(timeseries[:ir]) )
            surrogate_list.append( surrogate )
            
    if method == "sampling":
        for _ in range(nsurrogates):
            surrogate = np.random.choice(timeseries, size=n, replace=replace)
            surrogate_list.append( surrogate )
    
    if method == "fourier":
        # for _ in range(nsurrogates):
        #     fourier = np.fft.fft(timeseries)
        #     random_phases = 1j*np.random.choice(fourier.imag[:n//2], size=n//2, replace=replace)
        #     random_phases = np.concatenate(( random_phases, random_phases[::-1] ))
        #     fourier_new = fourier*random_phases
        #     surrogate = np.real( np.fft.ifft(fourier_new))
        #     surrogate_list.append( surrogate )
        
        for _ in range(nsurrogates):
            fourier = np.fft.rfft(timeseries)
            amplitude, phases= np.abs(fourier), np.angle(fourier)
            new_phases = copy.copy(phases)
            np.random.shuffle( new_phases )
            fourier_new = amplitude*new_phases
            surrogate = np.real( np.fft.ifft(fourier_new))
            surrogate_list.append( surrogate )

            
    if method == "block boostrapping":
        timeseries_splitted = np.array( np.array_split(timeseries, nsplits) ) 
        ns = len(timeseries_splitted)
        for _ in range(nsurrogates):
            indexes_surrogate = np.random.choice( np.arange(ns), size = ns, replace=replace)
            surrogate = np.concatenate( timeseries_splitted[indexes_surrogate] )
            surrogate_list.append( surrogate )
            
    return surrogate_list 

def surrogate_generation_list(timeseries, nsurrogates, method ="block boostrapping" , replace = False, nsplits = 99): 
    '''
    now timeseries is a list
    '''
    ntrials, n = np.shape( timeseries )
    surrogate_list = [ [] for k in range(ntrials) ] 

    if method == "2block boostraping":
        for k in range(ntrials):
            for _ in range(nsurrogates):
                ir = np.random.randint(low=0,high=len(timeseries[k]))
                surrogate = np.array( list(timeseries[k][ir:])+list(timeseries[k][:ir]) )
                surrogate_list[k].append( surrogate )
            
    if method == "sampling":
        for k in range(ntrials):
            for _ in range(nsurrogates):
                surrogate = np.random.choice(timeseries[k], size=n, replace=replace)
                surrogate_list[k].append( surrogate )
        
    if method == "fourier":
        for k in range(ntrials):
            for _ in range(nsurrogates):
                fourier = np.fft.rfft(timeseries[k])
                amplitude, phases= np.abs(fourier), np.angle(fourier)
                new_phases = copy.copy(phases)
                np.random.shuffle( new_phases )
                fourier_new = amplitude*new_phases
                surrogate = np.real( np.fft.ifft(fourier_new))
                surrogate_list[k].append( surrogate )

            
    if method == "block boostrapping":
        for k in range(ntrials):
            timeseries_splitted = np.array( np.array_split(timeseries[k], nsplits) ) 
            ns = len(timeseries_splitted)
            for _ in range(nsurrogates):
                indexes_surrogate = np.random.choice( np.arange(ns), size = ns, replace=replace)
                surrogate = np.concatenate( timeseries_splitted[indexes_surrogate] )
                surrogate_list[k].append( surrogate )
                
    return np.array(surrogate_list).T # (nsurrogates, ntrials)

def get_pvalue( samples,  popmean, method = "two-sided"):
    z0 = (popmean-samples.mean())/samples.std()
    if method == "two-sided":
        pvalue = scipy.stats.norm.sf(abs(z0))*2
    if method == "one-sided":
        pvalue = scipy.stats.norm.sf(abs(z0))
    return pvalue

def get_statistical_significance( samples, popmean, method = "two-sided"):
    
    if method == "two-sided":
        q1 = np.quantile( samples, 0.975, axis=0)
        q2 = np.quantile( samples, 0.025, axis=0)
        condition = np.logical_or( popmean >= q1, popmean <= q2)
        
    if method == "one-sided":
        q = np.quantile( samples , 0.95, axis = 0) # right
        condition = (popmean >= q)
    return condition

#################################################################################################
def compute_coherence_measurements(x, fs, nperseg, noverlap, ftheta_min=4, ftheta_max=20, fgamma_min=25, fgamma_max=100, dfgamma=2,
                                   norder=4, nfft=2048, surrogate_test=True, nsurrogates=100, surrogate_method = "block boostrapping", surrogate_nsplits=99):
    
    psi, nu = [],[] 
    psi_surrogates = [ [] for i in range(nsurrogates) ]
    nu_surrogates  = [ [] for i in range(nsurrogates) ]
    
    fgamma = np.arange(fgamma_min,fgamma_max,dfgamma)
    fr, pxx  = scipy.signal.csd( x, x, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
    
    np.random.seed(12)
    for fg in fgamma:
        yf, y, _ = bandpass_filter_and_hilbert_transform(x, norder=norder, fs=fs, f0=fg, df=10)
        fr, pyy  = scipy.signal.csd( y, y, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
        fr, pxy  = scipy.signal.csd( y, x, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap) 
        complex_coherence = pxy / np.sqrt( np.abs(pxx)*np.abs(pyy) ) 
        
        # compute of psi #################################################
        ibeta = 4
        imin = int( (ftheta_min-fr[0])/np.diff(fr)[0] )
        imax = int( (ftheta_max-fr[0])/np.diff(fr)[0] )
        cfd = []
        # module = np.abs(complex_coherence)
        # phase  = np.angle(complex_coherence)
        for i in range(imin, imax):  
            # factor1 = module[i-ibeta//2:i+ibeta//2]
            # factor2 = module[i+1-ibeta//2:i+1+ibeta//2]
            # factor3 = np.sin( phase[i+1-ibeta//2:i+1+ibeta//2]-phase[i-ibeta//2:i+ibeta//2])
            # cfd.append( np.sum(factor1*factor2*factor3) ) 
            factor1 = np.conj( complex_coherence[ i-ibeta//2: i+1+ibeta//2 ] )
            factor2 = complex_coherence[ i+1-ibeta//2 : i+2+ibeta//2 ]
            cfd.append(  np.imag( np.sum(factor1*factor2) ) )
            
        ##################################################################
        psi.append( np.array( cfd ) )
        nu.append( np.abs( complex_coherence[imin:imax] ) ) 

        if surrogate_test: 
            surrogates = surrogate_generation( y, nsurrogates = nsurrogates, method=surrogate_method, nsplits = surrogate_nsplits)
            for l, surrogate in enumerate(surrogates): 
                yf, y, _ = bandpass_filter_and_hilbert_transform(surrogate, norder=norder, fs=fs, f0=fg, df=10)
                
                fr, pyy  = scipy.signal.csd( y, y, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
                fr, pxy  = scipy.signal.csd( y, x, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap) 

                complex_coherence = pxy / np.sqrt( pxx*pyy ) 
                
                # compute psi ######################################################################
                cfd = []
                module = np.abs(complex_coherence)
                phase  = np.angle(complex_coherence)
                
                for i in range(imin, imax):  
                    factor1 = np.conj( complex_coherence[ i-ibeta//2: i+1+ibeta//2 ] )
                    factor2 = complex_coherence[ i+1-ibeta//2 : i+2+ibeta//2 ]
                    cfd.append(  np.imag( np.sum(factor1*factor2) ) )
                    
                    # factor1 = module[i-ibeta//2:i+ibeta//2]
                    # factor2 = module[i+1-ibeta//2:i+1+ibeta//2]
                    # factor3 = np.sin( phase[i+1-ibeta//2:i+1+ibeta//2]-phase[i-ibeta//2:i+ibeta//2])
                    # cfd.append( np.sum(factor1*factor2*factor3) ) 
        
                ####################################################################################
                psi_surrogates[l].append( np.array(cfd) )
                nu_surrogates[l].append( np.abs(complex_coherence[imin:imax]) )
    
    psi_surrogates = np.array(psi_surrogates) 
    nu_surrogates = np.array(nu_surrogates) 
    ftheta = fr[imin:imax]
    return ftheta, fgamma, psi, nu, psi_surrogates, nu_surrogates

#################################################################################################
    
def gamma_theta_phase_coupling( theta_phase, gamma_envelope, surrogate_test=True, nsurrogates=100, 
                              surrogate_method =  "block boostrapping"):
    
    height = []
    amplitude = [] 
    surrogates = []
    height_surrogates = []
    amplitude_surrogates = []
    np.random.seed(1)
    
    nbins = 20
    pbins = np.linspace(0,2*np.pi,nbins+1)
    for phi, amp in zip(theta_phase, gamma_envelope): 
        a_mean = np.zeros(len(pbins)-1)       
        p_mean = np.zeros(len(pbins)-1)       
        for k in range(len(pbins)-1):         
            pl = pbins[k]                     
            pr = pbins[k+1]                   
            indices=(phi>=pl) & (phi<pr)      
            a_mean[k] = np.mean(amp[indices])  
            p_mean[k] = np.mean([pL, pR]) 
            
        height.append( np.max(a_mean)-np.min(a_mean) )
        amplitude.append( a_mean )
        
        if surrogate_test:
            surrogates.append( [] )
            height_surrogates.append( [] )
            amplitude_surrogates.append( [] )

            surrogates = surrogate_generation( timeseries=amp, nsurrogates=nsurrogates, method = surrogate_method)
            
            for ii, surrogate in enumerate(surrogates):
                amean = np.zeros(len(pbins)-1)      
                pmean = np.zeros(len(pbins)-1)       
                for k in range(len(pbins)-1):         
                    pl = pbins[k]                    
                    pr = pbins[k+1]                   
                    indices=(phi>=pl) & (phi<pr)      
                    a_mean[k] = np.mean(surrogate[indices])  
                    p_mean[k] = np.mean([pL, pR]) 

                height_surrogates[-1].append( np.max(a_mean)-np.min(a_mean) )
                amplitude_surrogates[-1].append( a_mean )
                
    return amplitude, height, ampltiude_surrogates, height_surrogates, surrogates 

