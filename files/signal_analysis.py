import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy 

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
    #print("Len std: ",len(sig),std)
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


def compute_mi(timeseries, ftheta, fgamma, fs=1e3, nbins=20, surrogate_test=False, nsurrogates=100,choiseSurr=0,Filter_theta=1.0,Filter_gamma=10.0): 
    ntheta,  ngamma = len(ftheta), len(fgamma)
    
    mi = np.zeros((ntheta,ngamma))
    mi_surrogates = np.zeros((nsurrogates, ntheta, ngamma))

    sig_band1 = [ [] for i in range(ntheta)]
    sig_band2 = [ [] for i in range(ngamma)]
    phase     = [ [] for i in range(ntheta)] 
    envelope  = [ [] for i in range(ngamma)]
    
    norder = 4
    for i,fth in enumerate(ftheta):
        sig_band1[i], _, phase[i] = bandpass_filter_and_hilbert_transform3(timeseries, fs=fs, f0=fth, df=Filter_theta, norder=norder)
    for i,fg in enumerate(fgamma):
        sig_band2[i], envelope[i], _ = bandpass_filter_and_hilbert_transform3(timeseries, fs=fs, f0=fg, df=Filter_gamma, norder=norder)

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
                if(choiseSurr==0):
                    surrogates = surrogate_generation( gamma_amplitude, nsurrogates = nsurrogates)
                elif(choiseSurr==1):
                    surrogates = surrogate_generation( gamma_amplitude, nsurrogates = nsurrogates,nsplits=int(len(x)/nperseg))
                elif(choiseSurr==2):
                    surrogates = surrogate_generation( gamma_amplitude, nsurrogates = nsurrogates,method="2block boostrapping")
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
from scipy import signal
def round_up_to_odd(f):
    f = int(np.ceil(f))
    return f + 1 if f % 2 == 0 else f
def eegfilt(data,srate,locutoff,hicutoff,epochframes=0,filtorder=0,revfilt=0,firtype="firls",causal=0):
    
    frames = len(data)
    nyq            = srate*0.5;  # Nyquist frequency
    MINFREQ = 0;

    minfac         = 3
    min_filtorder  = 15
    trans          = 0.15

    if locutoff > 0 and hicutoff > 0 and locutoff > hicutoff:
        raise ValueError('locutoff > hicutoff ???')

    if locutoff < 0 or hicutoff < 0:
        raise ValueError('locutoff | hicutoff < 0 ???')

    if locutoff > nyq:
        raise ValueError('Low cutoff frequency cannot be > srate/2')

    if hicutoff > nyq:
        raise ValueError('High cutoff frequency cannot be > srate/2')

    if filtorder == 0:
        if locutoff > 0:
            filtorder = minfac * int(srate / locutoff)
            filtorder=round_up_to_odd(filtorder)
        elif hicutoff > 0:
            filtorder = minfac * int(srate / hicutoff)
            filtorder=round_up_to_odd(filtorder)
        min_filtorder = 15  # Define min_filtorder if it's not defined previously
        if filtorder < min_filtorder:
            filtorder = min_filtorder
            filtorder=round_up_to_odd(filtorder)

    if epochframes == 0:
        epochframes = frames  # default

    epochs = frames // epochframes

    if epochs * epochframes != frames:
        raise ValueError('epochframes does not divide frames.')

    if filtorder * 3 > epochframes:  # MATLAB filtfilt() restriction
        print(f'eegfilt(): filter order is {filtorder}.')
        raise ValueError('epochframes must be at least 3 times the filtorder.')

    if (1 + trans) * hicutoff / nyq > 1:
        raise ValueError('high cutoff frequency too close to Nyquist frequency')

    if locutoff > 0 and hicutoff > 0:  # bandpass filter
        if revfilt:
            # print("eegfilt() - performing {}-point notch filtering.".format(filtorder))
            pass
        else:
            # print("eegfilt() - performing {}-point bandpass filtering.".format(filtorder))
            pass
        # print("If a message, 'Matrix is close to singular or badly scaled,' appears,")
        # print("then MATLAB has failed to design a good filter. As a workaround,")
        # print("for band-pass filtering, first highpass the data, then lowpass it.")

        if firtype == 'firls':
            f = [MINFREQ, (1 - trans) * locutoff / nyq, locutoff / nyq, hicutoff / nyq, (1 + trans) * hicutoff / nyq, 1]
            # print("eegfilt() - low transition band width is {:.1f} Hz; high trans. band width, {:.1f} Hz.".format((f[2] - f[1]) * srate / 2, (f[4] - f[3]) * srate / 2))
            m = [0, 0, 1, 1, 0, 0]
        elif firtype == 'fir1':
            from scipy.signal import fir1
            filtwts = fir1(filtorder, [locutoff, hicutoff] / (srate / 2))
    elif locutoff > 0:  # highpass filter
        if locutoff / nyq < MINFREQ:
            raise ValueError("eegfilt() - highpass cutoff freq must be > {} Hz".format(MINFREQ * nyq))
        # print("eegfilt() - performing {}-point highpass filtering.".format(filtorder))
        if firtype == 'firls':
            f = [MINFREQ, locutoff * (1 - trans) / nyq, locutoff / nyq, 1]
            # print("eegfilt() - highpass transition band width is {:.1f} Hz.".format((f[2] - f[1]) * srate / 2))
            m = [0, 0, 1, 1]
        elif firtype == 'fir1':
            from scipy.signal import firwin
            filtwts = firwin(filtorder + 1, locutoff / (srate / 2), pass_zero=False)
    elif hicutoff > 0:  # lowpass filter
        if hicutoff / nyq < MINFREQ:
            raise ValueError("eegfilt() - lowpass cutoff freq must be > {} Hz".format(MINFREQ * nyq))
        # print("eegfilt() - performing {}-point lowpass filtering.".format(filtorder))
        if firtype == 'firls':
            f = [MINFREQ, hicutoff / nyq, hicutoff * (1 + trans) / nyq, 1]
            # print("eegfilt() - lowpass transition band width is {:.1f} Hz.".format((f[2] - f[1]) * srate / 2))
            m = [1, 1, 0, 0]
        elif firtype == 'fir1':
            from scipy.signal import firwin
            filtwts = firwin(filtorder + 1, hicutoff / (srate / 2))
    else:
        raise ValueError('You must provide a non-0 low or high cut-off frequency')


    if revfilt:
        if firtype == 'fir1':
            raise ValueError("Cannot reverse filter using 'fir1' option")
        else:
            m = [not x for x in m]
    print(f,m)
    if firtype == 'firls':
        from scipy.signal import firls
        filtwts = firls(filtorder, f, m)

    smoothdata = np.zeros(frames)
    print(epochframes)
    for e in range(epochs):  # Filter each epoch, channel
        if causal:
            smoothdata[ (e) * epochframes:(e+1) * epochframes] = signal.lfilter(filtwts, 1, data[(e ) * epochframes:(e+1) * epochframes])
        else:
            smoothdata[ (e) * epochframes:(e+1) * epochframes] = signal.filtfilt(filtwts, 1, data[(e) * epochframes:(e+1) * epochframes])
    
    return smoothdata 

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

def bandpass_filter_and_hilbert_transform2(sig, fs, f0, df, norder):
    scale = sig.mean() 
    f_range = (f0-df, f0+df)
    sig_band = bandpass_filter(sig, norder, f_range, fs)
    phase, envelope= hilbert_phase(sig_band)
    return sig_band+scale, envelope,  phase

def bandpass_filter_and_hilbert_transform3(sig, fs, f0, df, norder):
    scale = sig.mean() 
    f_range = (f0-df, f0+df)
    print(f_range)
    sig_band = eegfilt(sig, fs, f_range[0],f_range[-1])
    phase, envelope= hilbert_phase(sig_band)
    
    return sig_band+scale, envelope,  phase

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
    
    if method == "2block boostrapping":
        for _ in range(nsurrogates):
            ir = np.random.randint(low=0,high=len(timeseries))
            surrogate = np.array( list(timeseries[ir:])+list(timeseries[:ir]) )
            surrogate_list.append( surrogate )
            
    if method == "sampling":
        for _ in range(nsurrogates):
            surrogate = np.random.choice(timeseries, size=n, replace=replace)
            surrogate_list.append( surrogate )
    
    if method == "fourier":
        for _ in range(nsurrogates):
            fourier = np.fft.fft(timeseries)
            random_phases = 1j*np.random.choice(fourier.imag[:n//2], size=n//2, replace=replace)
            random_phases = np.concatenate(( random_phases, random_phases[::-1] ))
            fourier_new = fourier*random_phases
            surrogate = np.real( np.fft.ifft(fourier_new))
            surrogate_list.append( surrogate )
            
    if method == "block boostrapping":
        timeseries_splitted = np.array_split(timeseries, nsplits) 
        for _ in range(nsurrogates):
            indexes_surrogate = np.random.choice( np.arange(nsplits), size = nsplits, replace=replace)
            surrogate = []
            for index in indexes_surrogate: 
                surrogate.append(  timeseries_splitted[index] )
            surrogate_list.append( np.concatenate(surrogate) )
            
    return surrogate_list 

def get_pvalue( samples,  popmean, method = "two-sided"):
    z0 = (popmean-samples.mean(0))/samples.std(0)
    if method == "two-sided":
        pvalue = scipy.stats.norm.sf(abs(z0))*2
    if method == "one-sided":
        pvalue = scipy.stats.norm.sf(abs(z0))
    return pvalue

#################################################################################################
def compute_coherence_measurements(x, fs, nperseg, noverlap, ftheta_min=4, ftheta_max=20, fgamma_min=25, fgamma_max=100, dfgamma=2,Filter_gamma=10,
                                   norder=4, nfft=2048,surrogate_test=True, nsurrogates=100,choiseSurr = 0):
    
    psi, nu = [],[] 
    psi_surrogates = [ [] for i in range(nsurrogates) ]
    nu_surrogates  = [ [] for i in range(nsurrogates) ]
    
    fgamma = np.arange(fgamma_min,fgamma_max,dfgamma)
    #scipy.signal.csd: cross power spectral density 
    fr, pxx  = scipy.signal.csd( x, x, fs=fs, window="hann", nfft=nfft, nperseg=nperseg, noverlap = noverlap)
    
    np.random.seed(12)
    for fg in fgamma:
        _, y, _ = bandpass_filter_and_hilbert_transform3(x, norder=norder, fs=fs, f0=fg, df=Filter_gamma)
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
            if(choiseSurr==0):
                surrogates = surrogate_generation( x, nsurrogates = nsurrogates)
            elif(choiseSurr==1):
                surrogates = surrogate_generation( x, nsurrogates = nsurrogates,nsplits=int(len(x)/nperseg))
            elif(choiseSurr==2):
                surrogates = surrogate_generation( x, nsurrogates = nsurrogates,method="2block boostrapping")
            for l, surrogate in enumerate(surrogates): 
                yf, y, _ = bandpass_filter_and_hilbert_transform3(surrogate, norder=norder, fs=fs, f0=fg, df=10)
                
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
                
    ftheta = fr[imin:imax]
    return ftheta, fgamma, psi, nu, psi_surrogates, nu_surrogates

#################################################################################################
    
def gamma_theta_phase_coupling( theta_phase, gamma_envelope, surrogate_test=True, nsurrogates=100 ):
    
    height = []
    amplitude = [] 
    amplitude_error = [] 
    surrogates = []
    height_surrogates = []
    amplitude_surrogates = []
    np.random.seed(1)
    
    nbins = 20
    pbins = np.linspace(0,2*np.pi,nbins+1)
    for phi, amp in zip(theta_phase, gamma_envelope): 
        a_mean = np.zeros(len(pbins)-1)       
        a_err  = np.zeros(len(pbins)-1)
        p_mean = np.zeros(len(pbins)-1)       
        for k in range(len(pbins)-1):         
            pl = pbins[k]                     
            pr = pbins[k+1]                   
            indices=(phi>=pl) & (phi<pr)      
            a_mean[k] = np.mean(amp[indices])  
            a_err[k]  = np.std(amp[indices])/np.sqrt(len(amp[indices]))
            p_mean[k] = np.mean([pl, pr]) 
            
        height.append( np.max(a_mean)-np.min(a_mean) )
        amplitude.append( a_mean )
        amplitude_error.append( a_err ) 
        
        
        if surrogate_test:
            surrogates.append( [] )
            height_surrogates.append( [] )
            amplitude_surrogates.append( [] )

            surrogates = surrogate_generation( timeseries=amp, nsurrogates=nsurrogates, method = "block boostrapping")
            
            for ii, surrogate in enumerate(surrogates):
                amean = np.zeros(len(pbins)-1)      
                pmean = np.zeros(len(pbins)-1)       
                for k in range(len(pbins)-1):         
                    pl = pbins[k]                    
                    pr = pbins[k+1]                   
                    indices=(phi>=pl) & (phi<pr)      
                    a_mean[k] = np.mean(surrogate[indices])  
                    p_mean[k] = np.mean([pl, pr]) 

                height_surrogates[-1].append( np.max(a_mean)-np.min(a_mean) )
                amplitude_surrogates[-1].append( a_mean )
                
    return amplitude, amplitude_error, p_mean, height, amplitude_surrogates, height_surrogates, surrogates 

