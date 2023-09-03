import os 
import sys
import numpy as np 
sys.path.append('/home/jaime/Desktop/hippocampus/files_final_version/')
import file_management


def spiking_phases(folder,region,ax):
    data_folder = os.path.join(folder,"measurements")
    filename = os.path.join(data_folder,f"phases_spikes_distribution_{region}.lzma")
    data = file_management.load_lzma(filename)
    if region == "ca3":
        labels = ["pyr","bas","olm"]
    elif region == "ca1":
        labels = ["pyr","bas","olm","cck"]

    binphase = np.arange(0,360,15)
    binphase = binphase[1:]-np.diff(binphase)[0]/2
    bins = binphase.tolist() + (binphase+360).tolist() + [binphase[0]+720]

    for label in labels: 
        dens = data[f"{label}_mean"].values
        dens = dens.tolist() + dens.tolist() + [dens[0]]
        try:
            ax.step(bins,dens,label=label,where="mid")
            ax.fill_between(bins,dens,step="mid",alpha=0.25)
        except:
            print(f"{label} has not spikes")

def psd_voltage(folder,region,ax):
    data_folder = os.path.join(folder,"measurements")
    filename = os.path.join(data_folder,f"psd_voltage_{region}.lzma")
    data = file_management.load_lzma(filename)
    freq = data["freq"].values
    df = np.diff(freq)[0]
    w = (freq<100) & (freq>1)
    if region == "ca3": 
        labels = ["pyr","bas","olm"]
        lbs    = ["soma","bas","olm"]
    elif region == "ca1":
        labels = ["pyr","bas","olm","cck"]
        lbs    = ["soma","bas","olm","cck"]

    for label,lb in zip(labels,lbs):
        mu, sem = data[f"{lb}_mean"].values[w], data[f"{lb}_sem"].values[w]
        norm = np.sum(mu*df)
        mu /= norm
        sem /= norm
        ax.plot(freq[w],mu,label=label)
        ax.fill_between(freq[w],mu-sem,mu+sem,alpha=0.25)

def theta_gamma_coupling(folder, series, region, theta_component, gamma_component, ax,axt):
    data_folder = os.path.join(folder,"measurements")
    filename = os.path.join(data_folder,f"{series}_{region}_amplitude_coupling.lzma")

    if series == "synapses_pyr":
        ''' cross gamma-theta coupling not implemented.
        only the auto coupling is implemented '''
        data = file_management.load_lzma(filename)[f"i{theta_component}"]
        bins = data["bins"].values
        theta_mean, theta_sem = data["theta_amplitude_mean"].values, data["theta_amplitude_sem"].values
        gamma_mean, gamma_sem = data["gamma_amplitude_mean"].values, data["gamma_amplitude_sem"].values

        ax.step(bins, gamma_mean, label="gamma",where="mid")
        #ax.fill_between(bins, gamma_mean-gamma_sem, gamma_mean+gamma_sem, alpha=0.25,)
        ax.fill_between(bins, gamma_mean, alpha=0.25,step="mid")
        axt.plot(bins, theta_mean, label="theta", color="black")
        axt.fill_between(bins, theta_mean-theta_sem, theta_mean+theta_sem, alpha=0.25, color="black")
    else:
        data = file_management.load_lzma(filename)[f"{theta_component}_{gamma_component}"]
        bins = data["bins"].values
        theta_mean, theta_sem = data["theta_amplitude_mean"].values, data["theta_amplitude_sem"].values
        gamma_mean, gamma_sem = data["gamma_amplitude_mean"].values, data["gamma_amplitude_sem"].values

        ax.step(bins, gamma_mean, label="gamma",where="mid")
        #ax.fill_between(bins, gamma_mean-gamma_sem, gamma_mean+gamma_sem, alpha=0.25,)
        ax.fill_between(bins, gamma_mean, alpha=0.25,step="mid")
        axt.plot(bins, theta_mean, label="theta", color="black")
        axt.fill_between(bins, theta_mean-theta_sem, theta_mean+theta_sem, alpha=0.25, color="black")

def raster_plot(folder,region,ax,t1=10000,t2=12000):
    # spikes data is not in meausrements path
    if region == "ca3":
        files = ["external_inputs_pyr_ca3.lzma","spikes_bas_ca3.lzma","spikes_olm_ca3.lzma"]
        labels = ["pyr","bas","olm"]
    elif region == "ca1":
        files = ["spikes_pyr_ca1.lzma","spikes_bas_ca1.lzma","spikes_olm_ca1.lzma","spikes_cck_ca1.lzma"]
        labels = ["pyr","bas","olm","cck"]

    for file,label in zip(files,labels): 
        filename = os.path.join(folder,file)
        data = file_management.load_lzma(filename)
        w = (data["iseed"]==4) & (data["ibseed"]==4) & (data["tvec"]>t1) & (data["tvec"]<t2) 
        ax.plot(data["tvec"][w],data["idvec"][w],"o",markersize=1,label=label)


def gain_evolution(external_input,comp,axs):
    ''' evolution of the gamma and theta coherence as a function of the delay in the 
    different time series: membrare potential (volt), lfp, and ica'''

    folder = '/home/jaime/Desktop/hippocampus/external_inputs/baseline/'

    folders = [ os.path.join(folder,f"{external_input}CA3_{i}") for i in range(4)]
    folders.append(folder)
    gain = np.arange(0,1.1,0.25)

    filename = " psd_voltage_measurements_ca3.lzma"
    theta_power_mean = np.zeros(5)
    theta_power_sem  = np.zeros(5)
    gamma_power_mean = np.zeros(5)
    gamma_power_sem  = np.zeros(5)

    for i,folder in enumerate(folders):
        path = os.path.join(folder,"measurements")
        data = file_management.load_lzma(os.path.join(path,filename))

        theta_power_mean[i] = np.mean(data[f"{comp}"]["theta_power"].values)
        theta_power_sem[i]  = np.std(data[f"{comp}"]["theta_power"].values)/np.sqrt(len(data[f"{comp}"].values))
        gamma_power_mean[i] = np.mean(data[f"{comp}"]["gamma_power"].values)
        gamma_power_sem[i]  = np.std(data[f"{comp}"]["gamma_power"].values)/np.sqrt(len(data[f"{comp}"].values))

    theta_power_mean /= theta_power_mean[-1]
    theta_power_sem  /= theta_power_mean[-1]
    gamma_power_mean /= gamma_power_mean[-1]
    gamma_power_sem  /= gamma_power_mean[-1]

    axs[0].errorbar(gain,theta_power_mean,yerr=theta_power_sem,label="theta")
    axs[0].errorbar(gain,gamma_power_mean,yerr=gamma_power_sem,label="gamma")

    filename = "psd_lfp_measurements_ca3.lzma"

    for i,folder in enumerate(folders):
        path = os.path.join(folder,"measurements")
        data = file_management.load_lzma(os.path.join(path,filename))

        theta_power_mean[i] = np.mean(data[f"{comp}"]["theta_power"].values)
        theta_power_sem[i]  = np.std(data[f"{comp}"]["theta_power"].values)/np.sqrt(len(data[f"{comp}"].values))
        gamma_power_mean[i] = np.mean(data[f"{comp}"]["gamma_power"].values)
        gamma_power_sem[i]  = np.std(data[f"{comp}"]["gamma_power"].values)/np.sqrt(len(data[f"{comp}"].values))


    theta_power_mean /= theta_power_mean[-1]
    theta_power_sem  /= theta_power_mean[-1]
    gamma_power_mean /= gamma_power_mean[-1]
    gamma_power_sem  /= gamma_power_mean[-1]

    axs[1].errorbar(gain,theta_power_mean,yerr=theta_power_sem,label="theta")
    axs[1].errorbar(gain,gamma_power_mean,yerr=gamma_power_sem,label="gamma")

    filename = "psd_ica_measurements_ca3.lzma"

    for i,folder in enumerate(folders):
        path = os.path.join(folder,"measurements")
        data = file_management.load_lzma(os.path.join(path,filename))

        theta_power_mean[i] = np.mean(data[f"{comp}_{comp}"]["theta_power"].values)
        theta_power_sem[i]  = np.std(data[f"{comp}_{comp}"]["theta_power"].values)/np.sqrt(len(data[f"{comp}_{comp}"].values))
        gamma_power_mean[i] = np.mean(data[f"{comp}_{comp}"]["gamma_power"].values)
        gamma_power_sem[i]  = np.std(data[f"{comp}_{comp}"]["gamma_power"].values)/np.sqrt(len(data[f"{comp}_{comp}"].values))

    theta_power_mean /= theta_power_mean[-1]
    theta_power_sem  /= theta_power_mean[-1]
    gamma_power_mean /= gamma_power_mean[-1]
    gamma_power_sem  /= gamma_power_mean[-1]

    axs[2].errorbar(gain,theta_power_mean,yerr=theta_power_sem,label="theta")
    axs[2].errorbar(gain,gamma_power_mean,yerr=gamma_power_sem,label="gamma")

    filename = "psd_synaptic_measurements_ca3.lzma"

    for i,folder in enumerate(folders):
        path = os.path.join(folder,"measurements")
        data = file_management.load_lzma(os.path.join(path,filename))

        theta_power_mean[i] = np.mean(data[f"i{comp}"]["theta_power"].values)
        theta_power_sem[i]  = np.std(data[f"i{comp}"]["theta_power"].values)/np.sqrt(len(data[f"i{comp}"].values))
        gamma_power_mean[i] = np.mean(data[f"i{comp}"]["gamma_power"].values)
        gamma_power_sem[i]  = np.std(data[f"i{comp}"]["gamma_power"].values)/np.sqrt(len(data[f"i{comp}"].values))

    theta_power_mean /= theta_power_mean[-1]
    theta_power_sem  /= theta_power_mean[-1]
    gamma_power_mean /= gamma_power_mean[-1]
    gamma_power_sem  /= gamma_power_mean[-1]

    axs[3].errorbar(gain,theta_power_mean,yerr=theta_power_sem,label="theta")
    axs[3].errorbar(gain,gamma_power_mean,yerr=gamma_power_sem,label="gamma")
    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)
    axs[3].grid(True)


def main_phase_gain(folders, region, ax):
    gain = np.arange(0,1.1,0.25)
    phase_max = [[] for i in range(len(folders))]
    if region == "ca3":
        labels = ["pyr","bas","olm"]
    elif region == "ca1":
        labels = ["pyr","bas","olm","cck"]
    
    for i,folder in enumerate(folders):
        path = os.path.join(folder,"measurements")
        data = file_management.load_lzma(os.path.join(path,f"phases_spikes_distribution_{region}.lzma"))
        binphase = data["binphase"].values

        for j,label in enumerate(labels):
            counts = data[f"{label}_mean"].values
            phase_max[i].append( binphase[np.argmax(counts)] )
        
    phase_max = np.array(phase_max)
    for i,label in enumerate(labels):
        ax.plot(gain,phase_max[:,i],'o-',label=label)
        ax.set_ylim([0,360])
        ax.grid(True)

# main phase theta-gamma coupling 
def coupling_phase_gain(folders,series,region,ax):
    gain = np.arange(0,1.1,0.25)
    phase_max = [[] for i in range(len(folders))]
    amp_max   = [[] for i in range(len(folders))]
    labels = ["Bdend","soma","Adend1","Adend2","Adend3"]

    for i, folder in enumerate(folders):
        path = os.path.join(folder,"measurements")
        data = file_management.load_lzma(os.path.join(path,f"{series}_{region}_amplitude_coupling.lzma"))

        bins = data["soma_soma"]["bins"].values
        for j,label in enumerate(labels):

            amp = data[f"{label}_{label}"]["gamma_amplitude_mean"].values
            phase_max[i].append( bins[np.argmax(amp)] )
            amp_max[i].append( np.max(amp)-np.min(amp) )
        
    phase_max = np.array(phase_max)
    amp_max   = np.array(amp_max)
    for i,label in enumerate(labels):
        ax[0].plot(gain,phase_max[:,i]*180/np.pi,'o-',label=label)
        ax[1].plot(gain,amp_max[:,i]*180/np.pi,'o-',label=label)
    ax[0].grid(True)
    ax[1].grid(True)

def frequency_gain(folders,series,region,ax):
    gain = np.arange(0,1.1,0.25)
    gamma_freq_max = [[] for i in range(len(folders))]
    gamma_freq_sem = [[] for i in range(len(folders))]
    theta_freq_max = [[] for i in range(len(folders))]
    theta_freq_sem = [[] for i in range(len(folders))]
    gamma_amp_max = [[] for i in range(len(folders))]
    theta_amp_max = [[] for i in range(len(folders))]

    labels = ["Bdend","soma","Adend1","Adend2","Adend3"]
    for i, folder in enumerate(folders):
        path = os.path.join(folder,"measurements")
        data = file_management.load_lzma(os.path.join(path,f"psd_{series}_measurements_{region}.lzma"))

        if series == "ica":
            for j,label in enumerate(labels):
                gamma_freq_max[i].append( np.mean(data[f"{label}_{label}"]["gamma_freq"]) )
                gamma_freq_sem[i].append( np.std(data[f"{label}_{label}"]["gamma_freq"])/np.sqrt(len(data[f"{label}_{label}"]["gamma_freq"])) )
                theta_freq_max[i].append( np.mean(data[f"{label}_{label}"]["theta_freq"]) )
                theta_freq_sem[i].append( np.std(data[f"{label}_{label}"]["theta_freq"])/np.sqrt(len(data[f"{label}_{label}"]["theta_freq"])) )
        else:
            for j,label in enumerate(labels):
                gamma_freq_max[i].append( np.mean(data[f"{label}"]["gamma_freq"]) )
                gamma_freq_sem[i].append( np.std(data[f"{label}"]["gamma_freq"])/np.sqrt(len(data[f"{label}"]["gamma_freq"])) )
                theta_freq_max[i].append( np.mean(data[f"{label}"]["theta_freq"]) )
                theta_freq_sem[i].append( np.std(data[f"{label}"]["theta_freq"])/np.sqrt(len(data[f"{label}"]["theta_freq"])) )
        
    gamma_freq_max = np.array(gamma_freq_max)
    gamma_freq_sem = np.array(gamma_freq_sem)
    theta_freq_max = np.array(theta_freq_max)
    theta_freq_sem = np.array(theta_freq_sem)

    for i,label in enumerate(labels):
        ax[0].errorbar(gain,gamma_freq_max[:,i],yerr=gamma_freq_sem[:,i],label=label)
        ax[1].errorbar(gain,theta_freq_max[:,i],yerr=theta_freq_sem[:,i],label=label)
    ax[0].grid(True)
    ax[1].grid(True)


# new 
def theta_gamma_coupling_same_ref(folder, series, region, theta_component, gamma_component, ax,axt):
    '''
    Gamma and theta coupling are computed using the same reference, also interegionally 
    '''
    data_folder = os.path.join(folder,"measurements")
    filename = os.path.join(data_folder,f"{series}_{region}_amplitude_coupling.lzma")

    if series == "synapses_pyr":
        ''' cross gamma-theta coupling not implemented.
        only the auto coupling is implemented '''
        data = file_management.load_lzma(filename)[f"i{theta_component}"]
        bins = data["bins"].values
        theta_mean, theta_sem = data["theta_amplitude_mean"].values, data["theta_amplitude_sem"].values
        gamma_mean, gamma_sem = data["gamma_amplitude_mean"].values, data["gamma_amplitude_sem"].values

        ax.step(bins, gamma_mean, label="gamma",where="mid")
        #ax.fill_between(bins, gamma_mean-gamma_sem, gamma_mean+gamma_sem, alpha=0.25,)
        ax.fill_between(bins, gamma_mean, alpha=0.25,step="mid")
        axt.plot(bins, theta_mean, label="theta", color="black")
        axt.fill_between(bins, theta_mean-theta_sem, theta_mean+theta_sem, alpha=0.25, color="black")
    else:
        data = file_management.load_lzma(filename)[f"{theta_component}_{gamma_component}"]
        bins = data["bins"].values
        theta_mean, theta_sem = data["theta_amplitude_mean"].values, data["theta_amplitude_sem"].values
        gamma_mean, gamma_sem = data["gamma_amplitude_mean"].values, data["gamma_amplitude_sem"].values

        ax.step(bins, gamma_mean, label="gamma",where="mid")
        #ax.fill_between(bins, gamma_mean-gamma_sem, gamma_mean+gamma_sem, alpha=0.25,)
        ax.fill_between(bins, gamma_mean, alpha=0.25,step="mid")
        axt.plot(bins, theta_mean, label="theta", color="black")
        axt.fill_between(bins, theta_mean-theta_sem, theta_mean+theta_sem, alpha=0.25, color="black")