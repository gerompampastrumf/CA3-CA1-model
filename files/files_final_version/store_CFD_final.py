import sys
import os
import numpy as np 
sys.path.append('/home/dimitrios/Neurons/CA1model/final_files')
import file_management
import time as tm
from signal_analysis import compute_coherence_measurements
import pandas as pd
tick = tm.time()
#Calculate the avg modulation inxed over all_input2. Input3 and 4 are used in order to load the file.
def get_avg_cfd(fgamma,ftheta,fs,df,label,surrogate_test=False): 
    all_psi = []
    all_nu = []
    
    #for input1 in np.unique(df["iseed"].values):
    #    for input2 in np.unique(df["ibseed"].values):
    all_inputs2={}
    if(not surrogate_test):
        all_inputs1=np.unique(df["iseed"].values)
        for inp in all_inputs1:
            all_inputs2[inp]=np.unique(df[df["iseed"]==inp]["ibseed"].values)
    else:
        all_inputs1=np.unique(df["iseed"].values)[-6:]
        for inp in all_inputs1:
            all_inputs2[inp]=np.unique(df[df["iseed"]==inp]["ibseed"].values)[-6:]
    n_temps=0
    for input1 in all_inputs1:
        print(input1,all_inputs2[input1])
        for input2 in all_inputs2[input1]:
            print(input1,input2,len(df))
            n_temps+=1
            ts = df[label][(df['iseed'] == input1) & (df['ibseed'] == input2)].values
            ts=ts[len(ts)//15:]
            
            ftheta, fgamma, psi_temp, nu_temp, psi_sur_temp,nu_sur_temp = compute_coherence_measurements(ts, fs, nperseg, 0.75*nperseg, ftheta_min=ftheta_min, 
                                   ftheta_max=ftheta_max, fgamma_min=fgamma_min, fgamma_max=fgamma_max, dfgamma=dfgamma,
                                   norder=4, nfft=2048, surrogate_test=surrogate_test, nsurrogates=100,choiseSurr = 2)
            psi_temp= np.array(psi_temp)
            nu_temp=np.array(nu_temp)

            if(surrogate_test):
                inds_insig = np.logical_and(psi_temp>=np.quantile(np.array(psi_sur_temp),0.005,axis=0),psi_temp<=np.quantile(np.array(psi_sur_temp),0.995,axis=0))
                psi_temp[inds_insig]=0
                inds_insig = nu_temp<=np.quantile(np.array(nu_sur_temp),0.99,axis=0)
                nu_temp[inds_insig]=0
        
            all_psi.append(psi_temp)
            all_nu.append(nu_temp)
    print(n_temps)
    return all_psi,all_nu

current_folder = os.getcwd()
folder = os.path.join(current_folder, "lessDGburstCA1-12")
input1 = int(sys.argv[1])
ts_choise = ["lfp","ica","synapses_pyr"][input1]
if "CA1" in folder.split('/')[-1]:
    area="ca1"
else:
    area="ca3"
surrogate_test=False


df= file_management.load_lzma(os.path.join(folder,ts_choise+"_"+area+".lzma"))
ised = np.unique(df["iseed"].values)[0]
bsed = np.unique(df["ibseed"].values)[0]
if(ts_choise=="lfp"):
    elec_pos = np.unique(df["electrode"].values)[0]
    ts_ = df[ts_choise][(df["electrode"]==elec_pos)&(df["iseed"]==ised)&(df["ibseed"]==bsed)].values
    pos_vals = np.unique(df["electrode"].values)
elif(ts_choise=="ica"):
    elec_comp=np.unique(df["component"].values)[0]
    elec_pos = np.unique(df["electrode"].values)[0]
    ts_ = df[ts_choise][(df["electrode"]==elec_pos)&(df["component"]==elec_comp)&(df["iseed"]==ised)&(df["ibseed"]==bsed)].values
    pos_vals = ["Bdend","soma","Adend1","Adend2","Adend3"]
elif(ts_choise=="synapses_pyr"):
    all_comp=["Bdend","soma","Adend1","Adend2","Adend3"]
    df= file_management.load_lzma(os.path.join(folder,ts_choise+"_"+area+".lzma"))
    df = pd.concat([*[df.filter(like=comp).filter(like="mean").sum(axis=1) for comp in all_comp],*[df["iseed"],df["ibseed"]]],axis=1,keys=[*all_comp,"iseed","ibseed"])
    ts_ = df[all_comp[0]].values
    df.fillna(0,inplace=True)
    pos_vals = all_comp

fs = 500
ftheta_min,ftheta_max = 2,15
fgamma_min,fgamma_max,dfgamma = 20,90,1
nperseg = 500

all_psi= {}
all_psi_std={}
all_nu = {}
all_nu_std={}

ftheta, fgamma, _, _, _, _ = compute_coherence_measurements(ts_, fs, nperseg, 0.75*nperseg, ftheta_min=ftheta_min, ftheta_max=ftheta_max, fgamma_min=fgamma_min, fgamma_max=fgamma_max, dfgamma=dfgamma,norder=4, nfft=2048,surrogate_test=False, nsurrogates=1,choiseSurr = None)
for ind,pos in enumerate(pos_vals):
    if(ts_choise=="lfp"):
        df_pos = df[df["electrode"]==pos]
    elif(ts_choise=="ica"):
        df_pos = df[(df["electrode"]==elec_pos)&(df["component"]==pos)]
    elif(ts_choise=="synapses_pyr"):
        df_pos = df[[pos,"iseed","ibseed"]]
        df_pos.rename(columns={pos:ts_choise},inplace=True)

    psi_tmp,nu_tmp= get_avg_cfd(fgamma,ftheta,fs,df_pos,ts_choise,surrogate_test=surrogate_test)
    all_psi[pos]=np.mean(psi_tmp,axis=0)
    all_psi_std[pos]=np.std(psi_tmp,axis=0)
    all_nu[pos]=np.mean(nu_tmp,axis=0)
    all_nu_std[pos]=np.std(nu_tmp,axis=0)

if(surrogate_test):
    title = "CFD_sign_test"+ts_choise+"_"+area
else:
    title = "CFD_"+ts_choise+"_"+area

file_management.save_lzma([all_psi,all_psi_std,all_nu,all_nu_std,fgamma,ftheta],title,parent_dir=folder)

file_name = __file__
store = 'cp '+__file__.replace(current_folder+"/","") +" "+ folder.replace(current_folder+"/","") +__file__.replace(current_folder,"")

tock = tm.time()
diff = tock-tick
horas = int(diff/3600)
diff  = diff-horas*3600
mint  = int(diff/60)
diff  = diff-mint*60
seg   = int(diff)

print("All data saved after ", horas,'h', mint, 'min', seg, 's')
os.popen(store)

os.remove("CFD_out.dat")
os.remove("CFD_err.log")