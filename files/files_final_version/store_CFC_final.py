import sys
import os
import pandas as pd
import numpy as np 
sys.path.append('/home/dimitrios/Neurons/CA1model/final_files')
import file_management
import time as tm
from signal_analysis import compute_mi
from signal_analysis import compute_coherence_measurements

tick = tm.time()
#Calculate the avg modulation inxed over iseed and bseed
def get_avg_mi(fgamma,ftheta,fs,df,label,surrogate_test=False): 
    mi = np.zeros((len(fgamma),len(ftheta)))
    for input1 in np.unique(df["iseed"].values):
        for input2 in np.unique(df["ibseed"].values):
            print(input1,input2,len(df))
            ts = df[label][(df['iseed'] == input1) & (df['ibseed'] == input2)].values
            ts=ts[len(ts)//15:]
            ftheta, fgamma, mi_temp, mi_sur_temp = compute_mi(ts, ftheta, fgamma, fs=fs, 
                                          nbins=50, surrogate_test=surrogate_test, nsurrogates=10,choiseSurr = 2) 
            mi_temp = np.swapaxes(mi_temp,0,1)
            mi_sur_temp= np.swapaxes(mi_sur_temp,1,2)
            if(surrogate_test):
                inds_insig = mi_temp<=np.quantile(np.array(mi_sur_temp),0.995,axis=0)
                mi_temp[inds_insig] = 0
            mi=mi+mi_temp
    return mi/len(np.unique(df["iseed"].values))/len(np.unique(df["ibseed"].values))

current_folder = os.getcwd()
folder = os.path.join(current_folder, "baselineCA3")
ts_choise = ["lfp","ica","synapses_pyr"][1]
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
    elec_pos = np.unique(df["electrode"].values)[0]
    ts_ = df[ts_choise][(df["electrode"]==elec_pos)&(df["component"]==elec_pos)&(df["iseed"]==ised)&(df["ibseed"]==bsed)].values
    pos_vals = np.unique(df["electrode"].values)
elif(ts_choise=="synapses_pyr"):
    all_comp=["Bdend","soma","Adend1","Adend2","Adend3"]
    df= file_management.load_lzma(os.path.join(folder,ts_choise+"_"+area+".lzma"))
    df = pd.concat([*[df.filter(like=comp).mean(axis=1) for comp in all_comp],*[df["iseed"],df["ibseed"]]],axis=1,keys=[*all_comp,"iseed","ibseed"])
    ts_ = df[all_comp[0]].values
    df.fillna(0,inplace=True)
    pos_vals = all_comp

fs = 500
ftheta_min,ftheta_max = 4,20
fgamma_min,fgamma_max,dfgamma = 20,90,1
nperseg = 500
all_mi = {}

ftheta, fgamma, _, _, _, _ = compute_coherence_measurements(ts_, fs, nperseg, 0.75*nperseg, ftheta_min=ftheta_min, ftheta_max=ftheta_max, fgamma_min=fgamma_min, fgamma_max=fgamma_max, dfgamma=dfgamma,norder=4, nfft=2048,surrogate_test=False, nsurrogates=1,choiseSurr = None)
for ind,pos in enumerate(pos_vals):
    if(ts_choise=="lfp"):
        df_pos = df[df["electrode"]==pos]
    elif(ts_choise=="ica"):
        df_pos = df[(df["electrode"]==pos)&(df["component"]==pos)]
    elif(ts_choise=="synapses_pyr"):
        df_pos = df[[pos,"iseed","ibseed"]]
        df_pos.rename(columns={pos:ts_choise},inplace=True)
    all_mi[pos]= get_avg_mi(fgamma,ftheta,fs,df_pos,ts_choise,surrogate_test=False)

title = "CFC_"+ts_choise+"_"+area

file_management.save_lzma([all_mi,fgamma,ftheta],title,parent_dir=folder)

#file_name = __file__
#store = 'cp '+__file__.replace(current_folder+"/","") +" "+ folder.replace(current_folder+"/","") +__file__.replace(current_folder,"")

tock = tm.time()
diff = tock-tick
horas = int(diff/3600)
diff  = diff-horas*3600
mint  = int(diff/60)
diff  = diff-mint*60
seg   = int(diff)

print("All data saved after ", horas,'h', mint, 'min', seg, 's')
os.popen(store)

os.remove("CFC_out.dat")
os.remove("CFC_err.log")
