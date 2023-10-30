import numpy as np 
import matplotlib.pyplot as plt 
import os 
import sys 
sys.path.append('/home/dimitrios/Neurons/CA1model/final_files')

import file_management
import matplotlib.pyplot as plt
current_folder = os.getcwd()



def find_saturations(data_folder,area,thres=-40,window=500):
    seedlist=[]
    list_notBas=[]
    
    df_vbas = file_management.load_lzma(os.path.join(data_folder,"volt_bas_"+area+".lzma"))
    df_olm = file_management.load_lzma(os.path.join(data_folder,"volt_olm_"+area+".lzma"))
    df_pyr = file_management.load_lzma(os.path.join(data_folder,"volt_pyr_"+area+".lzma"))
    if(area=="ca1"):
        df_cck = file_management.load_lzma(os.path.join(data_folder,"volt_cck_"+area+".lzma"))
    #
    
    for iseed in np.unique(df_vbas["iseed"].values):
        df_vbas_temp = df_vbas[df_vbas["iseed"]==iseed]
        df_olm_temp = df_olm[df_olm["iseed"]==iseed]
        df_pyr_temp = df_pyr[df_pyr["iseed"]==iseed]
        if(area=="ca1"):
            df_cck_temp = df_cck[df_cck["iseed"]==iseed]
        
        for ibseed in np.unique(df_vbas_temp["ibseed"].values):
            vbas = df_vbas_temp[df_vbas_temp["ibseed"]==ibseed]["soma_volt_mean"].values
            volm = df_olm_temp[df_olm_temp["ibseed"]==ibseed]["soma_volt_mean"].values
            vpyr = df_pyr_temp[df_pyr_temp["ibseed"]==ibseed]["soma_volt_mean"].values
            flag_cck=False
            if(area=="ca1"):
                vcck = df_cck_temp[df_cck_temp["ibseed"]==ibseed]["soma_volt_mean"].values
                smoothed_vcck = np.convolve(vcck, np.ones(window) / window, mode='valid')
                flag_cck=np.sum(smoothed_vcck>thres)>0
            
            smoothed_bas = np.convolve(vbas, np.ones(window) / window, mode='valid')
            smoothed_olm = np.convolve(volm, np.ones(window) / window, mode='valid')
            smoothed_pyr = np.convolve(vpyr, np.ones(window) / window, mode='valid')
            
            if(np.sum(smoothed_bas>thres)>0 or np.sum(smoothed_olm>thres)>0 or np.sum(smoothed_pyr>thres)>0 or flag_cck):
                seedlist.append([iseed,ibseed])
            if(np.sum(smoothed_olm>thres)>0 or np.sum(smoothed_pyr>thres)>0 or flag_cck):
                list_notBas.append([iseed,ibseed])
            
    return seedlist,list_notBas

folders = ["baselineCA3","baselineCA1","folder52bCA3-0","folder52bCA3-1",
            "folder52bCA3-2","folder52bCA3-3","folder52bCA3-4"]
for f in folders:
    data_folder=os.path.join(os.getcwd(), f)
    if("CA3" in f):
        area = "ca3"
    else:
        area = "ca1"
    seedlist,list_notBas=find_saturations(data_folder,area)

    with open(os.path.join(data_folder,"saturation_lists.txt"), "w") as file:
        file.write("Saturated neuron seeds" + "\n")
        for item in seedlist:
            file.write(str(item) + "\n")
        file.write("Saturated NON bas neuron seeds" + "\n")
        for item in list_notBas:
            file.write(str(item) + "\n")