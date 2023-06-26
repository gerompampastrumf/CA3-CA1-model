'''#############################################################################
Date: 24/11/2022

version where I start removing everything related with basket and olm cells (by now)
in ca1. The rest is ok to keep
#############################################################################'''


''' ############################################################################
            Parameters: neurons -> neurons
#############################################################################'''
weights_neurons_neurons = {}
nsyns_neurons_neurons = {}
syn_neurons_neurons  = {}
delay_neurons_neurons = {}

scale = 0.87 # to compensate the fast part of NMDA that I do not understand
# Internal CA3
#weights_neurons_neurons["pyr_ca3_to_pyr_ca3"] = [[0.02348e-3, 0.004e-3]] # [0.02e-3+scale*0.004e-3, 0.004e-3]
#weights_neurons_neurons["bas_ca3_to_pyr_ca3"] = [[0.576e-3]] # [0.8*0.72e-3]
#weights_neurons_neurons["olm_ca3_to_pyr_ca3"] = [[57.6e-3]]  # [0.8*72.0e-3] #*0.9
#weights_neurons_neurons["pyr_ca3_to_bas_ca3"] = [[1.5606e-3, 1.38e-3]] #[0.36e-3+scale*1.38e-3, 1.38e-3]
#weights_neurons_neurons["bas_ca3_to_bas_ca3"] = [[4.5e-3]]
#eights_neurons_neurons["pyr_ca3_to_olm_ca3"] = [[0.969e-3, 0.7e-3]] #[0.36e-3+scale*0.7e-3, 0.7e-3]
# Internal CA1
#weights_neurons_neurons["bas_ca1_to_pyr_ca1"] = [[0.576e-3]] # [0.8*0.72e-3]
#weights_neurons_neurons["olm_ca1_to_pyr_ca1"] = [[57.6e-3]]  # [0.8*72.0e-3] #*0.9
#weights_neurons_neurons["pyr_ca1_to_bas_ca1"] = [[1.5606e-3, 1.38e-3]] # [0.36e-3+scale*1.38e-3, 1.38e-3]
#weights_neurons_neurons["bas_ca1_to_bas_ca1"] = [[4.5e-3]]
#weights_neurons_neurons["pyr_ca1_to_olm_ca1"] = [[0.969e-3, 0.7e-3]] # [0.36e-3+scale*0.7e-3, 0.7e-3]
#weights_neurons_neurons["cck_ca1_to_pyr_ca1"] = [[0e-3, 0e-3]] # to be optimize # updated to

# Schaffer colateral
# weights_neurons_neurons["pyr_ca3_to_pyr_ca1"] = [[0e-3, 0e-3]] # to be optimize
# weights_neurons_neurons["pyr_ca3_to_bas_ca1"] = [[0e-3, 0e-3]]# to be optimize
# weights_neurons_neurons["pyr_ca3_to_cck_ca1"] = [[0e-3, 0e-3]] # to be optimize

# Internal CA3
nsyns_neurons_neurons["pyr_ca3_to_pyr_ca3"] = [25]
nsyns_neurons_neurons["bas_ca3_to_pyr_ca3"] = [42]
nsyns_neurons_neurons["olm_ca3_to_pyr_ca3"] = [10]
nsyns_neurons_neurons["pyr_ca3_to_bas_ca3"] = [100]
nsyns_neurons_neurons["bas_ca3_to_bas_ca3"] = [60]
nsyns_neurons_neurons["pyr_ca3_to_olm_ca3"] = [10]
# Internal CA1
#nsyns_neurons_neurons["bas_ca1_to_pyr_ca1"] = [42]
#nsyns_neurons_neurons["olm_ca1_to_pyr_ca1"] = [10]
#nsyns_neurons_neurons["pyr_ca1_to_bas_ca1"] = [100]
#nsyns_neurons_neurons["bas_ca1_to_bas_ca1"] = [60]
#nsyns_neurons_neurons["pyr_ca1_to_olm_ca1"] = [10]
#nsyns_neurons_neurons["cck_ca1_to_pyr_ca1"] = [25,25] # to optimized # esto hay que tocarlo, habria que poner [ int(0.7*n), int(0.3*n)]

# Internal CA3
syn_neurons_neurons["pyr_ca3_to_pyr_ca3"] = [["BdendAMPA_pyr","BdendNMDA_pyr"]]  # [["BdendAMPAf","BdendNMDA"]]
syn_neurons_neurons["bas_ca3_to_pyr_ca3"] = [["somaGABA_bas"]]                   # [["somaGABA"]]
syn_neurons_neurons["olm_ca3_to_pyr_ca3"] = [["Adend3GABA_olm"]]                 # [["Adend3GABAfb"]]
syn_neurons_neurons["pyr_ca3_to_bas_ca3"] = [["somaAMPA_pyr","somaNMDA_pyr"]]    # [["somaAMPAf", "somaNMDA"]]
syn_neurons_neurons["bas_ca3_to_bas_ca3"] = [["somaGABA_bas"]]                   # [["somaGABAf"]]
syn_neurons_neurons["pyr_ca3_to_olm_ca3"] = [["somaAMPA_pyr","somaNMDA_pyr"]]     # [["somaAMPAf","somaNMDA"]]
# Internal CA1
#syn_neurons_neurons["bas_ca1_to_pyr_ca1"] = [["somaGABA_pyr"]]                   # [["somaGABAf"]]
#syn_neurons_neurons["olm_ca1_to_pyr_ca1"] = [["Adend3GABA_olm"]]                 # [["Adend3GABAfb"]]
#syn_neurons_neurons["pyr_ca1_to_bas_ca1"] = [["somaAMPA_bas","somaNMDA_bas"]]    # [["somaAMPAf", "somaNMDA"]]
#syn_neurons_neurons["bas_ca1_to_bas_ca1"] = [["somaGABA_bas"]]                   # [["somaGABAf"]]
#syn_neurons_neurons["pyr_ca1_to_olm_ca1"] = [["somaAMPA_pyr","somaNMDA_pyr"]]    # [["somaAMPAf","somaNMDA"]]
#syn_neurons_neurons["cck_ca1_to_pyr_ca1"] = [["Adend2GABA_cck", "somaGABA_cck"]] # [["Adend2GABAf", "somaGABAf"]] # updated on
# Schaffer colaterals
# syn_neurons_neurons["pyr_ca3_to_pyr_ca1"] = [["Adend1AMPAf","Adend1NMDA"]]
# syn_neurons_neurons["pyr_ca3_to_bas_ca1"] = [["somaAMPAf","somaNMDA"]]
# syn_neurons_neurons["pyr_ca3_to_cck_ca1"] = [["somaAMPAf","somaNMDA"]]

# Internal CA3 # Added in 15/11
delay_neurons_neurons["pyr_ca3_to_pyr_ca3"] = [2.0,1.0]
delay_neurons_neurons["bas_ca3_to_pyr_ca3"] = [2.0,1.0]
delay_neurons_neurons["olm_ca3_to_pyr_ca3"] = [2.0,1.0]
delay_neurons_neurons["pyr_ca3_to_bas_ca3"] = [2.0,1.0]
delay_neurons_neurons["bas_ca3_to_bas_ca3"] = [2.0,1.0]
delay_neurons_neurons["pyr_ca3_to_olm_ca3"] = [2.0,1.0]

# Internal CA1
#delay_neurons_neurons["bas_ca1_to_pyr_ca1"] = [2.0,1.0]
#delay_neurons_neurons["olm_ca1_to_pyr_ca1"] = [2.0,1.0]
#delay_neurons_neurons["pyr_ca1_to_bas_ca1"] = [2.0,1.0]
#delay_neurons_neurons["pyr_ca1_to_cck_ca1"] = [2.0,1.0]
#delay_neurons_neurons["bas_ca1_to_bas_ca1"] = [2.0,1.0]
#delay_neurons_neurons["pyr_ca1_to_olm_ca1"] = [2.0,1.0]
#delay_neurons_neurons["cck_ca1_to_pyr_ca1"] = [2.0,1.0]
# Schaffer colaterals
# syn_neurons_neurons["pyr_ca3_to_pyr_ca1"] = [2.0]
# syn_neurons_neurons["pyr_ca3_to_bas_ca1"] = [2.0]
# syn_neurons_neurons["pyr_ca3_to_cck_ca1"] = [2.0]

''' ############################################################################
            Parameters: inputs -> neurons
#############################################################################'''
weights_inputs_neurons = {}
nsyns_inputs_neurons   = {}
syn_inputs_neurons     = {}
delay_inputs_neurons   = {}

# to CA3
#weights_inputs_neurons["sep_180_to_bas_ca3"] = [[6.4e-4]]
#weights_inputs_neurons["sep_360_to_olm_ca3"] = [[3.2e-4]]
#weights_inputs_neurons["ec2_180_to_pyr_ca3"] = [[3.3e-4, 1.8e-4]] # AMPA/ NMDA
#weights_inputs_neurons["ec2_180_to_bas_ca3"] = [[3.3e-4, 1.8e-4]]
#weights_inputs_neurons["ec2_360_to_pyr_ca3"] = [[0.0,0.0]]
#weights_inputs_neurons["ec2_360_to_bas_ca3"] = [[0.0,0.0]]
#weights_inputs_neurons["dg_regular_to_pyr_ca3"] = [[0.9e-4,0.5e-4]]
#weights_inputs_neurons["dg_regular_to_bas_ca3"] = [[0.9e-4,0.5e-4]]
#weights_inputs_neurons["dg_burst_to_pyr_ca3"] = [[0.2*0.9e-4, 0.2*0.5e-4]]
#weights_inputs_neurons["dg_burst_to_bas_ca3"] = [[0.2*0.9e-4, 0.2*0.5e-4]]
# to CA1
#weights_inputs_neurons["sep_180_to_bas_ca1"] = [[6.4e-4]]
#weights_inputs_neurons["sep_360_to_olm_ca1"] = [[3.2e-4]]
#weights_inputs_neurons["ec3_180_to_pyr_ca1"] = [[0.0,0.0]]
#weights_inputs_neurons["ec3_180_to_bas_ca1"] = [[0.0,0.0]]
#weights_inputs_neurons["ec3_180_to_cck_ca1"] = [[0.0,0.0]]
#weights_inputs_neurons["ec3_360_to_pyr_ca1"] = [[3.3e-4, 1.8e-4]]
#weights_inputs_neurons["ec3_360_to_bas_ca1"] = [[3.3e-4, 1.8e-4]]
#weights_inputs_neurons["ec3_360_to_cck_ca1"] = [[3.3e-4, 1.8e-4]]

# to CA3
nsyns_inputs_neurons["sep_180_to_bas_ca3"] = [10]
nsyns_inputs_neurons["sep_360_to_olm_ca3"] = [10]
nsyns_inputs_neurons["ec2_180_to_pyr_ca3"] = [25]
nsyns_inputs_neurons["ec2_180_to_bas_ca3"] = [25]
nsyns_inputs_neurons["ec2_360_to_pyr_ca3"] = [25]
nsyns_inputs_neurons["ec2_360_to_bas_ca3"] = [25]
nsyns_inputs_neurons["dg_regular_to_pyr_ca3"] = [25]
nsyns_inputs_neurons["dg_regular_to_bas_ca3"] = [25]
nsyns_inputs_neurons["dg_burst_to_pyr_ca3"] = [10]
nsyns_inputs_neurons["dg_burst_to_bas_ca3"] = [10]
# to CA1
#nsyns_inputs_neurons["sep_180_to_bas_ca1"] = [10]
#nsyns_inputs_neurons["sep_360_to_olm_ca1"] = [10]
#nsyns_inputs_neurons["ec3_180_to_pyr_ca1"] = [25]
#nsyns_inputs_neurons["ec3_180_to_bas_ca1"] = [25]
#nsyns_inputs_neurons["ec3_180_to_cck_ca1"] = [25]
#nsyns_inputs_neurons["ec3_360_to_pyr_ca1"] = [25]
#nsyns_inputs_neurons["ec3_360_to_bas_ca1"] = [25]
#nsyns_inputs_neurons["ec3_360_to_cck_ca1"] = [25]

# to CA3
syn_inputs_neurons["sep_180_to_bas_ca3"]    = [["somaGABA_sep180"]]                          # [["somaGABAf"]]
syn_inputs_neurons["sep_360_to_olm_ca3"]    = [["somaGABA_sep360"]]                          # [["somaGABAf"]]
syn_inputs_neurons["ec2_180_to_pyr_ca3"]    = [["Adend3AMPA_ec2180", "Adend3NMDA_ec2180"]]   # [["Adend3AMPAf","Adend3NMDA"]]
syn_inputs_neurons["ec2_180_to_bas_ca3"]    = [["somaAMPA_ec2180",   "somaNMDA_ec2180"]]     # [["somaAMPAf", "somaNMDA"]]
syn_inputs_neurons["ec2_360_to_pyr_ca3"]    = [["Adend3AMPA_ec2360", "Adend3NMDA_ec2360"]]   # [["Adend3AMPAf","Adend3NMDA"]]
syn_inputs_neurons["ec2_360_to_bas_ca3"]    = [["somaAMPA_ec2360",   "somaNMDA_ec2360"]]     # [["somaAMPAf", "somaNMDA"]]
syn_inputs_neurons["dg_regular_to_pyr_ca3"] = [["Adend1AMPA_dgreg","Adend1NMDA_dgreg"]]      # [["Adend1AMPAf","Adend1NMDA"]]
syn_inputs_neurons["dg_regular_to_bas_ca3"] = [["somaAMPA_dgreg","somaNMDA_dgreg"]]          # [["somaAMPAf", "somaNMDA"]]
syn_inputs_neurons["dg_burst_to_pyr_ca3"]   = [["Adend1AMPA_dgburst","Adend1NMDA_dgburst"]]  # [["Adend1AMPAf","Adend1NMDA"]]
syn_inputs_neurons["dg_burst_to_bas_ca3"]   = [["somaAMPA_dgburst", "somaNMDA_dgburst"]]     # [["somaAMPAf", "somaNMDA"]]
# to CA1
#syn_inputs_neurons["sep_180_to_bas_ca1"] = [["somaGABA_sep180"]]                             # [["somaGABAf"]]
#syn_inputs_neurons["sep_360_to_olm_ca1"] = [["somaGABA_sep360"]]                             # [["somaGABAf"]]
#syn_inputs_neurons["ec3_180_to_pyr_ca1"] = [["Adend3AMPA_ec3180","Adend3NMDA_ec3180"]]       # [["Adend3AMPAf", "Adend3NMDA"]]
#syn_inputs_neurons["ec3_180_to_bas_ca1"] = [["somaAMPA_ec3180","somaNMDA_ec3180"]]           # [["somaAMPAf", "somaNMDA"]]
#syn_inputs_neurons["ec3_180_to_cck_ca1"] = [["somaAMPA_ec3180","somaNMDA_ec3180"]]           # [["somaAMPAf", "somaNMDA"]]
#syn_inputs_neurons["ec3_360_to_pyr_ca1"] = [["Adend3AMPA_ec3360","Adend3NMDA_ec3360"]]       # [["Adend3AMPAf", "Adend3NMDA"]]
#syn_inputs_neurons["ec3_360_to_bas_ca1"] = [["somaAMPA_ec3360","somaNMDA_ec3360"]]           # [["somaAMPAf", "somaNMDA"]]
#syn_inputs_neurons["ec3_360_to_cck_ca1"] = [["somaAMPA_ec3360","somaNMDA_ec3360"]]           # [["somaAMPAf", "somaNMDA"]]

# to CA3 # added in 15/11
delay_ext = 10.0
delay_inputs_neurons["sep_180_to_bas_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["sep_360_to_olm_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["ec2_180_to_pyr_ca3"] = [delay_ext,1.0] # AMPA/ NMDA
delay_inputs_neurons["ec2_180_to_bas_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["ec2_360_to_pyr_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["ec2_360_to_bas_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["dg_regular_to_pyr_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["dg_regular_to_bas_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["dg_burst_to_pyr_ca3"] = [delay_ext,1.0]
delay_inputs_neurons["dg_burst_to_bas_ca3"] = [delay_ext,1.0]
# to CA1
#delay_inputs_neurons["sep_180_to_bas_ca1"] = [2.0,1.0]
#delay_inputs_neurons["sep_360_to_olm_ca1"] = [2.0,1.0]
#delay_inputs_neurons["sep_180_to_cck_ca1"] = [2.0,1.0]
#delay_inputs_neurons["ec3_180_to_pyr_ca1"] = [2.0,1.0]
#delay_inputs_neurons["ec3_180_to_bas_ca1"] = [2.0,1.0]
#delay_inputs_neurons["ec3_180_to_cck_ca1"] = [2.0,1.0]
#delay_inputs_neurons["ec3_360_to_pyr_ca1"] = [2.0,1.0]
#delay_inputs_neurons["ec3_360_to_bas_ca1"] = [2.0,1.0]
#delay_inputs_neurons["ec3_360_to_cck_ca1"] = [2.0,1.0]


'''#############################################################################
                    Parameters: Background
#############################################################################'''

keys = ["pyr_ca3","bas_ca3","olm_ca3","pyr_ca1", "bas_ca1","olm_ca1","cck_ca1"]
weights_noise_neurons = dict.fromkeys(keys)
for key in keys:
    weights_noise_neurons[key] = {}

# weights_noise_neurons["pyr_ca3"]["somaAMPAf"]   = 0.05e-3
# weights_noise_neurons["pyr_ca3"]["somaGABAf"]   = 0.012e-3
weights_noise_neurons["pyr_ca3"]["somaAMPA_noise"]   = 0.05e-3
weights_noise_neurons["pyr_ca3"]["somaGABA_noise"]   = 0.012e-3
                                             
# weights_noise_neurons["pyr_ca3"]["Adend3AMPAf"] = 2*0.05e-3
# weights_noise_neurons["pyr_ca3"]["Adend3GABAf"] = 0.012e-3
# weights_noise_neurons["pyr_ca3"]["Adend3NMDA"]  = 0.0
weights_noise_neurons["pyr_ca3"]["Adend3AMPA_noise"] = 2*0.05e-3
weights_noise_neurons["pyr_ca3"]["Adend3GABA_noise"] = 0.012e-3
# weights_noise_neurons["pyr_ca3"]["Adend3NMDA_noise"] = 0.0
                                             
# weights_noise_neurons["bas_ca3"]["somaAMPAf"]   = 0.02e-3
# weights_noise_neurons["bas_ca3"]["somaGABAf"]   = 0.2e-3
weights_noise_neurons["bas_ca3"]["somaAMPA_noise"]   = 0.02e-3
weights_noise_neurons["bas_ca3"]["somaGABA_noise"]   = 0.2e-3
                                             
# weights_noise_neurons["olm_ca3"]["somaAMPAf"]   = 0.0625e-3
# weights_noise_neurons["olm_ca3"]["somaGABAf"]   = 0.2e-3
weights_noise_neurons["olm_ca3"]["somaAMPA_noise"]   = 0.0625e-3
weights_noise_neurons["olm_ca3"]["somaGABA_noise"]   = 0.2e-3
                                             
# weights_noise_neurons["pyr_ca1"]["somaAMPAf"]   = 0.05e-3
# weights_noise_neurons["pyr_ca1"]["somaGABAf"]   = 0.012e-3
weights_noise_neurons["pyr_ca1"]["somaAMPA_noise"]   = 0.05e-3
weights_noise_neurons["pyr_ca1"]["somaGABA_noise"]   = 0.012e-3

# weights_noise_neurons["pyr_ca1"]["Adend3AMPAf"] = 2*0.05e-3
# weights_noise_neurons["pyr_ca1"]["Adend3GABAf"] = 0.012e-3
# weights_noise_neurons["pyr_ca1"]["Adend3NMDA"]  = 0.0
weights_noise_neurons["pyr_ca1"]["Adend3AMPA_noise"] = 2*0.05e-3
weights_noise_neurons["pyr_ca1"]["Adend3GABA_noise"] = 0.012e-3
# weights_noise_neurons["pyr_ca1"]["Adend3NMDA_noise"]  = 0.0

# weights_noise_neurons["bas_ca1"]["somaAMPAf"]   = 0.02e-3
# weights_noise_neurons["bas_ca1"]["somaGABAf"]   = 0.2e-3
weights_noise_neurons["bas_ca1"]["somaAMPA_noise"]   = 0.02e-3
weights_noise_neurons["bas_ca1"]["somaGABA_noise"]   = 0.2e-3
                                             
# weights_noise_neurons["olm_ca1"]["somaAMPAf"]   = 0.0625e-3
# weights_noise_neurons["olm_ca1"]["somaGABAf"]   = 0.2e-3
weights_noise_neurons["olm_ca1"]["somaAMPA_noise"]   = 0.0625e-3
weights_noise_neurons["olm_ca1"]["somaGABA_noise"]   = 0.2e-3
                                             
# weights_noise_neurons["cck_ca1"]["somaAMPAf"]   = 0.045e-3 # updated on 15/11/2022
# weights_noise_neurons["cck_ca1"]["somaGABAf"]   = 0.012e-3
weights_noise_neurons["cck_ca1"]["somaAMPA_noise"]   = 0.045e-3 # updated on 15/11/2022
weights_noise_neurons["cck_ca1"]["somaGABA_noise"]   = 0.012e-3
