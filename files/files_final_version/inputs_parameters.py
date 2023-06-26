# dictionary for the generation of the STATIC external inputs 

# constants
time_initial = 50
theta_rythm = 125 

inputs = ["baseline","ec2_faster_1","ec2_faster_2","ec3_faster_1","ec3_faster_2","ec_slower"]
labels = ['sep_180', 'sep_360', 'ec2_180', 'ec2_360', 'ec3_180','ec3_360','dg_regular', 'dg_burst']
params = ["time_start", "ncells", "sigma","period"]
          
phase_shift = [0.1*i*theta_rythm for i in range(1,10)] # provisional values  

for i in range(len(phase_shift)):
    inputs.append("ec2_shifted_"+str(i+1))
for i in range(len(phase_shift)):
    inputs.append("dg_shifted_"+str(i+1))

inputs_params = dict.fromkeys(inputs)

for key1 in inputs_params.keys():
    inputs_params[key1] = dict.fromkeys(labels)
    for key2 in inputs_params[key1].keys():
        inputs_params[key1][key2] = dict.fromkeys(params)

# in all the inputs options, the number of cells and sigma is the same for all the labels
for key1 in inputs_params.keys():
    inputs_params[key1]["sep_180"]["ncells"]    = 10
    inputs_params[key1]["sep_360"]["ncells"]    = 10
    inputs_params[key1]["ec2_180"]["ncells"]    = 1000
    inputs_params[key1]["ec2_360"]["ncells"]    = 1000
    inputs_params[key1]["ec3_180"]["ncells"]    = 1000
    inputs_params[key1]["ec3_360"]["ncells"]    = 1000
    inputs_params[key1]["dg_regular"]["ncells"] = 1000
    inputs_params[key1]["dg_burst"]["ncells"]   = 1000

    inputs_params[key1]["sep_180"]["sigma"]    = 0.2*theta_rythm
    inputs_params[key1]["sep_360"]["sigma"]    = 0.2*theta_rythm
    inputs_params[key1]["ec2_180"]["sigma"]    = 0.2*theta_rythm
    inputs_params[key1]["ec2_360"]["sigma"]    = 0.2*theta_rythm
    inputs_params[key1]["ec3_180"]["sigma"]    = 0.2*theta_rythm
    inputs_params[key1]["ec3_360"]["sigma"]    = 0.2*theta_rythm
    inputs_params[key1]["dg_regular"]["sigma"] = 0.2*theta_rythm
    inputs_params[key1]["dg_burst"]["sigma"]   = 0.2*theta_rythm

# baseline parameters
inputs_params["baseline"]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["baseline"]["sep_360"]["time_start"]    = time_initial
inputs_params["baseline"]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["baseline"]["ec2_360"]["time_start"]    = time_initial
inputs_params["baseline"]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["baseline"]["ec3_360"]["time_start"]    = time_initial
inputs_params["baseline"]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm
inputs_params["baseline"]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm

inputs_params["baseline"]["sep_180"]["period"]    = theta_rythm
inputs_params["baseline"]["sep_360"]["period"]    = theta_rythm
inputs_params["baseline"]["ec2_180"]["period"]    = theta_rythm
inputs_params["baseline"]["ec2_360"]["period"]    = theta_rythm
inputs_params["baseline"]["ec3_180"]["period"]    = theta_rythm
inputs_params["baseline"]["ec3_360"]["period"]    = theta_rythm
inputs_params["baseline"]["dg_regular"]["period"] = theta_rythm
inputs_params["baseline"]["dg_burst"]["period"]   = theta_rythm
#####################################################################################

# EC2 with higer frequency +0.25 Hz
inputs_params["ec2_faster_1"]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec2_faster_1"]["sep_360"]["time_start"]    = time_initial
inputs_params["ec2_faster_1"]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec2_faster_1"]["ec2_360"]["time_start"]    = time_initial
inputs_params["ec2_faster_1"]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec2_faster_1"]["ec3_360"]["time_start"]    = time_initial
inputs_params["ec2_faster_1"]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm
inputs_params["ec2_faster_1"]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm

inputs_params["ec2_faster_1"]["sep_180"]["period"]    = theta_rythm
inputs_params["ec2_faster_1"]["sep_360"]["period"]    = theta_rythm
inputs_params["ec2_faster_1"]["ec2_180"]["period"]    = 121.2 #theta_rythm
inputs_params["ec2_faster_1"]["ec2_360"]["period"]    = 121.2 #theta_rythm
inputs_params["ec2_faster_1"]["ec3_180"]["period"]    = theta_rythm
inputs_params["ec2_faster_1"]["ec3_360"]["period"]    = theta_rythm
inputs_params["ec2_faster_1"]["dg_regular"]["period"] = theta_rythm
inputs_params["ec2_faster_1"]["dg_burst"]["period"]   = theta_rythm
##########################################################################################

# EC2 with higer frequency +0.5 Hz
inputs_params["ec2_faster_2"]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec2_faster_2"]["sep_360"]["time_start"]    = time_initial
inputs_params["ec2_faster_2"]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec2_faster_2"]["ec2_360"]["time_start"]    = time_initial
inputs_params["ec2_faster_2"]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec2_faster_2"]["ec3_360"]["time_start"]    = time_initial
inputs_params["ec2_faster_2"]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm
inputs_params["ec2_faster_2"]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm

inputs_params["ec2_faster_2"]["sep_180"]["period"]    = theta_rythm
inputs_params["ec2_faster_2"]["sep_360"]["period"]    = theta_rythm
inputs_params["ec2_faster_2"]["ec2_180"]["period"]    = 117.6 #theta_rythm
inputs_params["ec2_faster_2"]["ec2_360"]["period"]    = 117.6 #theta_rythm
inputs_params["ec2_faster_2"]["ec3_180"]["period"]    = theta_rythm
inputs_params["ec2_faster_2"]["ec3_360"]["period"]    = theta_rythm
inputs_params["ec2_faster_2"]["dg_regular"]["period"] = theta_rythm
inputs_params["ec2_faster_2"]["dg_burst"]["period"]   = theta_rythm
##########################################################################################

# EC3 with higer frequency +0.25 Hz
inputs_params["ec3_faster_1"]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec3_faster_1"]["sep_360"]["time_start"]    = time_initial
inputs_params["ec3_faster_1"]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec3_faster_1"]["ec2_360"]["time_start"]    = time_initial
inputs_params["ec3_faster_1"]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec3_faster_1"]["ec3_360"]["time_start"]    = time_initial
inputs_params["ec3_faster_1"]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm
inputs_params["ec3_faster_1"]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm

inputs_params["ec3_faster_1"]["sep_180"]["period"]    = theta_rythm
inputs_params["ec3_faster_1"]["sep_360"]["period"]    = theta_rythm
inputs_params["ec3_faster_1"]["ec2_180"]["period"]    = theta_rythm
inputs_params["ec3_faster_1"]["ec2_360"]["period"]    = theta_rythm
inputs_params["ec3_faster_1"]["ec3_180"]["period"]    = 121.2 #theta_rythm
inputs_params["ec3_faster_1"]["ec3_360"]["period"]    = 121.2 #theta_rythm
inputs_params["ec3_faster_1"]["dg_regular"]["period"] = theta_rythm
inputs_params["ec3_faster_1"]["dg_burst"]["period"]   = theta_rythm
##########################################################################################

# EC3 with higer frequency +0.5 Hz
inputs_params["ec3_faster_2"]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec3_faster_2"]["sep_360"]["time_start"]    = time_initial
inputs_params["ec3_faster_2"]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec3_faster_2"]["ec2_360"]["time_start"]    = time_initial
inputs_params["ec3_faster_2"]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec3_faster_2"]["ec3_360"]["time_start"]    = time_initial
inputs_params["ec3_faster_2"]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm
inputs_params["ec3_faster_2"]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm

inputs_params["ec3_faster_2"]["sep_180"]["period"]    = theta_rythm
inputs_params["ec3_faster_2"]["sep_360"]["period"]    = theta_rythm
inputs_params["ec3_faster_2"]["ec2_180"]["period"]    = theta_rythm
inputs_params["ec3_faster_2"]["ec2_360"]["period"]    = theta_rythm
inputs_params["ec3_faster_2"]["ec3_180"]["period"]    = 117.6 #theta_rythm
inputs_params["ec3_faster_2"]["ec3_360"]["period"]    = 117.6 #theta_rythm
inputs_params["ec3_faster_2"]["dg_regular"]["period"] = theta_rythm
inputs_params["ec3_faster_2"]["dg_burst"]["period"]   = theta_rythm

##########################################################################################

# EC with lower  frequency -0.25 Hz
inputs_params["ec_slower"]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec_slower"]["sep_360"]["time_start"]    = time_initial
inputs_params["ec_slower"]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec_slower"]["ec2_360"]["time_start"]    = time_initial
inputs_params["ec_slower"]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
inputs_params["ec_slower"]["ec3_360"]["time_start"]    = time_initial
inputs_params["ec_slower"]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm
inputs_params["ec_slower"]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm

inputs_params["ec_slower"]["sep_180"]["period"]    = theta_rythm
inputs_params["ec_slower"]["sep_360"]["period"]    = theta_rythm
inputs_params["ec_slower"]["ec2_180"]["period"]    = 129.0
inputs_params["ec_slower"]["ec2_360"]["period"]    = 129.0
inputs_params["ec_slower"]["ec3_180"]["period"]    = 129.0
inputs_params["ec_slower"]["ec3_360"]["period"]    = 129.0
inputs_params["ec_slower"]["dg_regular"]["period"] = theta_rythm
inputs_params["ec_slower"]["dg_burst"]["period"]   = theta_rythm

#############################################################################################
# EC2 shifted 
inputs_ec2 = [f"ec2_shifted_{i}" for i in range(1,len(phase_shift)+1)]
for i,keys in enumerate(inputs_ec2):
    inputs_params[keys]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
    inputs_params[keys]["sep_360"]["time_start"]    = time_initial
    inputs_params[keys]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0+phase_shift[i]
    inputs_params[keys]["ec2_360"]["time_start"]    = time_initial+phase_shift[i]
    inputs_params[keys]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
    inputs_params[keys]["ec3_360"]["time_start"]    = time_initial
    inputs_params[keys]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm
    inputs_params[keys]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm

    inputs_params[keys]["sep_180"]["period"]    = theta_rythm
    inputs_params[keys]["sep_360"]["period"]    = theta_rythm
    inputs_params[keys]["ec2_180"]["period"]    = theta_rythm
    inputs_params[keys]["ec2_360"]["period"]    = theta_rythm
    inputs_params[keys]["ec3_180"]["period"]    = theta_rythm
    inputs_params[keys]["ec3_360"]["period"]    = theta_rythm
    inputs_params[keys]["dg_regular"]["period"] = theta_rythm
    inputs_params[keys]["dg_burst"]["period"]   = theta_rythm
##########################################################################################

# DG shifted
inputs_dg = [f"dg_shifted_{i}" for i in range(1,len(phase_shift)+1)]
for i,keys in enumerate(inputs_dg):
    inputs_params[keys]["sep_180"]["time_start"]    = time_initial+theta_rythm/2.0
    inputs_params[keys]["sep_360"]["time_start"]    = time_initial
    inputs_params[keys]["ec2_180"]["time_start"]    = time_initial+theta_rythm/2.0
    inputs_params[keys]["ec2_360"]["time_start"]    = time_initial
    inputs_params[keys]["ec3_180"]["time_start"]    = time_initial+theta_rythm/2.0
    inputs_params[keys]["ec3_360"]["time_start"]    = time_initial
    inputs_params[keys]["dg_regular"]["time_start"] = time_initial+0.375*theta_rythm+phase_shift[i]
    inputs_params[keys]["dg_burst"]["time_start"]   = time_initial+0.125*theta_rythm+phase_shift[i]

    inputs_params[keys]["sep_180"]["period"]    = theta_rythm
    inputs_params[keys]["sep_360"]["period"]    = theta_rythm
    inputs_params[keys]["ec2_180"]["period"]    = theta_rythm
    inputs_params[keys]["ec2_360"]["period"]    = theta_rythm
    inputs_params[keys]["ec3_180"]["period"]    = theta_rythm
    inputs_params[keys]["ec3_360"]["period"]    = theta_rythm
    inputs_params[keys]["dg_regular"]["period"] = theta_rythm
    inputs_params[keys]["dg_burst"]["period"]   = theta_rythm
##########################################################################################



# future: DG regular and burst

# and more featured inputs