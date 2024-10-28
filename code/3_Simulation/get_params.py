import numpy as np
import pandas as pd
import itertools as it

num_decimal = 2
subj = [50]
prop1 = np.arange(0.1, 1, 0.2).round(decimals = num_decimal).tolist()
prop2 = np.arange(0.1, 1, 0.2).round(decimals = num_decimal).tolist()
effsize = [0, 0.01, 0.02, 0.03]

param_arr = pd.DataFrame(list(it.product(subj, prop1, prop2, effsize)))
param_arr["dirname"] = "sim_" + param_arr.astype(str).agg("_".join, axis = 1)
param_arr.to_csv("code/3_Simulation/sim_param.txt", sep = " ", header = False, index = False)
