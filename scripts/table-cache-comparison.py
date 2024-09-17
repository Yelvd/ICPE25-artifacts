import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tools import profile
from tools import misc
from os import walk
from os import listdir
from os.path import isfile, join

results_dir = snakemake.params["resultsDir"]
exp = snakemake.params["experiment"] 
exp = exp if not exp == "all" else None
scaling = 10**12

names = misc.names_for_hw_events
translate = misc.translate_hw_events

arm_metrics = ["instructions", "l1d_cache ", "l1d_cache_refill "]
snellius_rome_metrics  = ["RETIRED_INSTRUCTIONS", "DATA_CACHE_ACCESSES", "DATA_CACHE_REFILLS_ALL"]
archer_metrics         = ["ex_ret_instr", "all_dc_accesses", "l2_cache_accesses_from_dc_misses"]
snellius_genoa_metrics = ["ex_ret_instr", 'all_data_cache_accesses', 'l1_data_cache_fills_all']
AMD_bergamo = ["all_data_cache_accesses", "ex_ret_instr", 'l2_request_g1']
intel_metrics = ["L1-dcache-load-misses", "L1-dcache-loads", "instructions"]
metrics = {"AMD Rome 7H12": archer_metrics, "AMD Genoa": snellius_genoa_metrics, "AMD Rome 7742": archer_metrics, "Ampere Altra": arm_metrics, "AMD Bergamo": AMD_bergamo, "Intel SPR DDR": intel_metrics, "Intel SPR HBM": intel_metrics}

df = profile.load_profile(results_dir, metrics, exp, "perf")
df = df.loc[(df['benchmark'].isin(["S1", "S2", "W1"]))]

df = df.sort_values(by=['platform', 'tasks', 'benchmark'])
df['metric'] = [translate[m.strip()] for m in df['metric']]


new_df = []

df['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}" for _, r in df.iterrows()]

new_columns = ["benchmark", "platform", "tasks", names["l1_accesses"], names["l1_misses"], "miss-rate"]
for job in np.unique(df['unique']):
    tmp_df = df.loc[df['unique'] == job]

    i = tmp_df.loc[tmp_df['metric'] == names['instructions']]['sum'].values[0] / scaling

    ca = cr = si = miss_rate = 0
    if names["l1_accesses"] in tmp_df['metric'].values:
        ca = tmp_df.loc[tmp_df['metric'] == names['l1_accesses']]['sum'].values[0] / scaling
        cr = tmp_df.loc[tmp_df['metric'] == names['l1_misses']]['sum'].values[0] / scaling
    else:
        ca = tmp_df.loc[tmp_df['metric'] == names['l1_load_accesses']]['sum'].values[0] / scaling
        cr = tmp_df.loc[tmp_df['metric'] == names['l1_load_misses']]['sum'].values[0] / scaling

    si = (cr * scaling)*64 * 10**-9
    miss_rate = cr / ca 

    r = tmp_df.loc[tmp_df['metric'] == names['instructions']]
    new_df.append([r['benchmark'].values[0], r['platform'].values[0], r['tasks'].values[0], ca, cr, miss_rate])

new_df = pd.DataFrame(new_df, columns=new_columns)
new_df['platform'] = [p for p in new_df['platform']]
new_df = new_df.groupby(['platform', 'tasks', "benchmark"]).mean()

if exp:
    new_df = new_df.reset_index().set_index(['platform', 'tasks']).drop(columns=['benchmark'])



for c in np.array(new_df.columns):
    if names["instructions"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if "cycles" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2E}'.format(x))
    if names["l1_misses"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if names["l1_accesses"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if "miss-rate" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.3f}'.format(x))

plat = ""
new_df = new_df.reset_index()
new_platform = []
for p in new_df['platform']:
    if plat == p:
        new_platform.append(" ")
    else:
        plat = p
        new_platform.append(p if not "Intel" in p else f"{p}*")

new_df['platform'] = new_platform
# new_df = new_df.drop(names['instructions'], axis=1)


if exp:
    column_format = "lrrrrr"
else:
    column_format = "lrlrrrr"

latex_str = new_df.style.hide().to_latex(hrules=True, column_format='lrrrr')
with open(snakemake.output["tex"], "w") as file1:
    file1.write(latex_str)
