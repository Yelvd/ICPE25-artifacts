import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tools import profile
from os import listdir
from os.path import isfile, join
from os import walk
from tools import misc


names = misc.names_for_hw_events
translate = misc.translate_hw_events
scaling = 10**12

arm_metrics = ["cpu-cycles","instructions", "l1d_cache ", "l1d_cache_refill "]
archer_metrics         = ["cpu-cycles", "ex_ret_instr", "all_dc_accesses", "l2_cache_accesses_from_dc_misses"]
snellius_genoa_metrics = ["cpu-cycles", "ex_ret_instr", 'all_data_cache_accesses', 'l1_data_cache_fills_all']
AMD_bergamo = ["cpu-cycles", "all_data_cache_accesses", "ex_ret_instr", 'l2_request_g1']
intel_metrics = ["cpu-cycles", "L1-dcache-load-misses", "L1-dcache-loads", "instructions"]
metrics = {"AMD Rome 7H12": archer_metrics, "AMD Genoa": snellius_genoa_metrics, "AMD Rome 7742": archer_metrics, "Ampere Altra": arm_metrics, "AMD Bergamo": AMD_bergamo, "Intel SPR DDR": intel_metrics, "Intel SPR HBM": intel_metrics}

results_dir = snakemake.params["resultsDir"]
exp = snakemake.params["experiment"] 
exp = exp if not exp == "all" else None

df = profile.load_profile(results_dir, metrics, exp, "perf")

df = df.sort_values(by=['platform', 'tasks', 'benchmark'])
df['metric'] = [translate[m.strip()] for m in df['metric']]

new_columns = ["Benchmark", "Platform", "Processes", "time", names["instructions"], names["cycles"], "CPI", "mem-request/inst"]
new_df = []

df['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}" for _, r in df.iterrows()]

for job in np.unique(df['unique']):
    tmp_df = df.loc[df['unique'] == job]

    r = tmp_df.loc[tmp_df['metric'] == names['instructions']]

    platform = r['platform'].values[0]
    benchmark = r['benchmark'].values[0]
    if benchmark not in ['S1', "S2", "W1"]:
        continue
    tasks = r['tasks'].values[0]


    i = tmp_df.loc[tmp_df['metric'] == names['instructions']]['sum'].values[0] / scaling
    c = tmp_df.loc[tmp_df['metric'] == names['cycles']]['sum'].values[0] / scaling
    cpi = (c / i) 
    c = c / tasks

    err = 0
    try:
        t = tmp_df.loc[tmp_df['metric'] == names["runtime"]]['avg'].values[0]
    except:
        continue


    ca = dpi = 0
    if names["l1_accesses"] in tmp_df['metric'].values:
        ca = tmp_df.loc[tmp_df['metric'] == names['l1_accesses']]['sum'].values[0] / scaling
        dpi = ca / i
    else:
        ca = tmp_df.loc[tmp_df['metric'] == names['l1_load_accesses']]['sum'].values[0] / scaling
        dpi = ca / i

    # if platform in ["AMD Rome 7742", "AMD Rome 7H12", "AMD Genoa"]:
        # cpi = cpi * tasks
    # else:
        # c = c / tasks

    r = tmp_df.loc[tmp_df['metric'] == names['cycles']]
    new_df.append([benchmark, platform, tasks, t, i, c, cpi, dpi])

new_df = pd.DataFrame(new_df, columns=new_columns)
new_df = new_df.groupby(['Platform', 'Processes', "Benchmark"]).mean()

if exp:
    new_df = new_df.reset_index().set_index(['Platform', 'Processes']).drop(columns=['Benchmark'])

for c in np.array(new_df.columns):
    if names["instructions"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if names["cycles"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if "CPI" in c:
        new_df[c] = new_df[c].map(lambda x: '{0:.2f}'.format(x))
    if "mem-request/inst" in c:
        new_df[c] = new_df[c].map(lambda x: '{0:.2f}'.format(x))


plat = ""
new_df = new_df.reset_index()
new_platform = []
new_df.loc[new_df['Platform'] == "Intel SPR HBM", 'mem-request/inst'] = [f"*{v}" for v in new_df.loc[new_df['Platform'] == "Intel SPR HBM"]['mem-request/inst'].values]
new_df.loc[new_df['Platform'] == "Intel SPR DDR", 'mem-request/inst'] = [f"*{v}" for v in new_df.loc[new_df['Platform'] == "Intel SPR DDR"]['mem-request/inst'].values]
for p in new_df['Platform']:
    if plat == p:
        new_platform.append(" ")
    else:
        plat = p
        new_platform.append(p)

new_df['Platform'] = new_platform

latex_str = new_df.style.format().hide(axis="index").to_latex(hrules=True, column_format='lrrrrr')
with open(snakemake.output["tex"], "w") as file1:
    file1.write(latex_str)
