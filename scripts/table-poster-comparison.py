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
from tools import internal

results_dir = snakemake.params["resultsDir"]
exp = snakemake.params["experiment"] 
platform = snakemake.params["platform"]
exp = exp if not exp == "all" else None
scaling = 10**12

names = misc.names_for_hw_events
translate = misc.translate_hw_events

arm_metrics = ["INST_RETIRED", "CPU_CYCLES", "Memory data volume", "LD_SPEC", "ST_SPEC" ]
snellius_rome_metrics  = ["RETIRED_INSTRUCTIONS", "DATA_CACHE_ACCESSES", "DATA_CACHE_REFILLS_ALL"]
archer_metrics         = ["ex_ret_instr", "all_dc_accesses", "l2_cache_accesses_from_dc_misses"]
snellius_genoa_metrics = ['cpu-cycles', 'ex_ret_instr', 'all_data_cache_accesses']
intel_DDR_metrics = ["Runtime unhalted", "DDR data volume" , "INSTR_RETIRED_ANY", "CPU_CLK_UNHALTED_CORE", "MEM_INST_RETIRED_ALL_LOADS", "MEM_INST_RETIRED_ALL_STORES"]
intel_HBM_metrics = ["Runtime unhalted", "HBM data volume", "INSTR_RETIRED_ANY", "CPU_CLK_UNHALTED_CORE", "MEM_INST_RETIRED_ALL_LOADS", "MEM_INST_RETIRED_ALL_STORES"]
# metrics = {"Ampere_Q8030": arm_metrics,"Intel_SPR_HBM": intel_metrics, "Intel_SPR_DDR": intel_metrics, "snellius-rome": snellius_rome_metrics, "snellius-genoa": snellius_genoa_metrics, "archer_rome": archer_metrics}

if platform == "Intel SPR HBM":
    metrics = {"Intel SPR HBM": intel_HBM_metrics}
if platform == "Intel SPR DDR":
    metrics = {"Intel SPR DDR": intel_DDR_metrics}
if platform == "Ampere Altra":
    metrics = {"Ampere Altra": arm_metrics}
if platform == "AMD Genoa":
    metrics = {"AMD Genoa": snellius_genoa_metrics}


df = profile.load_profile(results_dir, metrics, exp, "perf")

df = df.sort_values(by=['platform', 'tasks', 'benchmark'])
df['metric'] = [translate[m] for m in df['metric']]

new_df = []

df_time = internal.load_time(results_dir, exp)
df_time = df_time.loc[df_time['platform'] == platform]
df_time = df_time.groupby(["platform", "tasks"], as_index=False).agg({"hemocell": ['mean', 'std']})
df_time = df_time.sort_values(by=["tasks"])
df_time = df_time.round(0)


df['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}" for _, r in df.iterrows()]

new_columns = ["benchmark", 
        "platform", 
        "tasks", 
        names["instructions"], 
        names["cycles"], 
        "CPI",
        names["runtime"], 
        names["data volume"], 
        names["mem bandwidth"]]


print(df)
for job in np.unique(df['unique']):
    tmp_df = df.loc[df['unique'] == job]

    r = tmp_df.loc[tmp_df['metric'] == names['instructions']]
    tasks = r['tasks'].values[0]

    v = tmp_df[tmp_df['metric'] == names['data volume']]['max'].values[0]*2
    # if platform == "AMD Genoa":
        # access = tmp_df[tmp_df['metric'] == names["l1_cache_accesses"]]['sum'].values[0]
    # else:
        # loads = tmp_df[tmp_df['metric'] == names['loads']]['sum'].values[0]
        # stores = tmp_df[tmp_df['metric'] == names['stores']]['sum'].values[0]
        # access = loads + stores

    # v = 1.0E-09 * (access) * 64 / 8
    v = (v)
    run = df_time.loc[df_time["tasks"] == tasks]['hemocell']['mean'].values[0]

    b = v / run
    tb = b

    i = tmp_df.loc[tmp_df['metric'] == names['instructions']]['sum'].values[0] / scaling
    c = tmp_df.loc[tmp_df['metric'] == names['cycles']]['max'].values[0] / scaling
    # cpi = (i / c) / tmp_df.loc[tmp_df['metric'] == 'cycles']['tasks'].values[0]
    cpi = (c / i) * tmp_df.loc[tmp_df['metric'] == names['cycles']]['tasks'].values[0]

    new_df.append([r['benchmark'].values[0], r['platform'].values[0], r['tasks'].values[0], i, c, cpi, run, (v)/1000, tb])

new_df = pd.DataFrame(new_df, columns=new_columns)
new_df['platform'] = [p for p in new_df['platform']]
new_df = new_df.groupby(['platform', 'tasks', "benchmark"]).mean()


if exp:
    new_df = new_df.reset_index().set_index(['platform', 'tasks']).drop(columns=['benchmark'])



for c in np.array(new_df.columns):
    if names["instructions"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if names["cycles"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if "CPI" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if names["runtime"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if "avg data volume" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if "mem bandwidth per task" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if names["mem bandwidth"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if names["data volume"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))

plat = ""
new_df = new_df.reset_index()
new_platform = []
for p in new_df['platform']:
    if plat == p:
        new_platform.append(" ")
    else:
        plat = p
        new_platform.append(p)

new_df = new_df.rename(columns={"tasks": "Processes"})
print(new_df)
new_df = new_df.drop(['platform'], axis=1)
new_df = new_df.set_index('Processes').T
new_df = new_df.reset_index()
print(new_df)

latex_str = new_df.style.hide().to_latex(hrules=True, column_format='lrr')
print(latex_str)
with open(snakemake.output["tex"], "w") as file1:
    file1.write(latex_str)
