import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tools import profile
from tools import misc
from tools import internal
from os import walk
from os import listdir
from os.path import isfile, join

results_dir = snakemake.params["resultsDir"]
exp = snakemake.params["experiment"] 
exp = exp if not exp == "all" else None
scaling = 10**12

names = misc.names_for_hw_events
translate = misc.translate_hw_events

# metrics = {"Ampere_Q8030": arm_metrics,"Intel_SPR_HBM": intel_metrics, "Intel_SPR_DDR": intel_metrics, "snellius-rome": snellius_rome_metrics, "snellius-genoa": snellius_genoa_metrics, "archer_rome": archer_metrics}
arm_metrics = ["instructions", "l3d_cache_refill"]
snellius_rome_metrics  = ["RETIRED_INSTRUCTIONS", "DATA_CACHE_ACCESSES", "DATA_CACHE_REFILLS_ALL"]
archer_metrics         = ["ex_ret_instr", "all_dc_accesses", "l2_cache_accesses_from_dc_misses"]
snellius_genoa_metrics = ["ex_ret_instr", 'all_data_cache_accesses', 'l1_data_cache_fills_all']
AMD_bergamo = ["all_data_cache_accesses", "ex_ret_instr", 'l2_request_g1']
intel_metrics = ["instructions", "LLC-load-misses", "LLC-store-misses"]
metrics = {"Intel SPR DDR": intel_metrics, "Intel SPR HBM": intel_metrics, "Ampere Altra": arm_metrics}

df = profile.load_profile(results_dir, metrics, exp, "perf")

df = df.sort_values(by=['platform', 'tasks', 'benchmark'])
df['metric'] = [translate[m] for m in df['metric']]

new_df = []

df['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}" for _, r in df.iterrows()]


df_time = internal.load_time(results_dir)
df_time = df_time.reset_index(drop=True)
df_time = df_time.groupby(["platform", "benchmark", "tasks"], as_index=False).agg({"iterate": ['mean', 'std', 'min', 'max']})


new_columns = ["benchmark", "platform", "tasks", names["runtime"], "volume [GB]", "Bandwidth [GB/s]"]
for job in np.unique(df['unique']):
    tmp_df = df.loc[df['unique'] == job]

    r = tmp_df.loc[tmp_df['metric'] == names['instructions']]

    platform = r['platform'].values[0]
    benchmark = r['benchmark'].values[0]
    tasks = r['tasks'].values[0]


    run = df_time.loc[(df_time['platform'] == platform) & (df_time['tasks'] == tasks) & (df_time["benchmark"] == benchmark)]['iterate']['mean'].values[0]
    v = 0
    if names["LLC_store_misses"] in tmp_df['metric'].values:
        v = tmp_df[tmp_df['metric'] == names['LLC_store_misses']]['avg'].values[0] + tmp_df[tmp_df['metric'] == names['LLC_load_misses']]['avg'].values[0]

    if names["LLC_misses"] in tmp_df['metric'].values:
        v = tmp_df[tmp_df['metric'] == names['LLC_misses']]['avg'].values[0]

    v = v*64 * 10**-9
    b = v / run
    tb = b * tasks

    run = tmp_df[tmp_df['metric'] == names['runtime']]['max'].values[0]

    new_df.append([benchmark, platform, tasks, run, v, b])

new_df = pd.DataFrame(new_df, columns=new_columns)
new_df['platform'] = [p for p in new_df['platform']]
new_df = new_df.groupby(['platform', 'tasks', "benchmark"]).mean()


if exp:
    new_df = new_df.reset_index().set_index(['platform', 'tasks']).drop(columns=['benchmark'])



for c in np.array(new_df.columns):
    if "Bandwidth [GB/s]" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if "volume [GB]" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if names["runtime"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if "avg data volume" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if "mem bandwidth per task" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if "mem bandwidth total" in c:
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

new_df['platform'] = new_platform

latex_str = new_df.style.hide().to_latex(hrules=True, column_format='lrrrr')
print(new_df)
print(latex_str)
# with open(snakemake.output["tex"], "w") as file1:
    # file1.write(latex_str)
