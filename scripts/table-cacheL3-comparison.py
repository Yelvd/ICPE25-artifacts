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

arm_metrics = ["instructions", "l3d_cache_refill ", "l3d_cache "]
intel_metrics = ["instructions", "LLC-loads ", "LLC-stores ", "LLC-load-misses ", "LLC-store-misses "]
metrics = {"Intel SPR DDR": intel_metrics, "Intel SPR HBM": intel_metrics, "Ampere Altra": arm_metrics}

df = profile.load_profile(results_dir, metrics, exp, "perf")
print(df)
df = df.loc[(df['benchmark'].isin(["S1", "S2", "W1"]))]

df = df.sort_values(by=['platform', 'tasks', 'benchmark'])
df['metric'] = [translate[m.strip()] for m in df['metric']]


new_df = []

df['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}" for _, r in df.iterrows()]
print(df)

new_columns = ["benchmark", "platform", "tasks", names["LLC_accesses"], names["LLC_misses"], "miss-rate"]
for job in np.unique(df['unique']):
    tmp_df = df.loc[df['unique'] == job]

    i = tmp_df.loc[tmp_df['metric'] == names['instructions']]['sum'].values[0] / scaling

    ca = cr = si = miss_rate = 0
    if names["LLC_accesses"] in tmp_df['metric'].values:
        ca = tmp_df.loc[tmp_df['metric'] == names['LLC_accesses']]['sum'].values[0] / scaling
        cr = tmp_df.loc[tmp_df['metric'] == names['LLC_misses']]['sum'].values[0] / scaling
    else:
        ca = tmp_df.loc[tmp_df['metric'] == names['LLC_load_accesses']]['sum'].values[0] / scaling
        ca = ca + tmp_df.loc[tmp_df['metric'] == names['LLC_store_accesses']]['sum'].values[0] / scaling
        cr = tmp_df.loc[tmp_df['metric'] == names['LLC_load_misses']]['sum'].values[0] / scaling
        cr = cr + tmp_df.loc[tmp_df['metric'] == names['LLC_store_misses']]['sum'].values[0] / scaling
    print(ca)

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
    if names["LLC_misses"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if names["LLC_accesses"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if "miss-rate" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))

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
# new_df = new_df.drop(names['instructions'], axis=1)


if exp:
    column_format = "lrrrrr"
else:
    column_format = "lrlrrrr"

print(new_df)
latex_str = new_df.style.hide().to_latex(hrules=True, column_format='lrrrr')
print(latex_str)
with open(snakemake.output["tex"], "w") as file1:
    file1.write(latex_str)
