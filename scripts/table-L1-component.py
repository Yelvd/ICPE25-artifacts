import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tools import profile
from tools import scorep
from os import listdir
from os.path import isfile, join
from os import walk
from tools import misc


names = misc.names_for_hw_events
translate = misc.translate_hw_events
scaling = 10**12

results_dir = snakemake.params["resultsDir"]
exp = snakemake.params["experiment"] 
exp = exp if not exp == "all" else None

df = scorep.load_scorep(results_dir, exp)
df = df.loc[df["region"].isin(["fluid", "particle", "couple"])]

df = df.sort_values(by=['platform', 'tasks', 'benchmark'])
f = scorep.load_scorep("./results-components/")
df.loc[df["region"].isin(["fluid", "particle", "couple"])]
df['metric'] = [translate[m.strip()] for m in df['metric']]

new_columns = ["Benchmark", "Platform", "Processes", "Region", names["l1_accesses"], names["l1_misses"], "miss-rate"]
new_df = []

df['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}-{r['region']}" for _, r in df.iterrows()]

for job in np.unique(df['unique']):
    tmp_df = df.loc[df['unique'] == job]

    r = tmp_df.loc[tmp_df['metric'] == "Instructions"]

    platform = r['platform'].values[0]
    benchmark = r['benchmark'].values[0]
    tasks = int(r['tasks'].values[0])
    region = r['region'].values[0]

    i = tmp_df.loc[tmp_df['metric'] == names['instructions']]['sum'].values[0] / scaling
    c = tmp_df.loc[tmp_df['metric'] == names['cycles']]['sum'].values[0] / scaling
    t = tmp_df.loc[tmp_df['metric'] == names['total']]['mean'].values[0]

    cpi = (c / i) 
    c = c / tasks

    ca = cr = si = miss_rate = 0
    ca = tmp_df.loc[tmp_df['metric'] == names['l1_accesses']]['sum'].values[0] / scaling
    cr = tmp_df.loc[tmp_df['metric'] == names['l1_misses']]['sum'].values[0] / scaling

    si = (cr * scaling)*64 * 10**-9
    miss_rate = cr / ca 

    new_df.append([benchmark, platform, tasks, region, ca, cr, miss_rate])

new_df = pd.DataFrame(new_df, columns=new_columns)
new_df = new_df.groupby(['Platform', 'Processes', "Benchmark", "Region"]).mean()

if exp:
    new_df = new_df.reset_index().set_index(['Platform', 'Processes']).drop(columns=['Benchmark'])

for c in np.array(new_df.columns):
    if names["instructions"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.0f}'.format(x))
    if names["cycles"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if names["l1_misses"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if names["l1_accesses"] in c:
        new_df[c] = new_df[c].map(lambda x: '{:.2f}'.format(x))
    if "CPI" in c:
        new_df[c] = new_df[c].map(lambda x: '{0:.2f}'.format(x))
    if "mem-request/inst" in c:
        new_df[c] = new_df[c].map(lambda x: '{0:.2f}'.format(x))
    if names["runtime"] in c:
        new_df[c] = new_df[c].map(lambda x: '{0:.0f}'.format(x))
    if "miss-rate" in c:
        new_df[c] = new_df[c].map(lambda x: '{:.3f}'.format(x))


plat = ""
new_df = new_df.reset_index()
new_platform = []
# new_df.loc[new_df['Platform'] == "Intel SPR HBM", 'mem-request/inst'] = [f"*{v}" for v in new_df.loc[new_df['Platform'] == "Intel SPR HBM"]['mem-request/inst'].values]
# new_df.loc[new_df['Platform'] == "Intel SPR DDR", 'mem-request/inst'] = [f"*{v}" for v in new_df.loc[new_df['Platform'] == "Intel SPR DDR"]['mem-request/inst'].values]
for p in new_df['Platform']:
    if plat == p:
        new_platform.append(" ")
    else:
        plat = p
        new_platform.append(p)

new_df['Platform'] = new_platform

print(new_df)
latex_str = new_df.style.format().hide(axis="index").to_latex(hrules=True, column_format='lrrrrr')
print(latex_str)
# with open(snakemake.output["tex"], "w") as file1:
    # file1.write(latex_str)
