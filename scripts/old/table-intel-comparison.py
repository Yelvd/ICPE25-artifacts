import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from os import walk
from tools import misc

metrics = ['Runtime', 'CPI', 'Memory bandwidth', 'HBM bandwidth']

results_dir = snakemake.params["resultsDir"]
exp = snakemake.params["experiment"]
platform = "Intel"

def load_experiment(experiment_dir):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    files = [f for f in files if "likwid_prof" in f]
    results = {}

    with open(f"{experiment_dir}/meta.yaml") as f:
        meta = yaml.load(f, Loader=yaml.FullLoader)

    if meta["General"].get("Benchmark") is None:
        print("WARNING FIXING TO FIXED")
        results["benchmark"] = "fixed"
    else:
        results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])


    data = {}
    for m in metrics:
        data[m] = [0] * results['tasks']
    for filename in files:
        taskId = int(filename.split('.')[0].split('_')[-1])

        with open(f"{experiment_dir}/{filename}") as f:
            for l in f:
                for m in metrics:
                    if m in l:
                        data[m][taskId] = float(l.strip().split('|')[-2].strip())


    for m in metrics:
        if m == "Runtime":
            results["runtime"] = np.max(data[m])
            continue

        results[f"avg {m}"] = np.mean(data[m])
        results[f"min {m}"] = np.min(data[m])
        results[f"max {m}"] = np.max(data[m])
        results[f"sum {m}"] = np.sum(data[m])

    for k in results.keys():
        results[k] = [results[k]]

    df = pd.DataFrame.from_dict(results)
    return df

def load_data(results_dir, experiment=None):
    data = []
    w = walk(results_dir)
    for (dirpath, dirnames, filenames) in w:
        if "meta.yaml" not in filenames:
            continue

        with open(f"{dirpath}/meta.yaml") as f:
            meta = yaml.load(f, Loader=yaml.FullLoader)

        if platform not in meta["General"]["Platform"]:
            continue

        if meta["General"]["Monitor-tool"] in ["profile"]:
            data.append(load_experiment(dirpath))

    df = pd.concat(data)

    if experiment is not None:
        match experiment:
            case "fixed":
                df = df.loc[df["benchmark"].isin(["fixed", "cube-fixed"])]
            case "fixed-large":
                df = df.loc[df["benchmark"].isin(["fixed-large", "cube-fixed-large"])]
            case "relative":
                df = df.loc[df["benchmark"].isin(["fixed-large", "cube-relative"])]

    return df 


df = load_data(results_dir, exp)

df = df.sort_values(by=['benchmark', 'tasks'])

df = df.drop(['sum CPI'], axis=1)
df['benchmark'] = [misc.translate_exp[b] for b in df['benchmark']]

#DEAL with 'Memory bandwidth'
#df = df.drop(['avg L3 bandwidth', 'avg Memory bandwidth'], axis=1)

for m in ['Memory bandwidth', "HBM bandwidth"]:
    df[f"min {m}"] =  df[f"min {m}"] / 1000
    df[f"max {m}"] =  df[f"max {m}"] / 1000
    df[f"sum {m}"] =  df[f"sum {m}"] / 1000
    df[f"avg {m}"] =  df[f"avg {m}"] / 1000

df = df.rename(columns={
    'sum Memory bandwidth': 'Sum Mem BW [GB/s]', 
    'min Memory bandwidth': 'Min Mem BW [GB/s]', 
    'max Memory bandwidth': 'Max Mem BW [GB/s]', 
    'avg Memory bandwidth': 'Avg Mem BW [GB/s]', 

    'sum HBM bandwidth': 'Sum HBM BW [GB/s]', 
    'min HBM bandwidth': 'Min HBM BW [GB/s]', 
    'max HBM bandwidth': 'Max HBM BW [GB/s]', 
    'avg HBM bandwidth': 'Avg HBM BW [GB/s]', 

    'sum L3 bandwidth': 'Sum L3 BW [MB/s]', 
    'min L3 bandwidth': 'Min L3 BW [MB/s]', 
    'max L3 bandwidth': 'Max L3 BW [MB/s]', 
    'avg L3 bandwidth': 'Avg L3 BW [MB/s]', 

    'sum CPI': 'Sum CPI', 
    'min CPI': 'Min CPI', 
    'max CPI': 'Max CPI', 
    'avg CPI': 'Avg CPI',

    'runtime': 'runtime [s]',
    })

for c in np.array(df.columns):

    if "Mem" in c:
        df[c] = df[c].map(lambda x: '{0:.0f}'.format(x))

    if "L3" in c:
        df[c] = df[c].map(lambda x: '{0:.0f}'.format(x))

    if "CPI" in c:
        df[c] = df[c].map(lambda x: '{0:.2f}'.format(x))

    if "HBM" in c:
        df[c] = df[c].map(lambda x: '{0:.0f}'.format(x))

    if "runtime" in c:
        df[c] = df[c].map(lambda x: '{0:.1f}'.format(x))

print(df.columns)
hidelist = ["Min CPI", "Max CPI"]
new_df = df.drop(hidelist, axis=1)
print(new_df.style.format().hide(axis="index").to_latex(hrules=True))

print(df)
