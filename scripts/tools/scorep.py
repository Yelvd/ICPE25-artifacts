import sys
import io
import os
import pandas as pd
import numpy as np
import yaml
import sys
import re
from . import misc
from os import listdir
from os.path import isfile, join
from os import walk


def load_snellius(experiment_dir, meta):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    results = {}
    for f in files:
        if re.match("results.csv", f):
            scorep_filename = f

    filter_list = ["Cnode ID", "Thread ID", "iteration", "regionName"]
    for f in files:
        if ".out" in f and "bull" not in f and "post" not in f:
            out_filename = f

    df = pd.read_csv(f"{experiment_dir}/{scorep_filename}")
    
    # Create rows for the components
    fluid_list = ["collideandstream", "setexternalvector"]
    particle_list = ["syncenvelopes", "advanceparticles", "applyconstitutivemodel", "deletenonlocalparticles"]
    couple_list = ["spreadparticleforce", "interpolatefluidvelocity"]
    
    new_df = []
    for r, l in [("fluid", fluid_list), ("particle", particle_list), ("couple", couple_list)]:
        tmp_df = df.loc[df["regionName"].isin(l)]

        for t in np.unique(tmp_df[" Thread ID"]):
            tmp2 = tmp_df.loc[tmp_df[" Thread ID"] == t]

            new_row = []
            for c in tmp2.columns:
                if c in ["Cnode ID", " Thread ID", "iteration"]:
                    new_row.append(tmp2[c].values[0])
                elif c == "regionName":
                    new_row.append(r)
                else:
                    new_row.append(np.sum(tmp2[c].values))

            new_df.append(new_row)


    df = pd.concat([df, pd.DataFrame(new_df, columns=df.columns)])
    regions = np.unique(df["regionName"])

    columnname = ["id", "platform", "benchmark", "tasks", "region", "metric", "min", "max", "mean", "std", "sum"]
    data = []
    for r in regions:
        tmp_df = df.loc[df["regionName"] == r]

        for c in df.columns:
            if c.strip() in filter_list:
                continue

            d = tmp_df[c].values
            data.append([meta["General"]["Id"], meta["General"]["Platform"] ,meta["General"]["Benchmark"], meta["General"]["Tasks"], r, c, np.min(d), np.max(d), np.mean(d), np.std(d), np.sum(d)])


    df = pd.DataFrame(data, columns=columnname)
    return df

    
def load_experiment(experiment_dir):
    
    with open(f"{experiment_dir}/meta.yaml") as f:
        meta = yaml.load(f, Loader=yaml.FullLoader)

    if meta["General"]["Platform"] in ["snellius_rome", "snellius_genoa", "snellius-genoa", "snellius-rome"]:
        return load_snellius(experiment_dir, meta)

    # elif meta["General"]["Platform"] in ["archer_rome", "archer"]:
        # return load_archer(experiment_dir, meta)

    # else:
        # return load_eviden(experiment_dir, meta)


def load_scorep(results_dir, experiment=None):
    experiment = misc.translate_exp_to_internal[experiment]
    data = []
    w = walk(results_dir)
    for (dirpath, dirnames, filenames) in w:
        if "meta.yaml" not in filenames:
            continue

        with open(f"{dirpath}/meta.yaml") as f:
            meta = yaml.load(f, Loader=yaml.FullLoader)

        if meta["General"].get("Monitor-tool") is None:
            continue

        if meta["General"]["Monitor-tool"] == "scorep":
            data.append(load_experiment(dirpath))


    df = pd.concat(data)

    match experiment:
        case "fixed":
            df = df.loc[df["benchmark"].isin(["fixed", "cube-fixed"])]
        case "fixed-large":
            df = df.loc[df["benchmark"].isin(["fixed-large", "cube-fixed-large"])]
        case "relative":
            df = df.loc[df["benchmark"].isin(["relative", "cube-relative"])]

    df["benchmark"] = [misc.translate_exp[b] for b in df["benchmark"]]
    df['platform'] = [misc.translate_platform[p] for p in df['platform']]
    return df 
