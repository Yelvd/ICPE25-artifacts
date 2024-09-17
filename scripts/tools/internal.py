import sys
import os
import pandas as pd
import numpy as np
import yaml
import sys
from . import misc
from os import listdir
from os.path import isfile, join
from os import walk

def load_time(results_dir, experiment=None):
    data = []
    # experiment = misc.translate_exp_to_internal[experiment]
    w = walk(results_dir)
    for (dirpath, dirnames, filenames) in w:
        if "meta.yaml" not in filenames:
            continue

        with open(f"{dirpath}/meta.yaml") as f:
            meta = yaml.load(f, Loader=yaml.FullLoader)

        if meta["General"].get("Monitor-tool") is None:
            if meta["Instrumentation"]["Collection tool"] == "internal":
                data.append(load_experiment(dirpath))
            continue

        if meta["General"]["Monitor-tool"] in ["time", "internal"]:
            data.append(load_experiment(dirpath))


    df = pd.concat(data)

    if experiment is not None:
        df = df.loc[df["benchmark"] == experiment]
        # match experiment:
            # case "fixed":
                # df = df.loc[df["benchmark"].isin(["fixed", "cube-fixed"])]
            # case "fixed-large":
                # df = df.loc[df["benchmark"].isin(["fixed-large", "cube-fixed-large"])]
            # case "relative":
                # df = df.loc[df["benchmark"].isin(["relative", "cube-relative"])]

    return df 


def load_experiment(experiment_dir):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    
    results = {}
    for f in files:
        if ".out" in f and "bull" not in f and "post" not in f:
            out_filename = f

    with open(f"{experiment_dir}/{out_filename}") as f:
        lines = [line.rstrip() for line in f]

        for l in lines:
            if "HemoCell:" in l:
                results["hemocell"] = float(l.split(" ")[-1])
            if "iterate:" in l:
                results["iterate"] = float(l.split(" ")[-1])
            if "Smallest atomic-block:" in l:
                results["atmoic-block"] = l.split(" ")[-1]
    
    with open(f"{experiment_dir}/meta.yaml") as f:
        meta = yaml.load(f, Loader=yaml.FullLoader)

    if meta["General"].get("Benchmark") is None:
        print("WARNING FIXING TO FIXED")
        print(experiment_dir)
        results["benchmark"] = "fixed"
    else:
        results["benchmark"] = meta["General"]["Benchmark"]
    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])


    for k in results.keys():
        results[k] = [results[k]]

    df = pd.DataFrame.from_dict(results)
    df["benchmark"] = [misc.translate_exp[b] for b in df["benchmark"]]
    df['platform'] = [misc.translate_platform[p] for p in df['platform']]
    return df
