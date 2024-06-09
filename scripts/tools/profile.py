import sys
import os
import pandas as pd
import numpy as np
import yaml
import sys
from . import misc

from enum import Enum
from os import listdir
from os.path import isfile, join
from os import walk

class DATA_TYPE(Enum):
    METRIC = 'Metric'
    RAW = 'Raw'

PROFILE_NAMES = ["likwid", "perf", "profile"]

def load_snellius_likwid_out(experiment_dir, meta, metrics):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    files = [f for f in files if "FLOPS_DP" in f]
    files = [f for f in files if "out" in f]
    files = [f for f in files if not f.startswith('.')]
    results = {}

    results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])
    results["id"] = meta["General"]["Id"]

    data = {}
    groups = 0
    for m in metrics:
        data[m] = [0] * results['tasks']

    for filename in files:
        taskId = 1

        with open(f"{experiment_dir}/{filename}") as f:
            for l in f:
                if "HemoCell:" in l:
                    data["hemocell"] = float(l.split(" ")[-1])
                if "iterate:" in l:
                    data["iterate"] = float(l.split(" ")[-1])
                if "Group" in l:
                    g = int(l.split(' ')[1][:-1])
                    groups = g if g > groups else groups
                for m in metrics:
                    if m in l:
                        data[m] = [float(x) for x in l.strip().split('|')[3:] if x != '' ]

    # Interpolation
    for m in metrics:
        data[m] = [d * groups for d in data[m]]

    columns = ["benchmark", "platform", "id", "tasks", "metric", "min", "max", "avg", "std", "sum"]

    new_results = []
    for m in metrics + ["hemocell", "iterate"]:
        new_results.append([ results['benchmark'], 
                results['platform'], 
                results['id'],
                results['tasks'], 
                m,
                np.min(data[m]),
                np.max(data[m]),
                np.mean(data[m]),
                np.std(data[m]),
                np.sum(data[m])])

    df = pd.DataFrame(new_results, columns=columns)

    return df

def load_eviden_likwid(experiment_dir, meta, metrics):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    files = [f for f in files if "likwid" in f]
    files = [f for f in files if "txt" in f]
    files = [f for f in files if not f.startswith('.')]
    results = {}

    results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])
    results["id"] = meta["General"]["Id"]

    data = {}
    groups = 0
    for m in metrics:
        data[m] = [0] * results['tasks']

    for filename in files:
        taskId = int(filename.split('.')[0].split('_')[-1])

        with open(f"{experiment_dir}/{filename}") as f:
            for l in f:
                if "Group" in l:
                    g = int(l.split(' ')[1][:-1])
                    groups = g if g > groups else groups
                for m in metrics:
                    if m in l:
                        data[m][taskId] = float(l.strip().split('|')[-2].strip())

    for m in metrics:
        data[m] = [d * groups for d in data[m]]
    columns = ["benchmark", "platform", "id", "tasks", "metric", "min", "max", "avg", "std", "sum"]

    new_results = []
    for m in metrics:
        new_results.append([ results['benchmark'], 
                results['platform'], 
                results['id'],
                results['tasks'], 
                m,
                np.min(data[m]),
                np.max(data[m]),
                np.mean(data[m]),
                np.std(data[m]),
                np.sum(data[m])])

    df = pd.DataFrame(new_results, columns=columns)

    return df

def load_snellius_rome(experiment_dir, meta, metrics):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    files = [f for f in files if "likwid" in f]
    files = [f for f in files if "csv" in f]
    files = [f for f in files if not f.startswith('.')]
    results = {}

    if meta["General"].get("Benchmark") is None:
        print("WARNING FIXING TO FIXED")
        results["benchmark"] = "fixed"
    else:
        results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])
    results["id"] = meta["General"]["Id"]

    data = {}
    for m in metrics:
        data[m] = [0] * results['tasks']

    for filename in files:
        taskId = int(filename.split('.')[0].split('-')[-1])

        table = False 
        with open(f"{experiment_dir}/{filename}") as f:
            for l in f:
                if "TABLE" in l:
                    if "Raw" in l:
                        table = 1
                    if "Metric" in l:
                        table = 2
                    continue

                if table: 
                    for m in metrics:
                        if m in l:
                            if table == 1:
                                data[m][taskId] = float(l.strip().split(',')[-1].strip())
                            if table == 2:
                                data[m][taskId] = float(l.strip().split(',')[-2].strip())


    columns = ["benchmark", "platform", "id", "tasks", "metric", "min", "max", "avg", "std", "sum"]

    new_results = []
    for m in metrics:
        new_results.append([ results['benchmark'], 
                results['platform'], 
                results['id'],
                results['tasks'], 
                m,
                np.min(data[m]),
                np.max(data[m]),
                np.mean(data[m]),
                np.std(data[m]),
                np.sum(data[m])])

    df = pd.DataFrame(new_results, columns=columns)

    return df


def load_perf_eviden(experiment_dir, meta, metrics):
    """
    Load perf data from eviden experiment.
    """
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    files = [f for f in files if "perf_report.log" in f]
    files = [f for f in files if not f.startswith('.')]
    results = {}

    results["benchmark"] = meta["General"]["Benchmark"]
    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])
    results["id"] = meta["General"]["Id"]

    metrics = metrics + ["runtime"]

    data = {}
    for m in metrics:
        data[m] = [0]

    for filename in files:
        taskId = 0

        with open(f"{experiment_dir}/{filename}") as f:
            for l in f:
                for m in metrics:
                    if m in l:
                        data[m][taskId] = float(l.strip().split(' ')[0].replace(',' , ''))
                    if "time elapsed" in l:
                        data["runtime"][taskId] = float(l.strip().split(' ')[0].replace(',' , ''))


    columns = ["benchmark", "platform", "id", "tasks", "metric", "min", "max", "avg", "std", "sum"]

    new_results = []
    for m in metrics:
        new_results.append([ results['benchmark'], 
                results['platform'], 
                results['id'],
                results['tasks'], 
                m,
                np.min(data[m]),
                np.max(data[m]),
                np.mean(data[m]),
                np.std(data[m]),
                np.sum(data[m])])

    df = pd.DataFrame(new_results, columns=columns)

    return df


def load_perf(experiment_dir, meta, metrics):
    """
    Loading Perf data where every task writes its own perf.data file
    """
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    files = [f for f in files if "perf.data" in f]
    files = [f for f in files if not f.startswith('.')]
    results = {}

    if meta["General"].get("Benchmark") is None:
        print("WARNING FIXING TO FIXED")
        results["benchmark"] = "fixed"
    else:
        results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])
    results["id"] = meta["General"]["Id"]

    data = {}
    for m in metrics:
        data[m] = [0] * results['tasks']

    for filename in files:
        taskId = int(filename.split('.')[-1])

        with open(f"{experiment_dir}/{filename}") as f:
            for l in f:
                for m in metrics:
                    if m in l:
                        data[m][taskId] = float(l.strip().split(' ')[0].replace(',' , ''))

    columns = ["benchmark", "platform", "id", "tasks", "metric", "min", "max", "avg", "std", "sum"]

    new_results = []
    for m in metrics:
        new_results.append([ results['benchmark'], 
                results['platform'], 
                results['id'],
                results['tasks'], 
                m,
                np.min(data[m]),
                np.max(data[m]),
                np.mean(data[m]),
                np.std(data[m]),
                np.sum(data[m])])

    df = pd.DataFrame(new_results, columns=columns)

    return df



def load_experiment(experiment_dir, meta, metrics):

    if not metrics.get("snellius-rome-likwid"):
        if not metrics.get(meta["General"]["Platform"]):
            return

    match meta["General"]["Platform"]:
        case "snellius-rome":
            if meta['General']["Monitor-tool"] == "likwid":
                if "snellius-rome-likwid" in metrics:
                    return load_snellius_likwid_out(experiment_dir, meta, metrics["snellius-rome-likwid"])
                else:
                    return
            return load_perf(experiment_dir, meta, metrics["snellius-rome"])
        case "snellius-genoa":
            if meta['General']["Monitor-tool"] == "likwid":
                return
            return load_perf(experiment_dir, meta, metrics["snellius-genoa"])
        case "archer_rome":
            return load_perf(experiment_dir, meta, metrics["archer_rome"])

    if meta["General"]["Platform"] in ["Intel_SPR_DDR", "Intel_SPR_HBM", "Ampere_Q8030", "AMD_Epyc"]:
        if meta["General"]["Monitor-tool"] == "perf":
            return load_perf_eviden(experiment_dir, meta, metrics[meta["General"]["Platform"]])
        if meta["General"]["Monitor-tool"] in ["profile", "likwid"]:
            return load_eviden_likwid(experiment_dir, meta, metrics[meta["General"]["Platform"]])


def load_profile(results_dir, metrics, experiment=None, tool=None):
    """ Load all the hardware counter data from the results_dir
    metircs: {platform_name: [list of hw counters to find]
    exp: filter based on experiment name
    tool: filter based on tool name
    """

    # Translate from external to internal exp names
    experiment = misc.translate_exp_to_internal[experiment]

    # Translate from external to internal patform names
    metrics = metrics.copy()
    keys = list(metrics.keys())
    for k in keys:
        metrics[misc.translate_platform_rev[k]] = metrics[k]
    data = []

    w = walk(results_dir)
    for (dirpath, dirnames, filenames) in w:
        if "meta.yaml" not in filenames:
            continue


        with open(f"{dirpath}/meta.yaml") as f:
            meta = yaml.load(f, Loader=yaml.FullLoader)


        if meta["General"].get("Monitor-tool") is None:
            continue
            if meta["Instrumentation"]["Collection tool"] in PROFILE_NAMES:
                data.append(load_experiment(dirpath, meta, metrics))

        if meta["General"]["Monitor-tool"] in PROFILE_NAMES:
            if tool and not meta["General"]["Monitor-tool"] == tool:
                continue
            data.append(load_experiment(dirpath, meta, metrics))

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
