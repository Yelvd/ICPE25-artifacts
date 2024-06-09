import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from os import walk

from scripts import tools
from tools import internal
from tools import energy
from tools import misc
from tools import misc

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = snakemake.params["resultsDir"]
exp = None

full_system = misc.platform_cores


df = internal.load_time(results_dir, exp)
df = df.reset_index(drop=True)

df.loc[df['platform'] == "AMD Genoa", 'nodes'] = df.loc[df['platform'] == "AMD Genoa"]['tasks'] / 192
df.loc[df['platform'] == "AMD Rome 7H12", 'nodes'] = df.loc[df['platform'] == "AMD Rome 7H12"]['tasks'] / 128
df.loc[df['platform'] == "AMD Genoa", 'efficiancy'] = np.mean(df.loc[((df['platform'] == "AMD Genoa") & (df['nodes'] == 1))]['iterate'].values) / df.loc[df['platform'] == "AMD Genoa"]['iterate'] 
df.loc[df['platform'] == "AMD Rome 7H12", 'efficiancy'] = np.mean(df.loc[((df['platform'] == "AMD Rome 7H12") & (df['nodes'] == 1))]['iterate'].values) / df.loc[df['platform'] == "AMD Rome 7H12"]['iterate'] 
# df = df.groupby(["benchmark", "platform", "tasks", "nodes"]).agg({'iterate':['mean','std']})

# coding: utf-8
def create_json_extraP(df, params, component, param_names=None, metric='iterate', filename=None):
    if filename is None:
        filename = "test.json"
        
    if param_names is None:
        param_names = params
    
    output={}
    output["parameters"] = param_names
    output["measurements"] = {}
    output["measurements"][component] = {}
    output["measurements"][component]["time"] = []
    
    measurements = []

    for n in np.unique(df["nodes"]):
        measurements.append({"point": [n], "values": list(df.loc[df['nodes'] == n][component])})

    
    output['measurements'][component]['time'] = measurements
    json_object = json.dumps(output, indent=4)
     
    with open(filename, "w") as outfile:
        outfile.write(json_object)


df = df.loc[df['atmoic-block'] == "25-by-25-by-25"]
print(df)
tmp_df = df.groupby(["platform", "nodes", "benchmark"]).agg({"iterate":['mean', 'std', 'min', 'max'], "efficiancy":['mean', 'std', 'min', 'max']})
print(tmp_df)
create_json_extraP(df.loc[df["platform"] != "AMD Genoa"], ["nodes"], "efficiancy")




# figsize = plt.rcParams['figure.figsize']
# figsize = (figsize[0], figsize[0] / 1.5)

# fig, ax = plt.subplots(2, figsize=figsize, layout='constrained', sharex=True)
# # fig, ax = plt.subplots(figsize=figsize, layout='constrained')

# plot_time(ax[0])
# plot_energy(ax[1])
# # match metric:
    # # case "time":
        # # plot_time()
    # # case "power":
        # # plot_power()
    # # case "energy":
        # # plot_energy()

# plt.savefig(snakemake.output["png"], dpi=300)
# plt.savefig(snakemake.output["pdf"])

def plot_energy(ax):
    df = energy.load_energy(results_dir, exp)
    df = df.reset_index(drop=True)

    # Filter to plot only for experiments wher the number of tasks is the number of cores
    for p in np.unique(df['platform']):
        index_not_full = df[(df['platform'] == p) & (df['tasks'] != full_system[p])].index
        df = df.drop(index_not_full)

    df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate-energy": ['mean', 'std', 'min', 'max']}).rename({"iterate-energy": "energy"}, axis=1)
    df = df.sort_values(by=["platform"])
    df['energy'] = df['energy'] / 1000
    df = df.round(0)


    x = np.arange(len(df["tasks"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        #color = color_cycle[index]
        color = color_cycle[1]
        rects = ax.bar(index, line["energy"]["mean"], 
                        width, color=color, yerr=line["energy"]["std"], label="{}".format(line['tasks']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Energy (kJ)')

    ax.set_xticks(x, df['platform'], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 1100)


def plot_power():
    df = energy.load_energy(results_dir, exp)
    df = df.reset_index(drop=True)

    # Filter to plot only for experiments wher the number of tasks is the number of cores
    for p in np.unique(df['platform']):
        index_not_full = df[(df['platform'] == p) & (df['tasks'] != full_system[p])].index
        df = df.drop(index_not_full)

    df = df.groupby(["platform", "tasks"], as_index=False).agg({"power": ['mean', 'std']})
    df = df.sort_values(by=["platform"])
    df = df.round(0)

    x = np.arange(len(df["tasks"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        #color = color_cycle[index]
        color = color_cycle[2]
        rects = ax.bar(index, line["power"]["mean"], 
                        width, color=color, yerr=line["power"]["std"], label="{}".format(line['tasks']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Avg power (W)')

    ax.set_xticks(x, df['platform'], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 1200)

def plot_time(ax):
    df = internal.load_time(results_dir, exp)
    df = df.reset_index(drop=True)

    # Filter to plot only for experiments wher the number of tasks is the number of cores
    for p in np.unique(df['platform']):
        index_not_full = df[(df['platform'] == p) & (df['tasks'] != full_system[p])].index
        df = df.drop(index_not_full)

    df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate": ['mean', 'std', 'min', 'max']})
    df = df.sort_values(by=["platform"])
    df = df.round(0)

    x = np.arange(len(df["platform"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        # color = color_cycle[index]
        color = color_cycle[0]
        rects = ax.bar(index, line["iterate"]["mean"], 
                        width, color=color, yerr=line["iterate"]["std"], label="{}".format(line['tasks']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Time (s)')

    #ax.set_xticks(x, df['platform'], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 2400)
