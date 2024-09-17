import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import walk

from scripts import tools
from tools import internal
from tools import energy
from tools import misc
from tools import misc

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = snakemake.params["resultsDir"]
exp = snakemake.params["experiment"]
platform = snakemake.params["platform"]
metric = snakemake.params["metric"]
snakemake.params["metric"]


def plot_energy():
    df = energy.load_energy(results_dir)
    df = df.loc[df['platform'] == platform]
    df = df.loc[df['tasks'] == misc.platform_cores[platform] ]
    df = df.loc[df['benchmark'].isin(["S1", "S2", "W1"])]
    df.loc[df['benchmark'] == "S1",'iterate-energy'] = df.loc[df['benchmark'] == "S1"]['iterate-energy'] / 20000 
    df.loc[df['benchmark'] == "S2",'iterate-energy'] = df.loc[df['benchmark'] == "S2"]['iterate-energy'] / 10000 
    df.loc[df['benchmark'] == "W1",'iterate-energy'] = df.loc[df['benchmark'] == "W1"]['iterate-energy'] / 20000 
    df = df.groupby(["platform", "benchmark"], as_index=False).agg({"iterate-energy": ['mean', 'std', 'min', 'max']}).rename({"iterate-energy": "energy"}, axis=1)
    df = df.sort_values(by=["benchmark"])
    df = df.round(0)


    x = np.arange(len(df["benchmark"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        #color = color_cycle[index]
        color = color_cycle[1]
        rects = ax.bar(index, line["energy"]["mean"], 
                        width, color=color, yerr=line["energy"]["std"], label="{}".format(line['benchmark']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Energy [J/ITER.]')

    ax.set_xticks(x, df['benchmark'])
    ax.set_ylim(0, 50)


def plot_power():
    df = energy.load_energy(results_dir)
    df = df.loc[df['platform'] == platform]
    df = df.loc[df['tasks'] == misc.platform_cores[platform] ]
    df = df.loc[df['benchmark'].isin(["S1", "S2", "W1"])]
    df = df.groupby(["platform", "benchmark"], as_index=False).agg({"power": ['mean', 'std']})
    df = df.sort_values(by=["benchmark"])
    df = df.round(0)

    x = np.arange(len(df["benchmark"]))  # the label locations
    width = .8
    width = .8

    for index, line in df.iterrows():
        #color = color_cycle[index]
        color = color_cycle[2]
        rects = ax.bar(index, line["power"]["mean"], 
                        width, color=color, yerr=line["power"]["std"], label="{}".format(line['benchmark']))

        padding = 2
        if index == 1:
            padding = 10

        ax.bar_label(rects, padding=padding)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Power [W]')

    ax.set_xticks(x, df['benchmark'])
    ax.set_ylim(0, 1500)

def plot_time():
    df = internal.load_time(results_dir)
    df = df.loc[df['platform'] == platform]
    df = df.loc[df['tasks'] == misc.platform_cores[platform] ]
    df = df.loc[df['benchmark'].isin(["S1", "S2", "W1"])]
    df.loc[df['benchmark'] == "S1",'iterate'] = df.loc[df['benchmark'] == "S1"]['iterate'] / 20000 * 1000
    df.loc[df['benchmark'] == "S2",'iterate'] = df.loc[df['benchmark'] == "S2"]['iterate'] / 10000 * 1000
    df.loc[df['benchmark'] == "W1",'iterate'] = df.loc[df['benchmark'] == "W1"]['iterate'] / 20000 * 1000
    df = df.groupby(["platform", "benchmark"], as_index=False).agg({"iterate": ['mean', 'std']})
    df = df.sort_values(by=["benchmark"])
    df = df.round(0)
    print(df)

    x = np.arange(len(df["benchmark"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        # color = color_cycle[index]
        color = color_cycle[0]
        rects = ax.bar(index, line["iterate"]["mean"], 
                        width, yerr=line["iterate"]["std"], color=color, label="{}".format(line['benchmark']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Time [ms/ITER.]')

    ax.set_xticks(x, df['benchmark'])
    ax.set_ylim(0, 50)


figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0] / 3, figsize[0] / 3)

fig, ax = plt.subplots(figsize=figsize, layout='constrained')

match metric:
    case "time":
        plot_time()
    case "power":
        plot_power()
    case "energy":
        plot_energy()

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])
