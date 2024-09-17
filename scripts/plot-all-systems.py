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

full_system = misc.platform_cores

systems = ["AMD Bergamo", "Intel SPR DDR", "Intel SPR HBM", "Ampere Altra", "AMD Genoa"]
systems = ["AMD Bergamo", "Intel SPR DDR", "Intel SPR HBM", "Ampere Altra", "AMD Genoa", "AMD Rome 7742", "AMD Rome 7H12"]

def plot_energy(ax):
    df = energy.load_energy(results_dir, exp)
    df = df.loc[df['platform'].isin(systems)]
    df = df.reset_index(drop=True)
    df['iterate-energy'] = df['iterate-energy'] / 20000 

    # Filter to plot only for experiments wher the number of tasks is the number of cores
    for p in np.unique(df['platform']):
        index_not_full = df[(df['platform'] == p) & (df['tasks'] != full_system[p])].index
        df = df.drop(index_not_full)
 
    df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate-energy": ['mean', 'std', 'min', 'max']}).rename({"iterate-energy": "energy"}, axis=1)
    df = df.sort_values(by=["platform"])
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
    ax.set_ylabel('Energy [J/ITER.]')

    ax.set_xticks(x, df['platform'], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 70)


def plot_power(ax):
    df = energy.load_energy(results_dir, exp)
    df = df.loc[df['platform'].isin(systems)]
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
    ax.set_ylabel('Power [W]')

    ax.set_xticks(x, df['platform'], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 1400)

def plot_time(ax):
    df = internal.load_time(results_dir, exp)
    df = df.loc[df['platform'].isin(systems)]
    df = df.reset_index(drop=True)

    # Filter to plot only for experiments wher the number of tasks is the number of cores
    for p in np.unique(df['platform']):
        index_not_full = df[(df['platform'] == p) & (df['tasks'] != full_system[p])].index
        df = df.drop(index_not_full)

    df['iterate'] = df['iterate'] / 20000 * 1000
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
    ax.set_ylabel('Time [ms/ITER.]')

    #ax.set_xticks(x, df['platform'], rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 150)


figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0], figsize[0] + 0.5)

fig, ax = plt.subplots(3, figsize=figsize, layout='constrained', sharex=True)
# fig, ax = plt.subplots(figsize=figsize, layout='constrained')

plot_time(ax[0])
plot_power(ax[1])
plot_energy(ax[2])

fig.align_ylabels()

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])
