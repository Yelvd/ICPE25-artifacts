import sys
import math
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

systems = ["AMD Bergamo", "Intel SPR DDR", "Intel SPR HBM", "Ampere Altra", "AMD Genoa", "AMD Rome 7742", "AMD Rome 7H12"]
systems = ["AMD Bergamo", "Intel SPR DDR", "Intel SPR HBM", "Ampere Altra", "AMD Genoa", "AMD Rome 7H12"]

def plot_energy(ax):
    df = energy.load_energy(results_dir, exp)
    df = df.loc[df['platform'] == platform]
    if platform == "AMD Genoa":
        df = df.loc[df['tasks'] != 32]
    df['iterate-energy'] = df['iterate-energy'] / 20000 
    df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate-energy": ['mean', 'std', 'min', 'max']}).rename({"iterate-energy": "energy"}, axis=1)
    df = df.sort_values(by=["tasks"])
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
    ax.set_title(f"{platform}")
    ax.set_xlabel('Processes')

    ax.set_xticks(x, df['tasks'])
    ax.set_ylim(0, 60)


def plot_power(ax):
    df = energy.load_energy(results_dir, exp)
    df = df.loc[df['platform'] == platform]
    if platform == "AMD Genoa":
        df = df.loc[df['tasks'] != 32]
    df = df.groupby(["platform", "tasks"], as_index=False).agg({"power": ['mean', 'std']})
    df = df.sort_values(by=["tasks"])
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
    ax.set_title(f"{platform}")
    ax.set_xlabel('Processes')

    ax.set_xticks(x, df['tasks'])
    ax.set_ylim(0, 1100)

def plot_time(ax):
    df = internal.load_time(results_dir, exp)
    df = df.loc[df['platform'] == platform]
    df['iterate'] = df['iterate'] / 20000 * 1000
    if platform == "AMD Genoa":
        df = df.loc[df['tasks'] != 32]
    
    df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate": ['mean', 'std']})
    df = df.sort_values(by=["tasks"])
    df = df.round(0)

    x = np.arange(len(df["tasks"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        # color = color_cycle[index]
        color = color_cycle[0]
        rects = ax.bar(index, line["iterate"]["mean"], 
                        width, yerr=line["iterate"]["std"], color=color, label="{}".format(line['tasks']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_title(f"{platform}")
    ax.set_xlabel('Processes')

    ax.set_xticks(x, df['tasks'])
    ax.set_ylim(0, 150)



figsize = plt.rcParams['figure.figsize']

if not platform == "all":
    figsize = (figsize[0] / 2, figsize[0] / 2)
    fig, ax = plt.subplots(figsize=figsize, layout='constrained')

    match metric:
        case "time":
            ax.set_ylabel('Time [ms/ITER.]')
            plot_time(ax)
        case "power":
            ax.set_ylabel('Power [W]')
            plot_power(ax)
        case "energy":
            ax.set_ylabel('Energy [J/ITER.]')
            plot_energy(ax)

    plt.savefig(snakemake.output["png"], dpi=300)
    plt.savefig(snakemake.output["pdf"])



if platform == "all":
    figsize = (figsize[0], figsize[1])
    ncolumns = 3
    fig, ax = plt.subplots(math.ceil(len(systems) / ncolumns), ncolumns, sharey=True,figsize=figsize, layout='constrained')

    flat_ax = ax.flatten()
    for i, platform in enumerate(systems):
        match metric:
            case "time":
                if i % ncolumns == 0:
                    flat_ax[i].set_ylabel('Time [ms/ITER.]')
                plot_time(flat_ax[i])
            case "power":
                if i % ncolumns == 0:
                    flat_ax[i].set_ylabel('Power [W]')
                plot_power(flat_ax[i])
            case "energy":
                if i % ncolumns == 0:
                    flat_ax[i].set_ylabel('Energy [J/ITER.]')
                plot_energy(flat_ax[i])

    plt.savefig(snakemake.output["png"], dpi=300)
    plt.savefig(snakemake.output["pdf"])
