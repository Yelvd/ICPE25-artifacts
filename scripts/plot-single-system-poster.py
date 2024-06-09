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


def plot_energy(ax):
    df = energy.load_energy(results_dir, exp)
    df = df.loc[df['platform'] == platform]
    df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate-energy": ['mean', 'std', 'min', 'max']}).rename({"iterate-energy": "energy"}, axis=1)
    df = df.sort_values(by=["tasks"])
    df['energy'] = df['energy'] / 1000
    df = df.round(0)


    x = np.arange(len(df["tasks"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        #color = color_cycle[index]
        color = color_cycle[1]
        rects = ax.bar(index, line["energy"]["mean"], 
                        width, color=color, label="{}".format(line['tasks']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Energy (kJ)')
    ax.set_xlabel('Processes')

    ax.set_xticks(x, df['tasks'])
    ax.set_ylim(0, 1200)


def plot_power(ax):
    df = energy.load_energy(results_dir, exp)
    df = df.loc[df['platform'] == platform]
    df = df.groupby(["platform", "tasks"], as_index=False).agg({"power": ['mean', 'std']})
    df = df.sort_values(by=["tasks"])
    df = df.round(0)

    x = np.arange(len(df["tasks"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        #color = color_cycle[index]
        color = color_cycle[2]
        rects = ax.bar(index, line["power"]["mean"], 
                        width, color=color, label="{}".format(line['tasks']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Power (W)')

    ax.set_xticks(x, df['tasks'])
    ax.set_ylim(0, 1300)

def plot_time(ax):
    df = internal.load_time(results_dir, exp)
    df = df.loc[df['platform'] == platform]
    df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate": ['mean', 'std']})
    df = df.sort_values(by=["tasks"])
    df = df.round(0)

    x = np.arange(len(df["tasks"]))  # the label locations
    width = .8

    for index, line in df.iterrows():
        # color = color_cycle[index]
        color = color_cycle[0]
        rects = ax.bar(index, line["iterate"]["mean"], 
                        width, color=color, label="{}".format(line['tasks']))

        ax.bar_label(rects, padding=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Time (s)')

    ax.set_xticks(x, df['tasks'])
    ax.set_ylim(0, 3100)


figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0], figsize[1])

fig, ax = plt.subplots(3, figsize=figsize, layout='constrained', sharex=True)

plot_time(ax[0])
plot_power(ax[1])
plot_energy(ax[2])

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])
