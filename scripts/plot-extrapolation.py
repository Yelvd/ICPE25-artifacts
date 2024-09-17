# %% imports
import sys
import math
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from os import walk

from tools import internal
from tools import energy
from tools import misc
from tools import misc

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = snakemake.params["resultsDir"]


# %% loaddata
plt.style.use("../matplotlib-style.rc")
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = "../../ICPP24/results-atomic-block-size/"
exp = None
iterations = 1000

full_system = misc.platform_cores
df = internal.load_time(results_dir, exp)
df = df.reset_index(drop=True)
work = [x.split('-by-') for x in df['atmoic-block'].values]
df['wpn'] = [(int(x[0]) * int(x[1]) * int(x[2])) for x in work]
df['wpn'] = df['wpn'] * df['tasks']
df['iterate-scaled'] = df['iterate'] / df['wpn'] / iterations
df = df.groupby(["atmoic-block", "benchmark", "platform", "tasks", "wpn"], as_index=False).agg({"iterate": ['mean', 'std', 'min', 'max'], "iterate-scaled": ['mean', 'std', 'min', 'max', 'count']})
workunit = 25*25*25*192


# %% extrapolation
weakEfficiancy = lambda n : 0.9801932153425763 + -0.03638144889683877 * np.log2(n)
# weakEfficiancy = lambda n : 1

x = np.arange(1, 512, 2)



figsize = (3, 3)


# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)
plt.plot(x, [weakEfficiancy(s) for s in x])
plt.xlabel("Nodes")
plt.ylabel("Weak-scaling Efficiency")
# plt.xticks(x)
plt.savefig("f{results_dir}/model.pdf")

# %% performance bar plot
new_df = df.loc[df['platform'] == "AMD Genoa"]
xs = np.arange(0, new_df['atmoic-block'].values.size, 1)

figsize = (4, 3)

# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)
plt.errorbar(xs, new_df['iterate']['mean'], yerr=new_df['iterate']['std'])
plt.ylabel("Time per W")
ax.set_xticks(xs, new_df['wpn'] / workunit, rotation=45, ha='right', rotation_mode='anchor')
plt.ylim(0*10**-8)
plt.savefig("f{results_dir}/Time_per_LU.pdf")
# plt.show()



# %% setup work
works = []
for index, row in df.iterrows():

    power = 1
    match row['platform'].values[0]:
        case "AMD Genoa":
            power = 743
        case "AMD Rome 7H12":
            power = 547

    works.append({"atomic": row["atmoic-block"].values[0], "T": row['iterate']['mean'] / iterations, "power": power, "wpn": row['wpn'].values[0], "cpu": row['platform'].values[0]})
        
# %% plot extrapolation time
xs = np.arange(10 * workunit,  128 * workunit, 1 * workunit)

figsize = (4, 3)

# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)
for i, w in enumerate(works):
    if not w['cpu'] == "AMD Genoa":
        continue
    
    nodes = [ math.ceil(x/w['wpn']) for x in xs ]
    if(w['atomic'] == "25-by-25-by-25"):
        time_1 = np.array([(w['T'] / weakEfficiancy(x)) for x in nodes])
        # time_1 = 1

# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)
prev = {}
for i, w in enumerate(works):
    if not w['cpu'] == "AMD Genoa":
        continue

    nodes = [ math.ceil(x/w['wpn']) for x in xs ]
    # nodes = [ (x/w['wpn']) for x in xs ]
    time = np.array([w['T'] / weakEfficiancy(x) for x in nodes])
    
    if(w['atomic'] == "25-by-25-by-25"):
        plt.plot(xs, time / time_1, label=f"wpn: {w['wpn'] / workunit:.0f}")

        # for x, y in zip(xs, time):
            # if int((x / w['wpn'])) in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]:
                # plt.text(x, y, str(int(x / w['wpn'])))
    else:
        if round(w['wpn'] / workunit) in [4, 8, 22, 55]:
            plt.plot(xs, time / time_1, label=f"wpn: {w['wpn'] / workunit:.0f}")


plt.vlines(40 * workunit, 0, 19, linestyle='--', label="Puncture simulation", color='g')
plt.xlabel("Work [$LU^{3}$]")
plt.ylabel("Normalized Time")
plt.legend()
plt.savefig("../artifacts-presentation/time.png", dpi=300)
# plt.show()

# %% plot extrapolation time
xs = np.arange(10 * workunit,  128 * workunit, 1 * workunit)


prev = {}
figsize = (4, 3)

# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)
for i, w in enumerate(works):
    if not w['cpu'] == "AMD Genoa":
        continue
    
    nodes = [ math.ceil(x/w['wpn']) for x in xs ]
    if(w['atomic'] == "25-by-25-by-25"):
        energy_1 = np.array([(w['T'] / weakEfficiancy(x)) * w['power'] * (x) / 1000 for x in nodes])
        # energy_1 = 1

for i, w in enumerate(works):
    if not w['cpu'] == "AMD Genoa":
        continue
    
    nodes = [ math.ceil(x/w['wpn']) for x in xs ]
    # nodes = [ (x/w['wpn']) for x in xs ]
    # print(nodes)
    energy = [(w['T'] / weakEfficiancy(x)) * w['power'] * (x) / 1000 for x in nodes]
    if(w['atomic'] == "25-by-25-by-25"):
        plt.plot(xs, np.array(energy) / energy_1, label=f"wpn: {w['wpn'] / workunit:.0f}")

        # for x, y in zip(xs, energy):
            # if int((x * workunit / w['wpn'])) in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]:
                # plt.text(x-5, y, str(int(x  * workunit/ w['wpn'])))
    else:
        if round(w['wpn'] / workunit) in [4, 8, 22, 55]:
            plt.plot(xs,  np.array(energy) / energy_1, label=f"wpn: {w['wpn'] / workunit:.0f}")

plt.vlines(40 * workunit, 0, 1.9, linestyle='--', label="Puncture simulation", color='g')
plt.xlabel("Work [$LU^3$]")
plt.ylabel("Normalized Energy Cost")
plt.legend()
plt.savefig(f"{results_dir}/energy.pdf")
# plt.show()


# %% performance bar plot
new_df = df.loc[df['platform'] == "AMD Genoa"]
xs = new_df['wpn']
xs2 = np.arange(workunit, max(new_df['wpn'].values), 10)

print(xs2)

figsize = (4, 3)

# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)
plt.plot(xs2, xs2 / workunit, label="Linear scaling")
# plt.errorbar(xs, new_df['iterate']['mean'] / new_df['iterate']['mean'].values[0], yerr=new_df['iterate']['std'] / new_df['iterate']['mean'].values[0], label="Measured")
plt.ylabel("Normalized Time")
plt.xlabel("Work [$LU^3$]")
# ax.set_xticks(xs, new_df['wpn'] / workunit, rotation=45, ha='right', rotation_mode='anchor')
# ax.set_xticks(xs, new_df['wpn'] / workunit, rotation=45, ha='right', rotation_mode='anchor')
plt.legend()
plt.ylim(0*10**-8)
plt.savefig(f"{results_dir}/Time_per_LU.pdf")
# plt.show()
