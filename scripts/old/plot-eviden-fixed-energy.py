import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

CPUs = ["AMD", "ARM", "INTEL", "INTEL_HBM", "Archer2", "ShowCees"]

df = pd.read_csv(snakemake.params["resultsDir"] + "/Eviden-results.csv", sep=';')

# PLOT FIXED
df_fixed = df.loc[df["benchmark"] == "fixed"]
df_fixed["energy"] = np.round((df_fixed['time'] * df_fixed['power']) / 1000).astype(int)
df_fixed_energy = df_fixed[["cpu", "threads", "energy"]]


plt.show()

x = np.arange(len(CPUs))  # the label locations
width = 1 / len(np.unique(np.array(df_fixed_energy['threads']))) # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for threads in np.sort(np.unique(np.array(df_fixed_energy['threads']))):
    offset = width * multiplier
    color = color_cycle[multiplier]
    first = True
    for i, line in df_fixed_energy.loc[df_fixed_energy['threads'] == threads].iterrows():
        if first:
            rects = ax.bar((CPUs.index(line["cpu"]) * 1.5) + offset, line["energy"], width, color=color , label="{}".format(threads))
            first = False
        else:
            rects = ax.bar((CPUs.index(line["cpu"]) * 1.5) + offset, line["energy"], width, color=color)
        ax.bar_label(rects, padding=3)
    multiplier += 1

#Archer
archer = {"Time": 697, "Energy": 308, "Power": 442}
color_i = list(np.sort(np.unique(np.array(df_fixed_energy['threads'])))).index(128)
rects = ax.bar((CPUs.index("Archer2") * 1.5) + .5, archer['Energy'], width, color=color_cycle[color_i])
ax.bar_label(rects, padding=3)

#ShowCees
showcees = {"Time": 983, "Energy": (300*983) / 1000, "Power": 300}
color_i = list(np.sort(np.unique(np.array(df_fixed_energy['threads'])))).index(128)
rects = ax.bar((CPUs.index("ShowCees") * 1.5) + .5, showcees['Energy'], width, color=color_cycle[color_i])
ax.bar_label(rects, padding=3)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Kilo Joules')
ax.set_title('Results for fixed size benchmark')
ax.set_xticks(x*1.5 + .5, CPUs)
# ax.legend(loc='upper left', ncols=3)
ax.legend(ncols=4)
ax.set_ylim(0)

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"], dpi=300)
