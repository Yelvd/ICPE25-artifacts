import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import walk

sys.path.append('./scripts/tools')
import internal

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = snakemake.params["resultsDir"]

df = internal.load_time(results_dir, "fixed")
df = df.groupby(["platform", "tasks"], as_index=False).agg({"iterate": ['mean', 'std']}).round(2)

platforms = list(np.unique(df['platform']))

x = np.arange(len(platforms))  # the label locations
width = 1 / len(np.unique(np.array(df['tasks']))) # the width of the bars
multiplier = 0

fig, ax = plt.subplots(figsize=(10, 5), layout='constrained')

plt.plot([], [], ' ', label="Processes:")
for tasks in np.sort(np.unique(np.array(df['tasks']))):
    offset = width * multiplier
    color = color_cycle[multiplier]
    first = True
    for i, line in df.loc[df['tasks'] == tasks].iterrows():
        if first:
            rects = ax.bar((platforms.index(np.array(line["platform"])[0]) * 1.5) + offset, line["iterate"]["mean"], 
                            width, color=color, yerr=line["iterate"]["std"], label="{}".format(tasks))
            first = False
        else:
            rects = ax.bar((platforms.index(np.array(line["platform"])[0]) * 1.5) + offset, line["iterate"]["mean"], 
                            width, yerr=line["iterate"]["std"], color=color)

        ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('seconds')
ax.set_title('Results for fixed size benchmark')
ax.set_xticks(x*1.5 + .5, platforms)
# ax.legend(loc='upper left', ncols=3)
ax.legend(ncol=4)
ax.set_ylim(0)

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])

with open(snakemake.output["tex"], 'w') as tf:
     tf.write(df.style.to_latex())
