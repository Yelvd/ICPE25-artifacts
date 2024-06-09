import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import walk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from matplotlib.ticker import FormatStrFormatter

from scripts import tools
from tools import internal
from tools import energy
from tools import misc
from tools import misc
from tools import profile
from tools import scorep

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = snakemake.params["resultsDir"]
full_system = misc.platform_cores

names = misc.names_for_hw_events
translate = misc.translate_hw_events

roofline_metrics = ["RETIRED_SSE_AVX_FLOPS_ALL |"]
metrics = {"rome-likwid": roofline_metrics}
df = profile.load_profile(results_dir, metrics, None, "likwid")
df = df.sort_values(by=['platform', 'tasks', 'benchmark'])
df['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}-{r['tasks']}" for _, r in df.iterrows()]

df_score = scorep.load_scorep(results_dir, None)
df_score = df_score.loc[df_score["region"].isin(["fluid", "particle", "couple"])]

df_score = df_score.sort_values(by=['platform', 'tasks', 'benchmark'])
f = scorep.load_scorep("./results-components/")
df_score.loc[df_score["region"].isin(["fluid", "particle", "couple"])]
df_score['metric'] = [translate[m.strip()] for m in df_score['metric']]

platform = "AMD Genoa"
df_time = internal.load_time("./results")
df_time = df_time.loc[df_time['tasks'] == misc.platform_cores[platform] ]
df_time = df_time.loc[df_time['benchmark'].isin(["S1", "S2", "W1"])]
df_time = df_time.groupby(["platform", "benchmark"], as_index=False).agg({"iterate": ['mean', 'std']})
df_time = df_time.sort_values(by=["benchmark"])
df_time = df_time.round(0)



df_score['unique'] = [f"{r['benchmark']}-{r['platform']}-{r['id']}" for _, r in df_score.iterrows()]

def plot_roofline(name, bandwidth, flops, points):
    xRange = [0.01, 10]
    yRange = (10, 10000)

    peakFlops = np.max(list(flops.values()))
    peakBw    = np.max(list(bandwidth.values()))
    linewidth = 1

    plt.style.use('./matplotlib-style.rc')
    color_cycle = matplotlib.pyplot.rcParams['axes.prop_cycle'].by_key()['color']
    linestyles = ['-', ':']


    for i, f in enumerate(flops.keys()):
        FLOPS = flops[f]
        flopIntersect = FLOPS / peakBw
        plt.plot([flopIntersect, xRange[1]], [FLOPS, FLOPS], linestyles[0], linewidth=linewidth,color=color_cycle[i], label="{}: {} GFLOP/s".format(f, FLOPS))

    # Plot slanted bandwidht lines
    iOffset = len(flops.keys())
    for i, b in enumerate(bandwidth.keys()):
        bw = bandwidth[b]
        flopIntersect = peakFlops / bw
        plt.plot([0, flopIntersect], [0, bw * flopIntersect], linestyles[1], linewidth=linewidth,color=color_cycle[i + iOffset], label="{}: {} GB/s".format(b, bw))

    
    iOffset = len(flops.keys()) + len(bandwidth.keys())


    for i, p in enumerate(points.keys()):
        marker='x' if "fluid" in p else '+'
        
        tmp = 9 if "L1" in p else 7
        color = color_cycle[tmp]
        plt.plot(points[p][0], points[p][1], color=color, marker=marker, linestyle='None', markersize=5, label="{}: {:.2f} FLOP/Byte".format(p, points[p][0]))


    ax=plt.gca()
    plt.grid(True, which="both",linestyle = '--', linewidth = 0.2)
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("GFLOP/s")
    plt.xlabel("FLOP/Byte")
    plt.ylim(yRange[0], yRange[1])
    plt.xlim(xRange[0], xRange[1])
    plt.xticks([0.01, 0.1, 1, 10])
    plt.yticks([10, 100, 1000, 10000])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.50),
          ncol=2, fancybox=True, shadow=True)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1g'))



new_columns = ["Benchmark", "Platform", "Processes", "Time", "FLOPS"]
new_df = []


# Extract FLOP per iteraion from AMD rome snellius likwid results
for job in np.unique(df['unique']):
    tmp_df = df.loc[df['unique'] == job]

    r = tmp_df.loc[tmp_df['metric'] == "iterate"] 
    time = tmp_df.loc[tmp_df['metric'] == "iterate"]['avg'].values[0]
    platform = r['platform'].values[0]
    benchmark = r['benchmark'].values[0]
    tasks = int(r['tasks'].values[0])

    flop = tmp_df.loc[tmp_df['metric'] == "RETIRED_SSE_AVX_FLOPS_ALL |"]['sum'].values[0]
    
    if tasks == 1:
        flop = flop / 100
        time = time / 100
    # else:
        # if "S2" in benchmark:
            # flop = flop / 10000
        # else:
            # flop = flop / 20000

    new_df.append([benchmark, platform, tasks, time, flop])

new_df = pd.DataFrame(new_df, columns=new_columns)
flop_df = new_df.groupby(['Platform', 'Processes', "Benchmark"]).mean().reset_index()
flop_df = flop_df.loc[flop_df['Processes'] == 128]

new_columns = ["Benchmark", "Platform", "Region", "Time", "FLOP/s", "FLOPS/Byte", "FLOPS/Byte-miss"]
# Calculate FLOPS/s and FLOPS/byte
new_df = []
for job in np.unique(df_score['unique']):
    tmp_df = df_score.loc[df_score['unique'] == job]
    tmp_f = tmp_df.loc[tmp_df['region'] == "fluid"]
    tmp_p = tmp_df.loc[tmp_df['region'].isin(["particle", "couple"])]

    r = tmp_df.loc[tmp_df['metric'] == "total"] 
    platform = r['platform'].values[0]
    benchmark = r['benchmark'].values[0]
    tasks = int(r['tasks'].values[0])
    
    # Time spend in each kernel
    total_time = flop_df.loc[flop_df['Benchmark'] == f"{benchmark}"]["Time"].values[0]
    f_time = np.sum(tmp_f.loc[tmp_f['metric'] == "total"]["mean"])
    p_time = np.sum(tmp_p.loc[tmp_p['metric'] == "total"]["mean"]) 
    f_ratio = (f_time / total_time)
    p_ratio = (p_time / total_time)
    total_time = df_time.loc[df_time["benchmark"] == benchmark]["iterate"]["mean"].values[0]
    f_time = df_time.loc[df_time["benchmark"] == benchmark]["iterate"]["mean"].values[0] * f_ratio
    p_time = df_time.loc[df_time["benchmark"] == benchmark]["iterate"]["mean"].values[0] * p_ratio


    # Total bytes of Data requested load/store
    total_req = np.sum(tmp_df.loc[tmp_df['metric'] == names['l1_accesses']]["sum"]) * 8
    f_req = np.sum(tmp_f.loc[tmp_f['metric'] == names['l1_accesses']]["sum"]) * 8
    p_req = np.sum(tmp_p.loc[tmp_p['metric'] == names['l1_accesses']]["sum"]) * 8 

    total_miss = np.sum(tmp_df.loc[tmp_df['metric'] == names['l1_misses']]["sum"]) * 8
    f_miss = np.sum(tmp_f.loc[tmp_f['metric'] == names['l1_misses']]["sum"]) * 8
    p_miss = np.sum(tmp_p.loc[tmp_p['metric'] == names['l1_misses']]["sum"]) * 8 

    # Total FLOPS for both kernels
    f_floppi = flop_df.loc[flop_df['Benchmark'] == f"{benchmark}-no-RBC"]["FLOPS"].values[0]
    p_floppi = flop_df.loc[flop_df['Benchmark'] == f"{benchmark}"]["FLOPS"].values[0] - f_floppi

    # if benchmark == "S2":
        # f_flop = f_floppi * 10000
        # p_flop = p_floppi * 10000
    # else:
        # f_flop = f_floppi * 20000
        # p_flop = p_floppi * 20000
    f_flop = f_floppi
    p_flop = p_floppi

    
    new_df.append([benchmark, platform, "full"    , total_time, (f_flop + p_flop) / total_time, (f_flop + p_flop)/total_req, (f_flop + p_flop)/total_miss])
    new_df.append([benchmark, platform, "fluid"   ,     f_time,                f_flop / f_time,                f_flop/f_req,                f_flop/f_miss])
    new_df.append([benchmark, platform, "particle",     p_time,                p_flop / p_time,                p_flop/p_req,                p_flop/p_miss])

new_df = pd.DataFrame(new_df, columns=new_columns)
new_df = new_df.loc[new_df["Platform"] == "AMD Genoa"]
new_df = new_df.groupby(['Platform', "Benchmark", "Region"]).mean().reset_index()



figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0], figsize[1])


# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)


kernels = {
        "fluid-L1": (new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "fluid"))]['FLOPS/Byte'].values[0], 
                  new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "fluid"))]['FLOP/s'].values[0] / 10**9),
        "particle-L1": (new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "particle"))]['FLOPS/Byte'].values[0], 
                     new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "particle"))]['FLOP/s'].values[0] / 10**9),

        "fluid-L2": (new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "fluid"))]['FLOPS/Byte-miss'].values[0], 
                  new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "fluid"))]['FLOP/s'].values[0] / 10**9),
        "particle-L2": (new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "particle"))]['FLOPS/Byte-miss'].values[0], 
                     new_df.loc[((new_df['Benchmark'] == "S1") & (new_df['Region'] == "particle"))]['FLOP/s'].values[0] / 10**9),
        }
plot_roofline("Hemocell AMD Genoa", {"L1": 38798, "L2": 5055, "LLC": 1268, "DRAM": 700}, {"peak": 4852}, kernels)

plt.savefig(snakemake.output["S1_png"], dpi=300)
plt.savefig(snakemake.output["S1_pdf"])

# Code to plot the rooline
fig, ax = plt.subplots(1, figsize=figsize, layout='constrained', sharex=True)


kernels = {
        "fluid-L1": (new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "fluid"))]['FLOPS/Byte'].values[0], 
                  new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "fluid"))]['FLOP/s'].values[0] / 10**9),
        "particle-L1": (new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "particle"))]['FLOPS/Byte'].values[0], 
                     new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "particle"))]['FLOP/s'].values[0] / 10**9),

        "fluid-L2": (new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "fluid"))]['FLOPS/Byte-miss'].values[0], 
                  new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "fluid"))]['FLOP/s'].values[0] / 10**9),
        "particle-L2": (new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "particle"))]['FLOPS/Byte-miss'].values[0], 
                     new_df.loc[((new_df['Benchmark'] == "S2") & (new_df['Region'] == "particle"))]['FLOP/s'].values[0] / 10**9),
        }
plot_roofline("Hemocell AMD Genoa", {"L1": 38798, "L2": 5055, "LLC": 1268, "DRAM": 700}, {"peak": 4852}, kernels)

plt.savefig(snakemake.output["S2_png"], dpi=300)
plt.savefig(snakemake.output["S2_pdf"])
